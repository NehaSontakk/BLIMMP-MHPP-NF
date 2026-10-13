# Porting blimmp-nf to a new cluster

This is written for an agent working with an operator on a cluster this pipeline has never run
on. The goal is a real run of real data, and a pull request carrying whatever had to change to
get there.

The pipeline runs today on two Slurm clusters, University of Arizona HPC and CU Boulder Alpine.
Everything that differs between those two is listed in [Running on HPC](hpc.md#what-differs-and-what-it-costs-to-get-wrong),
and that table is the shape of what will differ on a third.

## Ground rules

**Your operator has a login and you do not.** Anything needing interactive authentication, a
password, or a browser is theirs to run. Hand them the command and wait.

**Confirm before submitting to a shared queue.** Batch jobs bill an account and occupy nodes
that other people are waiting for. Show the script and the expected cost before the first
`sbatch`, and once the operator is comfortable with the pattern, keep going without asking each
time.

**Prefer a site config over editing the pipeline.** `conf/site_template.config` exists so a new
cluster needs no change to `main.nf` or `nextflow.config`. The `slurm` profile plus a `-c` file
is the shape that fits most clusters, though not all: a site that does not share nodes wants the
executor inside one allocation instead, and `conf/stampede3.config` is what that looks like. Reach for a code change only when you
have a concrete failure that configuration cannot express, and say what that failure was.

**Report what happened, including what failed.** A port that needed four attempts is more useful
written down than a port that claims to have worked first time.

## Phase 0: reconnaissance

Answer these before running anything. Most come from the cluster's own documentation, which
outlasts anything written here.

| Question | How to find out |
|---|---|
| Partition names | `sinfo -s` |
| Account to bill | the operator, or `sacctmgr show assoc user=$USER` |
| Is a QoS mandatory, and which | cluster docs; a wrong QoS is rejected at submission |
| How memory must be requested | cluster docs. This is the dangerous one, see below |
| Container engine | `which singularity apptainer docker`, and `module avail` for all three |
| Nextflow available | `module avail nextflow`; if none, it installs into a directory you own |
| Java version | `java -version`, and `module avail java openjdk` |
| Do compute nodes reach the internet | submit a one-line job that curls something |
| Shared scratch that compute nodes can read | cluster docs; `$HOME` is often too small or not mounted |
| Per-user job limit | cluster docs; it sets `--max_concurrent` |

**Memory is the question that punishes you silently.** Some Slurm sites honour `--mem`, some
only `--mem-per-cpu`, and sbatch rejects a job that carries both. Nextflow renders its own
`memory` directive as `--mem`. On UA HPC, `--mem` is accepted and then ignored, allocating 1 MB,
so jobs die partway through a long run with no error that points at memory. Every other row in
the table above fails loudly at submission. This one does not.

If your cluster ignores `--mem`, the fix is the `uahpc` profile's shape: clear `memory` and put
the request in `clusterOptions` instead. Read that profile in `nextflow.config` before writing
your own.

## Phase 1: setup

```bash
git clone https://github.com/NehaSontakk/BLIMMP-MHPP-NF.git
cd blimmp-nf
```

**Nextflow.** The floor is 22.10, and CI covers 22.10.7, 24.04.4, 25.10.4 and 26.04.0. Nextflow
24.11 and later require Java 17; earlier releases run on Java 11. The installer does not install
Java, which is the step that catches people out. If no module is close enough:

```bash
curl -s https://get.nextflow.io | bash
```

**Container cache.** Point it somewhere shared and persistent, not `$HOME`, and make sure
compute nodes can read it:

```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/shared/singularity_cache
```

Seven images are pulled: HMMER, pyrodigal, tantan, BUSCO, BLIMMP, `python:3.11` and
`ubuntu:22.04`. If compute nodes have no outbound network, pull them from a login node first by
running the bundled test there, or by pulling each one by hand.

**Reference data.** `fetch_reference_data.sh` pulls it from the places it comes from:

```bash
./fetch_reference_data.sh --kofam   /path/to/shared/kofam
./fetch_reference_data.sh --genomes /path/to/genomes
```

Each part skips what is already there, so an interrupted download resumes by re-running.

KOfam is 1.6 GB compressed and about 7.2 GB extracted, 27,262 profiles from
[GenomeNet](https://www.genome.jp/ftp/db/kofam/). Put it on shared storage rather than in a home
directory; one copy serves everyone on the cluster who uses this.

`--genomes` fetches the nine closed *Prochlorococcus* assemblies the pipeline was developed
against, about 17 MB, from NCBI by accession. They are useful as a known-answer set: the same
input has been run here, so a result that differs is about the cluster rather than the data.
Edit the accession list in the script for your own.

Every input is public and versioned upstream, so reproducing our tests takes no transfer from
anyone: run the same command and you get the same bytes. Nothing here depends on a person being
available, and the nine genomes are verifiable, since a run that differs from ours differs
because of the cluster.

Run the fetch wherever this cluster has outbound network, which is nearly always the login node,
and write it to storage the compute nodes can read. They never need network of their own: the
containers come from the cache you set with `NXF_SINGULARITY_CACHEDIR`, and a BUSCO lineage
comes from `--busco_downloads`. If even the login node is walled off, your site will have a
documented route for pulling public reference data, and this is ordinary public reference data.

**A BUSCO lineage, only if you need one.** Read
[how complete are your genomes](overview.md#how-complete-are-your-genomes) first: for closed
assemblies `--sigma 1.0` is the right answer and BUSCO is the wrong one. For drafts and MAGs:

```bash
busco --list-datasets
./fetch_reference_data.sh --busco bacteria_odb12 /path/to/busco_downloads
```

**Scaling the test down.** You do not need a smaller dataset to do a smaller run.
`--subsample_hmms 2000` and `--subsample_genomes 2` cut the work without changing anything about
the inputs, which is how to exercise chunking and the scheduler before committing a full run.

## Phase 2: the test ladder

Climb it in order. Each rung isolates a different thing, and skipping to the top makes a failure
much harder to attribute.

**1. Stub run. No containers, no data, no queue.**

```bash
nextflow run . -profile test -stub-run
```

This parses `main.nf` and walks the whole workflow with empty files. It catches a Nextflow
version incompatibility and nothing else. It runs on a login node.

**2. Bundled dataset, locally.**

```bash
nextflow run . -profile test,singularity
assets/test/assert.sh results_test
```

Two 150 kb genome slices against 20 KOfam profiles, about 30 seconds. This proves the containers
pull and run. Run it wherever containers work, which on most clusters means a compute node.

`assert.sh` matters more than the exit code. BLIMMP 0.1.2 once shipped without its KEGG data and
produced structurally perfect, entirely empty output for four months with every task green. Every
check in that script asserts a number.

**3. Bundled dataset, through the scheduler.**

```bash
cp conf/site_template.config conf/<cluster>.config
# fill in queue, clusterOptions, and memory as your cluster spells it
nextflow run . -profile test,slurm,singularity -c conf/<cluster>.config
assets/test/assert.sh results_test
```

Same 30 seconds of work, now submitted as jobs. This is where partition names, accounts, QoS and
memory spelling fail. Expect to iterate here, and expect the failures to be uninformative:
`sacct -j <id> --format=JobID,State,ExitCode,Reason` and the `.command.err` in the task's work
directory are where the real message is.

**4. A handful of real genomes.**

Copy `run_test.sh`, fill in the placeholders, submit it. Ten or so of the operator's own genomes
against the real KOfam. This is the first rung that exercises chunking, the queue at depth, and
storage under load.

**5. A real analysis.**

`run_analysis.sh` is the template for this: offline BUSCO, a taxonomy, a pinned resume session.
Read its comments before running it, particularly the one about pinning `-resume` to a session
id rather than using it bare.

## The known-answer run

Rung 4 has a right answer. The nine genomes `fetch_reference_data.sh` downloads have been run on
two clusters, University of Arizona HPC and CU Boulder Alpine, and both produced these counts
exactly:

```
nextflow run main.nf \
    --genomes '<genomes>/*/*.fna' \
    --hmms <kofam> \
    --sigma 1.0 \
    --max_concurrent 50 \
    --outdir results_compare \
    -profile slurm,singularity -c conf/<cluster>.config
```

| genome | count<1 | count<1e-6 |
|---|---|---|
| MED4 | 36736 | 22538 |
| MIT0801 | 39609 | 23639 |
| MIT9211 | 35621 | 21497 |
| MIT9301 | 36960 | 22446 |
| MIT9303 | 52400 | 31265 |
| MIT9312 | 37329 | 22884 |
| MIT9313 | 48954 | 29769 |
| NATL2A | 38339 | 23203 |
| SS120 | 37175 | 22226 |

Four things have to match for that table to mean anything, and the first is the one that bites:

- **KOfam release 2026-07-02, 27,816 profiles.** This is the default in `fetch_reference_data.sh`.
  Count them before comparing: `find <kofam> -maxdepth 1 -name '*.hmm' | wc -l`. The rolling
  `profiles.tar.gz` at the top of the KOfam tree gave 27,262, 27,816 and 27,864 over a few weeks
  in 2026, and a KO count from a different release tells you nothing about your cluster.
- `--sigma 1.0` and no `--busco_lineage`. These are closed assemblies.
- Masking left at its default, which is on. `--no_masking` changes module calls.
- BLIMMP_MHPP `26sep2026-mhpp1` and HMMER 3.4, which is what the containers pin.

Compare sorted, because chunk boundaries depend on `--max_concurrent` and row order differs
legitimately between clusters:

```bash
diff <(sort results_compare/summary.tsv) <(sort docs/reference_summary.tsv)
```

Not mismatches: wall time, CPU hours, task count, chunk count, row order. Those are scheduling.
A count that differs on the same pinned release is real, and worth chasing before you trust
anything else the cluster produces.

## Phase 3: known landmines

Ordered by how much time they cost when hit.

| Symptom | Cause |
|---|---|
| Jobs killed partway through a long run, no memory error | `--mem` ignored by the scheduler. See Phase 0 |
| `--mem, --mem-per-cpu and --mem-per-gpu are mutually exclusive` | The profile sets `memory` and the site config also sets `--mem-per-cpu`. Clear `memory` |
| A resume recomputes everything | Bare `-resume` picked up a different run. Pin the session id |
| `command not found` inside one task, hours in | A run naming no container profile. The preflight check catches this at startup; do not bypass it |
| BUSCO fails or hangs on every genome | Compute nodes have no outbound network. Use `--busco_downloads` |
| Container cannot see a path you passed | `autoMounts` covers `$HOME`, `/tmp` and the work directory. Anything else needs `-B` in `singularity.runOptions` |
| Work directory vanishes between tasks | `/tmp` is node-local on most clusters. Put `-w` on shared storage |
| Thousands of tiny jobs, each waiting far longer than it runs | Chunks too small. The pipeline autosizes from `--max_concurrent`; check what it logged |
| Nextflow refuses to start, complains about Java | A `module load` that silently did nothing. Check `java -version` after loading, not before |

## Phase 4: contributing the fixes back

Branch off `main`, which is the default branch and where the work is:

```bash
git checkout main && git pull
git checkout -b port-<cluster>
```

What a port usually contributes, in descending order of value:

**A site config.** `conf/<cluster>.config`, following `conf/site_template.config`. Keep the
comments explaining why each value is what it is; the next person needs the reasoning more than
the value.

**A row in the cluster comparison table.** `docs/hpc.md` sets UA HPC and Alpine side by side
under "What differs, and what it costs to get wrong". A third column is the single most useful
thing a port can add, because it turns two data points into a pattern. Cite the cluster's own
documentation pages, which will outlast the table.

**A worked example.** `docs/hpc.md` carries a complete sbatch script per cluster. Add yours if it
differs structurally rather than only in values.

**A profile in `nextflow.config`.** Only if the cluster needs something a `-c` file cannot
express, and say what that was. `slurm` plus a site config is the intended path.

**A trap, if you found one.** Anything that failed silently or misleadingly belongs in
`docs/reporting-problems.md` or as a comment where the code does the surprising thing.

### Before opening the PR

```bash
python3 .github/scripts/param_table.py --check   # after any nextflow_schema.json edit, run --write first
python3 .github/scripts/check_links.py
nextflow run . -profile test -stub-run
nextflow run . -profile test,singularity && assets/test/assert.sh results_test
```

CI runs six jobs on every push and pull request: a docs check, the bundled dataset under
Nextflow 22.10.7, 24.04.4, 25.10.4 and 26.04.0 with Docker, and once more under Apptainer. A
change that only touches `conf/` and `docs/` will not break them, but check anyway.

Commit messages in this repo are prose, not bullet lists. They say what changed, what it fixes,
and what was measured rather than assumed. `git log` shows the register. State which cluster,
which Slurm version, and what you actually ran.

## What to report back to the operator

- Which rung of the ladder each attempt reached, and the error that stopped it.
- The Phase 0 table filled in for this cluster.
- Anything that failed silently, meaning it produced no error and the wrong result.
- Wall clock and queue time for the real run, so the next person can size an allocation.
- What you changed, and whether it is in the PR.
