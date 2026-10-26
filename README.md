# BLIMMP-MHPP-NF

A Nextflow pipeline for identifying KEGG metabolic modules in bacterial genomes with
**BLIMMP_MHPP** — the MHPP001 (`K99001`/`pme`) module addition to
[BLIMMP](https://github.com/TravisWheelerLab/BLIMMP) (Bayesian Likelihood Inference of Metabolic
Module Presence), maintained at
[NehaSontakk/BLIMMP_MHPP](https://github.com/NehaSontakk/BLIMMP_MHPP).

This is a fork of [TravisWheelerLab/blimmp-nf](https://github.com/TravisWheelerLab/blimmp-nf)
(BSD-3-Clause, © 2026 Travis Wheeler — see [LICENSE](LICENSE)), rewired to run the MHPP fork of
BLIMMP instead of upstream. Everything about the Nextflow workflow itself — gene calling,
masking, HMMER search, chunk sizing, BUSCO completeness — is unchanged from upstream; only the
`RUN_BLIMMP` step and its container differ. See
[docs/developer.md § The BLIMMP container](docs/developer.md#the-blimmp-container) for exactly
what's different and why.

Give it genome assemblies. It calls genes with pyrodigal, masks low-complexity regions in the
predicted proteins with tantan, searches them against the KOfam HMM database (plus `K99001` for
the MHPP module), optionally measures genome completeness with BUSCO, and runs BLIMMP_MHPP to
infer the probability that each KEGG module — including MHPP001 — is present, accounting for
what a draft assembly is missing.

## Installation

1. **Install Nextflow** (22.10 or later; needs Java 11, or Java 17 for Nextflow 24.11+):

   ```bash
   curl -s https://get.nextflow.io | bash
   mkdir -p ~/bin && mv nextflow ~/bin/   # make sure ~/bin is on your PATH
   nextflow -version
   ```

   On macOS, `brew install nextflow` pulls in a working Java automatically. On a cluster, try
   `module load nextflow` first. Full details — the Java-version pairing, installing Java with
   SDKMAN, running a specific Nextflow version, cluster module systems — are in
   [docs/configuration.md § Getting Nextflow](docs/configuration.md#getting-nextflow).

2. **Install Docker**, or **Singularity/Apptainer** on HPC. Every tool this pipeline uses —
   pyrodigal, tantan, HMMER, BUSCO, and BLIMMP_MHPP — runs in a public container image that's
   pulled automatically, so this is the only other thing to install:
   - Docker Desktop: <https://www.docker.com/products/docker-desktop/>
   - Singularity/Apptainer: <https://apptainer.org/>

   **On HPC, don't actually install this yourself** — no root access, and the cluster almost
   certainly already has it. `module load singularity` (or `apptainer`) failing with "unknown
   module" doesn't mean it's missing; check whether it's just on `PATH` already instead:

   ```bash
   which singularity apptainer
   module avail 2>&1 | grep -i -E "singular|apptain"
   ```

   Many clusters ship it as a system-wide binary (no module needed) rather than, or in addition
   to, a module — worth checking before assuming it needs installing.

3. **Clone this repository**:

   ```bash
   git clone https://github.com/NehaSontakk/BLIMMP-MHPP-NF.git
   cd BLIMMP-MHPP-NF
   ```

Nothing else to build or configure — continue to **Check that it works** below.

## Check that it works, in about 30 seconds

First, check which one you actually have — don't guess:

```bash
docker --version          # prints a version if you have Docker
singularity --version     # prints a version if you have Singularity
apptainer --version       # prints a version if you have Apptainer
```

Whichever one printed a version (not "command not found") is the one you have. Run the matching
block below exactly as written — nothing in it needs editing, and nothing about this step needs a
scheduler or a specific cluster, even if you're on one:

If `docker --version` worked:

```bash
nextflow run main.nf -profile test,docker
assets/test/assert.sh results_test
```

If `singularity --version` or `apptainer --version` worked:

```bash
nextflow run main.nf -profile test,singularity
assets/test/assert.sh results_test
```

**On a shared or HPC machine, a plain interactive shell may only give you 1 CPU**, which is below
what one of the steps asks for (2). If you see `Process requirement exceeds available CPUs`, it's
this, not a real failure — request more CPUs for your interactive session (the exact command
varies by cluster; on UA HPC it's `idev -c 2`) and rerun the same command with `-resume` appended,
which picks up from the step that failed instead of starting over.

Nothing to build first: `ghcr.io/nehasontakk/blimmp_mhpp` is a public image, and Docker/Singularity
pulls it automatically the first time `RUN_BLIMMP` needs it, same as the public images every other
step in this pipeline uses.

Bundled data: two 150 kb slices of real *Prochlorococcus* genomes against 20 KOfam profiles.
`assert.sh` then checks that the run found the KOs it should have, because a pipeline can exit 0
having produced nothing at all.

## Quick start

Pick a profile: where the work runs, then where the tools come from. Nothing to copy or edit.

| your machine | `-profile` |
|---|---|
| Laptop or workstation | `laptop,docker` |
| Large shared server | `server,docker` |
| Slurm cluster | `slurm,singularity`, plus a [small site config](docs/configuration.md#another-slurm-cluster) naming your account and partition |
| University of Arizona HPC | `uahpc,singularity` |
| TACC Stampede3 | `singularity` with [`conf/stampede3.config`](conf/stampede3.config), and no `slurm`: jobs get whole nodes, so Nextflow schedules inside one allocation |
| Tools already on `PATH` | just `laptop` |

Three things in the command below need replacing with your own values; everything else — every
flag name, `bacteria_odb12`, the whole shape of the command — stays exactly as written:

```bash
nextflow run main.nf \
    --genomes <path to your genome files> \
    --hmms <path to your KOfam profiles> \
    --busco_lineage bacteria_odb12 \
    -profile <profile from the table above>
```

For example, if your genomes are in `~/genomes/`, your KOfam profiles in `~/kofam/`, and you're
on a laptop with Docker, the command you'd actually type is:

```bash
nextflow run main.nf \
    --genomes ~/genomes/ \
    --hmms ~/kofam/ \
    --busco_lineage bacteria_odb12 \
    -profile laptop,docker
```

`--genomes` takes a directory of `.fna`/`.fa`/`.fasta` files or a glob (point it at one file, or
a directory with exactly one, to run a single genome); `--proteins` takes protein FASTAs instead
and skips gene calling. `--hmms` takes one concatenated `.hmm` file or a directory of them.
`--busco_lineage` measures how complete each sample is, which is what lets BLIMMP tell an absent
pathway from an unsequenced one — it works on `--proteins` too (BUSCO's `proteins` mode, instead
of `genome` mode), not just `--genomes`; `--busco_lineage auto_prok` lets BUSCO choose the
lineage for each sample, and `--sigma 1.0` states the fraction yourself instead (use this for
closed, single-contig assemblies — there's nothing for BUSCO to measure).

Results land in `results/`, one subdirectory per genome, holding the HMMER tables and
BLIMMP_MHPP's per-module probabilities (`*_BLIMMP_substituted_*` files — see
[docs/running.md](docs/running.md) for the full output layout). `nextflow run main.nf --help`
lists every parameter.

### Testing the MHPP module specifically

The MHPP001 module needs its marker KO, `K99001` (placeholder for the `pme` enzyme), which is
**not yet packaged anywhere versioned** — you have to supply `K99001.hmm` yourself and drop it
into whatever `--hmms` directory you point at, alongside the regular KOfam profiles. Three things
to replace: where you want the KOfam profiles downloaded to, where your own copy of `K99001.hmm`
currently is, and your genome file:

```bash
./fetch_reference_data.sh --kofam <where to put the KOfam profiles>   # full KOfam, ~7.2 GB extracted
cp <path to your copy of K99001.hmm> <same KOfam path as above>

nextflow run main.nf \
    --genomes <path to one genome file> \
    --hmms <same KOfam path as above> \
    --sigma 1.0 \
    -profile <profile from the table above>
```

Worked example — KOfam going to `~/kofam/`, your own `K99001.hmm` sitting in `~/Downloads/`, one
genome at `~/genomes/strain1.fna`, laptop with Docker:

```bash
./fetch_reference_data.sh --kofam ~/kofam/
cp ~/Downloads/K99001.hmm ~/kofam/

nextflow run main.nf \
    --genomes ~/genomes/strain1.fna \
    --hmms ~/kofam/ \
    --sigma 1.0 \
    -profile laptop,docker
```

Then check the module's confidence (replace `<sample>` with your genome's name, which matches the
filename you gave `--genomes` minus its extension):

```bash
grep MHPP001 results/<sample>/<sample>_BLIMMP_substituted_module_probabilities.csv
```

**Plain Docker needs an extra bind mount for a directory-mode `--hmms`.** `GROUP_HMMS` receives a
directory path as a raw string rather than as a staged Nextflow `path` input (staging ~28k
individual symlinks isn't practical), so Docker never automatically sees it — unlike single
staged files, which Nextflow auto-binds by resolving their symlink targets. Singularity's
`autoMounts` doesn't cover this either (see the "needs a bind of your own" note in
[docs/hpc.md § BUSCO lineages, and why to download one first](docs/hpc.md#busco-lineages-and-why-to-download-one-first),
which covers the same gap for `--busco_downloads`). On plain Docker Desktop, add a one-off config
— replace both copies of the path with your actual `--hmms` path (it must be the same path on
both sides of the colon):

```groovy
// local_mount.config
docker {
    runOptions = '-u $(id -u):$(id -g) -v <your --hmms path>:<your --hmms path>'
}
```

Then add `-c local_mount.config` to whichever `nextflow run` command you were already running —
nothing else about it changes:

```bash
nextflow run main.nf --genomes ... --hmms ... -profile laptop,docker -c local_mount.config
```

## Documentation

- [What the pipeline does](docs/overview.md): the steps, completeness, taxonomic priors, and re-scoring without re-searching
- [Configuration](docs/configuration.md): installing Nextflow, and the profile for your machine
- [Running the workflow](docs/running.md): inputs, output layout, resuming, dry runs
- [Running at scale](docs/hpc.md): submitting to a scheduler, chunk size, concurrency, cleanup
- [Parameters](docs/parameters.md): every flag, generated from the schema
- [Samplesheets](docs/samplesheets.md): a completeness and taxonomy per genome, with ten worked examples
- [When something goes wrong](docs/reporting-problems.md): diagnosing a failure, and what to include in a report
- [Porting to a new cluster](docs/porting.md): reconnaissance, a test ladder, and how to send the fixes back
- [Developer guide](docs/developer.md): architecture, channel patterns, adding a step, **and how the BLIMMP_MHPP container is built and released**

## Requirements

- [Nextflow](https://www.nextflow.io/) 22.10 or later. Nextflow 24.11 and later need Java 17;
  earlier releases run on Java 11. The Nextflow installer does not install Java, which is the
  step that catches people out. See [getting Nextflow](docs/configuration.md#getting-nextflow).
  (`brew install nextflow` on macOS pulls in a working Java automatically.)
- [Docker](https://www.docker.com/), or [Singularity/Apptainer](https://apptainer.org/) on HPC.
  HMMER, pyrodigal, tantan, BUSCO and BLIMMP_MHPP all run in public containers — nothing to build.
  Without a container profile you'd need HMMER 3.x, pyrodigal 3.x, tantan, BUSCO 6.x and
  BLIMMP_MHPP installed yourself; BUSCO and BLIMMP_MHPP have no supported non-container path here.
- The KOfam HMM profiles, 27,816 of them (release `2026-07-02`), about 7.2 GB extracted.
  `./fetch_reference_data.sh --kofam /path/to/kofam` downloads and unpacks them from
  [GenomeNet](https://www.genome.jp/ftp/db/kofam/). The same script fetches example genomes and a
  BUSCO lineage; `--help` lists what it does.
- `K99001.hmm`, to detect the MHPP001 module specifically — see
  [Testing the MHPP module specifically](#testing-the-mhpp-module-specifically) above. Not
  bundled anywhere yet; you'll need a copy from wherever BLIMMP_MHPP work currently keeps it.

## If it does not work on your system

Nextflow's reach is the reason this pipeline exists in this form: the same workflow runs on a
laptop, a Slurm cluster, or AWS Batch. That same reach is what makes it hard to test. Upstream has
covered the combinations it can get at, and HPC environments differ in ways that can't be
predicted from here. How memory must be requested, whether compute nodes can reach the network,
which module names exist, how Singularity was built: all of it varies.

If something does not work, open an issue on this repo with the details of your system and the
error; [when something goes wrong](docs/reporting-problems.md) lists what to include and gives
you a command that collects most of it.

## License

See [LICENSE](LICENSE) — BSD-3-Clause, © 2026 Travis Wheeler, carried over unmodified from
upstream `blimmp-nf`. The BLIMMP_MHPP source this pipeline runs is separately licensed in its own
repository.
