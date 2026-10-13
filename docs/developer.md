# Developer guide

This guide is for contributors adding new processing steps to blimmp-nf. It covers how the workflow is structured, the Nextflow DSL2 patterns it uses, and two concrete examples of appending downstream analysis steps.

---

## Repository layout

```
blimmp-nf/
├── main.nf              # All process definitions and the workflow block
├── nextflow.config      # Parameters, resource labels, execution profiles
├── nextflow_schema.json # Every parameter, its type and its help text
├── bin/
│   ├── split_hmms.py    # SPLIT_HMMS: one .hmm file into chunks
│   ├── group_hmms.py    # GROUP_HMMS: a directory of .hmm files into chunks
│   ├── postprocess.py   # POSTPROCESS: hit counts per genome
│   └── module_matrix.py # MODULE_MATRIX: module confidences, genomes as columns
├── conf/
│   ├── test.config      # The bundled dataset, behind -profile test
│   ├── alpine.config    # CU Boulder Alpine
│   ├── stampede3.config # TACC Stampede3, which allocates whole nodes
│   └── site_template.config # Starting point for a cluster with no profile
├── assets/
│   ├── schema_input.json # The samplesheet columns, in the same format
│   └── test/            # The bundled dataset behind -profile test, and assert.sh
├── .github/
│   ├── workflows/       # CI
│   ├── scripts/         # param_table.py, check_links.py, run by the Docs job
│   └── ISSUE_TEMPLATE/  # what a bug report should carry
└── docs/                # Documentation
```

Scripts placed in `bin/` are automatically added to `PATH` inside every process at runtime — Nextflow handles this for both local and container execution. Add new helper scripts there rather than referencing them via `${projectDir}`.

---

## How the workflow is structured

### The two axes of parallelism

The core design is a cross-product over two independent dimensions:

```
genome_1 ──┐
genome_2 ──┤   × [chunk_001, chunk_002, ..., chunk_N]  →  HMMSEARCH (G × N jobs)
...        ─┘
```

- **Genome axis**: `PREDICT_PROTEINS` runs once per genome, producing a `.faa` protein file.
- **HMM axis**: the HMM database is split into N chunks by `GROUP_HMMS` or `SPLIT_HMMS`.
- `HMMSEARCH` runs for every (genome, chunk) pair — `G × N` jobs in total.
- `MERGE_RESULTS` collects all N chunk results for each genome back into a single `.tblout`/`.domtblout` — one job per genome, running as soon as all its chunks finish.

`N` is chosen by `autoChunkCount()` unless `--hmm_chunks` or `--hmm_chunk_size` is given. It reads the executor and queue size out of `workflow.session.config.navigate()` rather than inferring them from profile names, so it behaves the same for a site config it has never seen: on a grid executor a slot is a queue slot, and on a local one it is a core divided by what `process_medium` asks for in whatever profile is in force. The target is `slots / genomes`, enough tasks to fill the slots once. `minProfilesPerChunk()` in `main.nf` mirrors `MIN_CHUNK_SIZE` in both chunking scripts and stops a run with many slots and one genome from producing chunks too small to be worth scheduling; the two constants have to move together, since the first is what the log advertises and the second is what actually happens.
- `POSTPROCESS` runs once, after all `MERGE_RESULTS` jobs complete.

### The BUSCO branch

`RUN_BUSCO` reads `genomes_ch` directly, so it does not wait for pyrodigal or for anything in the HMM search. It runs alongside that whole branch, drawing from the same pool of scheduler jobs, so BUSCO and `HMMSEARCH` are submitting against each other and `--max_concurrent` covers both. The BUSCO branch requires `--busco_lineage`. `RUN_BLIMMP` no longer does: it takes sigma as a value, which `--busco_lineage` supplies per genome by reading it out of BUSCO's TSV in the workflow, and `--sigma` supplies as one number for every input. The two are mutually exclusive at validation, and with neither the workflow skips BLIMMP and logs why.

`BUSCO_SUMMARY` collects every per-genome completeness value into one `busco_summary.tsv`. `RUN_BLIMMP` instead joins the per-genome completeness against the merged `.domtblout` on `sample_id`, giving one BLIMMP job per genome.

### Why the BLIMMP check is its own process

`CHECK_BLIMMP` does nothing but count rows in `dk.csv`, which looks like it belongs at the end of `RUN_BLIMMP`. It is separate because Nextflow publishes nothing at all from a task that exits non-zero. Fail inside `RUN_BLIMMP` and the debug log, the one artifact that explains the failure, stays in the work directory. Splitting them means BLIMMP's own exit status and our verdict on its output are two different tasks: `RUN_BLIMMP` succeeds and publishes, `CHECK_BLIMMP` stops the run.

`RUN_BLIMMP` also installs an `EXIT` trap that copies the debug log to `--outdir` as a failing task exits, which covers a BLIMMP crash rather than a rejected result. The trap needs `--outdir` reachable from inside the container. That held under `autoMounts` with `--cleanenv` on a compute node here, and will not hold everywhere, so a failed copy reports where the log is and leaves the underlying error to propagate.

It costs one thing worth knowing. The trap names an absolute destination, which puts `--outdir` into the command string and therefore into the task hash, so changing `--outdir` re-runs the BLIMMP tasks where it would normally only republish. Measured: 9 `RUN_BLIMMP` tasks re-run, all 2457 `HMMSEARCH` tasks still cached. Optional flags elsewhere in this process carry their own leading space (`' -l'`, `" -t ${...}"`) precisely so that leaving one unset does not perturb the command string for everyone else.

### Entry points

Two workflows live in `main.nf`. The unnamed `workflow { }` is what runs by default and covers the whole pipeline. `workflow blimmp_only` starts from finished `.domtblout` files and calls only `RUN_BLIMMP` and `CHECK_BLIMMP`, reached with `-entry blimmp_only`.

The pipeline was not carved into sub-workflows to make this work. Adding a named workflow beside the default one is enough, and the processes already take what the second entry point can supply: a domtblout, a number, and an optional taxonomy name. That is what decoupling sigma from BUSCO bought.

Three functions are shared between the two, so a rule cannot hold in one entry point and not the other:

| Function | What it does |
|---|---|
| `checkSigma(value, source)` | Range-checks one completeness fraction. `source` names where it came from, so the message points at the flag, the sheet row, or the completeness table |
| `validatedSigma()` | Checks `--sigma` and returns it as a number, or null if unset |
| `expandInputs(spec, extensions)` | Turns a directory or a glob into a file list |
| `readCompleteness(path)` | Reads a sample-ID-to-completeness table, skipping any row whose second column is not a number, which covers the `busco_summary.tsv` header and accepts headerless per-genome TSVs |

`blimmp_only` refuses `--genomes`, `--proteins`, `--hmms`, and `--busco_lineage` rather than ignoring them, and refuses `--sigma` together with `--completeness`. It also checks up front that every sample has a row in the completeness table, since these inputs are plain lists at that point rather than channels, so a missing value can be reported by name instead of surfacing later as an unset sigma.

Equivalence is the property worth protecting when changing either path: a full run with `--busco_lineage` and a `blimmp_only` run over its outputs produce byte-identical BLIMMP files, verified across 36 files on 9 genomes.

### Three ways in, two channels out

`--genomes`, `--proteins` and `--input` all end up feeding the same two channels: `genomes_ch` for nucleotides that need gene calling, and `ready_proteins_ch` for proteins that do not. `PREDICT_PROTEINS` has one call site, fed by a channel that is simply empty on a run that starts from proteins.

`readSamplesheet(path)` parses the sheet into one map per row and does every check there, before a task is submitted, so a bad row is reported by name and line number rather than surfacing later as a missing file. `resolveInputCell()` turns one cell into a sorted file list, resolving relative paths against the sheet's own directory rather than the launch directory.

Blank `sigma` cells stay null through parsing rather than being filled from `--sigma` on the spot. Sigma has three sources with a rank order (the cell, then BUSCO's measurement, then the flag), and collapsing the cell and the flag early would let `--sigma` outrank a measurement it should lose to. The resolution happens in the workflow, where all three are in view.

**Why `CONCAT_INPUTS` is a separate process.** A samplesheet row may name several files, which have to become one sequence before gene calling. The obvious implementation is to widen `PREDICT_PROTEINS` to take a file list and concatenate inside its script. Do not do this. The output name is derived from the input file, so widening the process rewrites its command string, which changes its task hash, which produces a new `.faa`, which re-runs every `HMMSEARCH` downstream. On nine genomes that is 2,457 searches and 8.3 CPU hours, for a feature nobody using the old input flags asked for. `RUN_BUSCO` needs a single file too, and would have needed its own copy of the same logic.

Joining in a separate process keeps single-file rows on exactly the path they were on before samplesheets existed, with the hashes they already had. Only multi-file rows touch `CONCAT_INPUTS`.

**Why `MASK_REPEATS` is a separate process.** The `CONCAT_INPUTS` argument applies again. tantan rewrites the `.faa`, and the `.faa` is an input to `HMMSEARCH`, so masking inside `PREDICT_PROTEINS` would change every search hash for everyone, including runs that want no masking. As its own process, the unmasked path that `--no_masking` selects keeps hashes of its own, and neither setting disturbs the other.

It sits after `PREDICT_PROTEINS.out.mix(ready_proteins_ch)`, which is the one point where all three input routes have converged, so a single call site covers proteins pyrodigal predicted and proteins you supplied. Switching the flag either way re-searches everything, which is correct: the sequences are different, and so are the results.

It is invoked once, over both genome and protein multi-file rows together, because Nextflow allows a process one invocation per workflow. The two kinds are separated again afterwards by name, and an `extension` value in the input tuple decides what the joined file is called.

Files are sorted before joining. Glob expansion order is not guaranteed across filesystems, and concatenation order reaches the task hash.

### Channel flow

```
samplesheet ─→ readSamplesheet ─┬─→ (rows naming >1 file) ─→ CONCAT_INPUTS ─┐
                                └─────────────────────────────────────────┐ │
                                                                          ↓ ↓
genomes_ch ─┬─→ PREDICT_PROTEINS ─→ proteins_ch ─┐
            │                    (MASK_REPEATS) ─┤
            │                                    ├─→ HMMSEARCH ─→ MERGE_RESULTS ─┬─→ POSTPROCESS
            │  hmm_dir/file ─→ GROUP/SPLIT_HMMS ─┘        (flattened)            │
            │                                                                    │ (domtblout)
            └─→ RUN_BUSCO ─┬─→ BUSCO_SUMMARY                                     │
                           │                                                     │
                           └─────→ join on sample_id ←───────────────────────────┘
                                          │
                                          └─→ RUN_BLIMMP ─┬─→ (dk.csv) ─→ CHECK_BLIMMP
                                                          └─→ (probabilities) ─→ MODULE_MATRIX
```

In Nextflow DSL2, channels carry data between processes. Key operators used here:

| Operator | Where used | What it does |
|---|---|---|
| `.combine()` | `proteins_ch.combine(hmm_chunks)` | Cross product — every protein file paired with every HMM chunk |
| `.groupTuple()` | Feeding `MERGE_RESULTS` | Groups per-chunk results by `sample_id` so all chunks for one genome arrive together |
| `.join()` | Joining tblout and domtblout groups | Aligns two channels by key (`sample_id`) |
| `.collect()` | Feeding `POSTPROCESS` | Waits for all items and passes them as a list to one job |
| `.flatten()` | After `GROUP/SPLIT_HMMS` | Turns a list of chunk files into individual channel items |
| `.map()` | Various | Transforms or selects fields from a tuple |

### Resource labels

All processes are tagged with one of two labels defined in `nextflow.config`:

| Label | CPUs | Time | Used for |
|---|---|---|---|
| `process_low` | 1 | 1h | PREDICT_PROTEINS, CONCAT_INPUTS, MASK_REPEATS, GROUP/SPLIT_HMMS, MERGE_RESULTS, POSTPROCESS, MODULE_MATRIX, BUSCO_SUMMARY, CHECK_BLIMMP |
| `process_medium` | 4 | 8h | RUN_BUSCO, RUN_BLIMMP |
| `withName: HMMSEARCH` | 2 | 8h | HMMSEARCH, which takes `process_medium`'s time and its own cpus |

New processes should use one of these labels. If a new step is genuinely heavier, add a `process_high` label to every execution profile in `nextflow.config`, and to any site
config of your own.

HMMSEARCH is the exception. If a cpus setting seems to be ignored, this is why: its cpus come
from a `withName` selector, which outranks every `withLabel`, so a config that raises
`process_medium.cpus` moves BUSCO and BLIMMP and leaves the search alone. It is set by name
because hmmsearch is the one process whose thread scaling we timed, and two is where the gain
stops.

Both `time` and `memory` are written as closures multiplied by `task.attempt`, so the second
attempt at a task gets twice the first one's allocation and the third gets three times.
`errorStrategy` retries only the exit codes a signal produces (130 to 145, plus 104), which is
where the out-of-memory and walltime kills live. Anything else
is a tool rejecting its input, and it will reject it again, so the run stops with `finish` and
the message stays readable.

Memory scales only where it is set. The `uahpc` profile clears `memory` and passes
`--mem-per-cpu` through `clusterOptions` instead, because Slurm there ignores `--mem`, and a
string in `clusterOptions` cannot see `task.attempt`. Retries there get more time but the same
memory. A site config that owns memory should scale it the same way the local profiles do.

---

## Adding a new parameter — checklist

`nextflow_schema.json` describes every parameter once. `--help` prints from it and the
unknown-parameter check validates against it, so a parameter cannot be documented in one place
and checked against another.

1. **Add it to `params` in `nextflow.config`** with its default.
2. **Add it to `nextflow_schema.json`**, inside whichever group in `definitions` it belongs to.
   Give it a `type`, a one-line `description`, and a `help_text` carrying the reasoning that
   would otherwise live only in the docs. Add `format` (`file-path`, `directory-path`, `path`)
   when it names something on disk, which is what makes `--help` print `<file>` rather than
   `<string>`.
3. **Miss either one and it breaks.** A parameter in the config but not the schema is absent
   from `--help`; a parameter in the schema but not the config is never set. `run_schemacheck.sh`
   compares the two lists and reports anything in one and not the other.
4. **Regenerate the parameter table**: `.github/scripts/param_table.py --write`. It rewrites the
   generated block in `docs/parameters.md` from the schema. CI runs the same script with
   `--check` and fails if you forget, so the docs cannot drift from what the pipeline accepts.
   Anything needing more than the one-line `description` goes in prose in `docs/running.md`.

The file follows the nf-core schema format, so the `nf-schema` plugin could take over `--help`
and validation unchanged. It is read directly instead, for two reasons. The plugin is fetched
at runtime, so a cluster with no outbound network could not run the pipeline at all. And
`nf-schema` requires Nextflow 25.10 or later, against this pipeline's floor of 22.10; its
predecessor `nf-validation` is deprecated. Adopting the plugin therefore means moving the
Nextflow floor forward by three years and past the point where Java 17 became mandatory,
which is a much larger change than the roughly ninety lines it would replace.

`assets/schema_input.json` describes the samplesheet columns in the same format. Nothing reads
it today; `readSamplesheet()` in `main.nf` does its own checking. It exists so the columns have
a machine-readable definition, and it is what `nf-schema` would validate sheets against.

### Nextflow versions

Runs on 22.10.7, 24.04.4, 25.10.4 and 26.04.0. All four are in the CI matrix, so all four are
checked on every push. The floor stays at 22.10 because nothing in the pipeline needs a later
release.

Two things to know before raising it. Java 17 became mandatory in 24.11, so any floor at or
above that release asks users to change their Java as well, and the Nextflow installer does not
install Java. And 26.04 turns the strict syntax parser on by default, which rejects imports,
class declarations, `for` and `while` loops, `switch`, the spread operator, and statements at
the top level outside a declaration. `NXF_SYNTAX_PARSER=v1` opts out.

Expect two kinds of breakage on 26.04, and only one of them shows up in `nextflow lint`. The
strict parser rejects the constructs listed above and the linter reports every one. It stays
silent on the second kind: 26.04 evaluates directive strings eagerly, so
`publishDir "${params.outdir}/${sample_id}"` looks up an input variable during `--help`, before
any process runs. A directive that builds a value out of process inputs needs the closure form,
`publishDir path: { ... }`.

---

## Continuous integration

`.github/workflows/ci.yml` runs on every push to `main` and on every pull request:

| job | what it covers |
|---|---|
| Docs | the parameter table still matches the schema, and no internal link is broken |
| Nextflow 22.10.7 | the supported floor |
| Nextflow 24.04.4 | the last widely deployed release that still runs on Java 11 |
| Nextflow 25.10.4 | what `module load nextflow` gives on a current cluster |
| Nextflow 26.04.0 | the strict syntax parser |
| Apptainer | the engine used on HPC, where Docker is unavailable |

The Docs job needs no containers and finishes in seconds. The other five all run the same
30-second bundled dataset.

The four Nextflow jobs each print `--help`, confirm an unknown parameter is rejected, run
`-profile test,docker`, and then run `assets/test/assert.sh`. The Apptainer job runs
`-profile test,singularity` and the same assertion.

A run can exit 0 having produced nothing, so the workflow does not treat a green pipeline as a
passing test. `assert.sh` requires that at least half the bundled KOfam profiles hit each genome, that
all four BLIMMP outputs are non-empty, that the per-KO table has rows, and that `summary.tsv`
names both samples.

On failure a job prints the tail of `.nextflow.log` and of every non-empty `.command.err`, and
uploads both as an artifact kept for seven days. To reproduce one locally:

```bash
NXF_VER=26.04.0 nextflow run . -profile test,docker
assets/test/assert.sh results_test
```

The repository is private, so these are billed Actions minutes. The first run measured 1.7
minutes per job, about 7 minutes per push, most of it spent pulling containers. The version
matrix is the expensive part and shrinking it is a one-line edit.

---

## Adding a new process — checklist

1. **Define the process** in `main.nf` above the `workflow {}` block.
2. **Assign a label** (`process_low` or `process_medium`).
3. **Assign a container** — use an existing image if the tool is already available, or add a new `quay.io/biocontainers/` image.
4. **Declare inputs and outputs** using the appropriate channel types (`val`, `path`, `tuple`).
5. **Call it from the `workflow {}` block**, wiring it to the correct upstream channel.
6. **Add `publishDir`** if the output should land in `--outdir`. Use `"${params.outdir}/${sample_id}"` for per-genome outputs so the subdirectory structure is preserved.
7. **Add a `stub:` block** so the process takes part in `-stub-run`. It only has to create the files the `output:` block declares; a stub run never inspects their contents.
   ```groovy
   script:
   """
   my_tool ${input_file} > ${output_file}
   """

   stub:
   """
   touch ${output_file}
   """
   ```
   A process with no stub aborts the whole stub run, so this is not optional. Resist the urge to echo the real command into the stub: two copies of a command line drift, which is why the `--dry_run` branches this replaced are gone.
8. **Do not add skip logic.** Re-run behaviour is Nextflow's job: `-resume` reuses a task when its inputs, script, and container are unchanged. Hand-written file-existence checks cannot see parameter or container changes, and they treat the partial output of an interrupted run as complete.
9. **Test with subsampling** before running at scale: `--subsample_genomes 2 --subsample_hmms 50`.
10. **Test with a stub run** to verify the output directory structure: `-stub-run -profile standard`.

---

## Example 1: single-threaded post-processing step

This is the pattern to use when your analysis needs **all genomes' results at once** — for example, computing KEGG pathway completeness across the full dataset, or building a presence/absence matrix.

The key is `.collect()`, which holds the channel open until every upstream item has arrived, then emits them all as a single list to one job.

<details><summary>Full implementation</summary>
<p>

**1. Add the script to `bin/`**

Create `bin/my_analysis.py` (make it executable: `chmod +x bin/my_analysis.py`):

```python
#!/usr/bin/env python3
"""
my_analysis.py — example post-processing script.
Reads all .tblout files and writes a combined result.

Usage: my_analysis.py <file1.tblout> [file2.tblout ...] > output.tsv
"""
import sys

def main(paths):
    print("genome\thits_1e-6")
    for path in sorted(paths):
        name = path.removesuffix('.tblout')
        count = sum(
            1 for line in open(path)
            if not line.startswith('#') and float(line.split()[4]) < 1e-6
        )
        print(f"{name}\t{count}")

if __name__ == '__main__':
    main(sys.argv[1:])
```

**2. Define the process in `main.nf`**

Add this after `POSTPROCESS` and before the `workflow {}` block:

```groovy
process MY_ANALYSIS {
    label 'process_low'
    container 'python:3.11'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path tblouts   // all per-genome .tblout files, collected into one job

    output:
    path 'my_analysis.tsv'

    script:
    """
    my_analysis.py *.tblout > my_analysis.tsv
    """
}
```

**3. Wire it into the `workflow {}` block**

After the existing `POSTPROCESS` call, add:

```groovy
// collect all merged tblouts — same pattern as POSTPROCESS
all_tblouts = MERGE_RESULTS.out
    .map { sample_id, tblout, domtblout -> tblout }
    .collect()

MY_ANALYSIS(all_tblouts)
```

If `POSTPROCESS` already computed `all_tblouts`, reuse the same channel — don't recompute it:

```groovy
// define once, pass to both
all_tblouts = MERGE_RESULTS.out
    .map { sample_id, tblout, domtblout -> tblout }
    .collect()

POSTPROCESS(all_tblouts)
MY_ANALYSIS(all_tblouts)
```

Nextflow channels can be consumed by multiple processes without re-triggering upstream work.

</p>
</details>

---

## Example 2: fan-out per-genome step (parallel)

This is the pattern to use when each genome can be **processed independently** as soon as its results are ready — for example, filtering hits, annotating each genome, or computing per-genome statistics. These jobs launch immediately after each `MERGE_RESULTS` job completes, without waiting for other genomes.

The key difference from Example 1: no `.collect()`. The channel carries one `(sample_id, tblout, domtblout)` tuple per genome, so one job runs per genome.

<details><summary>Full implementation</summary>
<p>

**1. Add the script to `bin/`**

Create `bin/filter_hits.py`:

```python
#!/usr/bin/env python3
"""
filter_hits.py — filter a tblout file to hits above an E-value threshold.

Usage: filter_hits.py --evalue 1e-6 <input.tblout> > output.tblout
"""
import sys, argparse

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--evalue', type=float, default=1e-6)
    p.add_argument('tblout')
    args = p.parse_args()

    with open(args.tblout) as fh:
        for line in fh:
            if line.startswith('#'):
                print(line, end='')
            elif float(line.split()[4]) < args.evalue:
                print(line, end='')

if __name__ == '__main__':
    main()
```

**2. Define the per-genome process**

```groovy
process FILTER_HITS {
    tag "${sample_id}"
    label 'process_low'
    container 'python:3.11'
    publishDir "${params.outdir}/filtered", mode: 'copy'

    input:
    tuple val(sample_id), path(tblout), path(domtblout)

    output:
    tuple val(sample_id), path("${sample_id}.filtered.tblout")

    script:
    """
    filter_hits.py --evalue 1e-6 ${tblout} > ${sample_id}.filtered.tblout
    """
}
```

**3. Optionally, add a downstream single-threaded summary**

Once all per-genome jobs finish, you can collect their outputs for a final summary — combining the fan-out pattern with Example 1's collect pattern:

```groovy
process FILTERED_SUMMARY {
    label 'process_low'
    container 'python:3.11'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path filtered_tblouts   // all filtered files, collected

    output:
    path 'filtered_summary.tsv'

    script:
    """
    my_analysis.py *.filtered.tblout > filtered_summary.tsv
    """
}
```

**4. Wire both into the `workflow {}` block**

```groovy
// Fan-out: one FILTER_HITS job per genome, starts as soon as that genome's
// MERGE_RESULTS completes — no waiting for other genomes
FILTER_HITS(MERGE_RESULTS.out)

// Collect all filtered outputs, then run one summary job
all_filtered = FILTER_HITS.out
    .map { sample_id, filtered -> filtered }
    .collect()

FILTERED_SUMMARY(all_filtered)
```

The resulting execution graph looks like:

```
MERGE_RESULTS (genome_1) ──→ FILTER_HITS (genome_1) ──┐
MERGE_RESULTS (genome_2) ──→ FILTER_HITS (genome_2) ──┤
...                                                    ├──→ FILTERED_SUMMARY
MERGE_RESULTS (genome_N) ──→ FILTER_HITS (genome_N) ──┘
```

`FILTERED_SUMMARY` starts only after all `FILTER_HITS` jobs complete. `FILTER_HITS` jobs start as soon as each individual genome's `MERGE_RESULTS` finishes.

</p>
</details>

---

## The BLIMMP container

`RUN_BLIMMP` uses a custom container built from the [BLIMMP_MHPP repository](https://github.com/NehaSontakk/BLIMMP_MHPP) — the MHPP001 fork of BLIMMP, not the upstream [TravisWheelerLab/BLIMMP](https://github.com/TravisWheelerLab/BLIMMP):

```
ghcr.io/nehasontakk/blimmp_mhpp:26sep2026-mhpp1
```

Unlike upstream BLIMMP, this fork is not pip-packaged: `blimmp_src/` is copied into the image at `/opt/blimmp/`, and `blimmp_container_wrapper.py` loads `BLIMMP_Version_26Sep2026.py` directly via `importlib` and wires up its `Paths()` relative to its own location. The `BLIMMP` command on `PATH` is a one-line shim that execs the wrapper. There is no pip distribution, so there's no `importlib.metadata` version to read — the Dockerfile stamps one into `/opt/blimmp/VERSION` at build time instead, which `RUN_BLIMMP` reads directly.

The cooccurrence lookup is too large for git and ships as a GitHub Release asset on that repo (tag `v0.1-mhpp-data`, asset `BLIMMP_COOCCURRENCE_MHPP.tar.gz`): the Dockerfile downloads and extracts it over a placeholder directory of the same name. The build fails, on purpose, if any required data path is missing afterward — except `ONE_HOP_NEIGHBOR_DATA-2`/`TWO_HOP_NEIGHBOR_DATA-2`, which this MHPP addition deliberately omits (its scoping is domain-level only; BLIMMP only ever `.exists()`-checks those two paths to decide on a lineage-level fallback, never reads their contents, so pointing them at nothing just forces that fallback).

Both the CLI and the output files differ from upstream BLIMMP: BLIMMP_MHPP runs a substitution pass by default and writes only `*_BLIMMP_substituted_*` files from it; the plain pass-1 files and the debug log only appear with `-v`/`--verbose` (which then writes both). `RUN_BLIMMP` always reads the substituted files, and `no_blimmp_log` now toggles `-v` rather than a log-only flag.

### Releasing a new BLIMMP_MHPP version

There is no CI workflow for this fork yet (unlike upstream BLIMMP's `.github/workflows/docker.yml`) — building and pushing the image is a manual step:

1. **Make changes** in the BLIMMP_MHPP repo on a feature branch and open a PR against `main`.

2. **Merge to main** once the PR is reviewed.

3. **Publish the data release**, if `blimmp_src/BLIMMP_COOCCURRENCE_MHPP` changed — tag the tarball as a release asset under `v0.1-mhpp-data` (or a new tag; update the Dockerfile's `BLIMMP_MHPP_DATA_RELEASE` build arg to match).

4. **Build and push the container**:
   ```bash
   cd BLIMMP_MHPP
   docker build -t ghcr.io/nehasontakk/blimmp_mhpp:<version> .
   docker push ghcr.io/nehasontakk/blimmp_mhpp:<version>
   ```

5. **Update the container tag in blimmp-nf** — edit the `container` directive in the `RUN_BLIMMP` process in `main.nf`:
   ```groovy
   container 'ghcr.io/nehasontakk/blimmp_mhpp:<version>'
   ```
   Commit and push the change to blimmp-nf.

6. **Test** with a small run before deploying to production:
   ```bash
   nextflow run main.nf \
       --genomes /path/to/genomes/ \
       --hmms /path/to/kofam/ \
       --busco_lineage bacteria_odb12 \
       --subsample_genomes 2 \
       --outdir test_results/ \
       -profile uahpc,singularity
   ```

   Or, faster and with an assertion attached, `nextflow run . -profile test,singularity`
   followed by `assets/test/assert.sh results_test`. The bundled dataset exercises BLIMMP
   end to end, so a container that lost its KEGG data files fails there in 30 seconds.

> **Note:** The container tag in `main.nf` is pinned to an exact version (not `latest`) to keep runs reproducible. Always update the pin explicitly rather than relying on `latest`.

---

## Tips

**Use `tag` for readable logs.** Adding `tag "${sample_id}"` to a process makes it easy to identify which genome a job belongs to in Nextflow's output and in `work/` directory names.

**Resuming during development.** Add `-resume` when re-running after a code change. Nextflow will re-run only processes whose code or inputs changed, reusing cached results for everything else. This makes iterating on a new step fast.

**Inspecting a failed job.** When a process fails, Nextflow prints the work directory path. That directory contains `.command.sh` (the exact script that ran), `.command.log`, `.command.err`, and all staged input files. You can reproduce the failure with `bash .command.run` from that directory.

**Container note.** If your new process uses a tool not already in an existing container, find the appropriate image on [quay.io/biocontainers](https://quay.io/organization/biocontainers) or [Docker Hub](https://hub.docker.com). Verify the image tag exists and that the tool runs before committing (see the testing approach in the commit history).
