# Running the workflow

What to give the pipeline, what comes back, and how to resume a run that stopped partway.
The steps themselves and the science behind them are in [what the pipeline
does](overview.md), getting it running on your machine is [configuration](configuration.md), and
a scheduler with a lot of genomes is [running at scale](hpc.md).

Every tool runs in a container, so nothing but Nextflow needs installing. Images are pulled on
first use and cached.

---

## Before your first real run

```bash
nextflow run main.nf -profile test,docker
assets/test/assert.sh results_test
```

About 30 seconds on bundled data. It answers whether Nextflow can reach a container engine, and
whether BLIMMP scores anything once it does, which between them is most of what goes wrong with
an installation. Swap `docker` for `singularity` on a cluster. It says nothing about your own
genomes. See [the test dataset](../assets/test/README.md).

---

## Parameters

Start with exactly one input parameter and one HMM parameter:

| | |
|---|---|
| `--input FILE` | a [samplesheet](samplesheets.md), one row per genome, each carrying its own completeness and taxonomy |
| `--genomes PATH` | a directory of nucleotide FASTA files, or a glob. Genes are called with pyrodigal |
| `--proteins PATH` | a directory of protein FASTA files, or a glob. Skips pyrodigal |

and

| | |
|---|---|
| `--hmms PATH` | one concatenated `.hmm` file such as KOfam `profiles.hmm`, or a directory of individual `.hmm` files |

Everything else has a working default. The full list, with types and defaults, is in
[parameters](parameters.md), which is generated from `nextflow_schema.json` so it cannot drift
from what the pipeline accepts. The flags that need more than a line get one below, or in [what the pipeline
does](overview.md) for completeness and priors.

---

## Discovering parameters

```bash
nextflow run main.nf --help
nextflow run main.nf --help_all
```

Both read `nextflow_schema.json`, which describes every parameter the pipeline takes, grouped
the same way [parameters](parameters.md) is. `--help_all` adds the longer explanation that would
otherwise send you here.

A parameter this pipeline does not recognise stops the run, with the nearest real one suggested,
so `--outdr` never writes your results to `results/` while you wait for them elsewhere.
[One dash or two](parameters.md#one-dash-or-two) covers that and the `-resume` versus `--resume`
trap behind it.

---

## Basic usage examples

```bash
# One genome per subdirectory (data/genomes/MED4/*.fna): IDs are the strain names
nextflow run main.nf \
    --genomes 'data/genomes/*/*.fna' \
    --hmms data/kofam/profiles.hmm

# All genomes in one flat directory: IDs come from the file names
nextflow run main.nf \
    --genomes data/genomes/ \
    --hmms data/kofam/profiles.hmm \
    --id_from_filename

# One row per genome, each with its own completeness and priors
nextflow run main.nf \
    --input samples.csv \
    --hmms data/kofam/profiles.hmm

# Already have protein FASTAs, so skip pyrodigal
nextflow run main.nf \
    --proteins data/proteins/ \
    --hmms data/kofam/profiles.hmm \
    --id_from_filename

# Custom output directory
nextflow run main.nf \
    --genomes 'data/genomes/*/*.fna' \
    --hmms data/kofam/profiles.hmm \
    --outdir my_results/
```

A bare directory such as `data/genomes/` matches only the FASTA files directly inside it, not files in subdirectories. Use a glob like `'data/genomes/*/*.fna'` for the one-genome-per-subdirectory layout, and quote it so the shell passes it through.

---

## Output

Each input file gets a subdirectory named after the directory holding it, so a layout of `genomes/<strain>/*.fna` comes back mirrored:

```
results/
  <strain>/
    <sample_id>.tblout                  Per-sequence HMMER hit table
    <sample_id>.domtblout               Per-domain HMMER hit table
    busco_<strain>/                     Full BUSCO output (only if --busco_lineage provided)
    <sample_id>.busco_completeness.tsv  BUSCO completeness fraction and the dataset that measured it (only if --busco_lineage provided)
    <sample_id>_BLIMMP_substituted_modules.json     BLIMMP module presence calls (only if --busco_lineage provided)
    <sample_id>_BLIMMP_substituted_module_probabilities.csv   Per-module probabilities
    <sample_id>_BLIMMP_substituted_module_steps.csv           Per-step presence
    <sample_id>_BLIMMP_substituted_dk.csv                     Per-KO probability estimates
    <sample_id>_BLIMMP_version.txt                The BLIMMP version that produced them
    <sample_id>_debug.log                         BLIMMP debug log (unless --no_blimmp_log)
  module_matrix.tsv         Module confidence per genome, modules as rows and genomes as columns
  summary.tsv               Hit counts per genome at two E-value thresholds
  busco_summary.tsv         BUSCO completeness per genome, with the dataset that measured each (only if --busco_lineage provided)
  pipeline_info/
    timeline.html           Wall-time breakdown per process
    report.html             Resource usage report
```

For example, genomes in `genomes/MED4/GCF_000011465.fna` and `genomes/MIT9313/GCF_000011485.fna` would produce:

```
results/
  MED4/
    GCF_000011465.tblout
    GCF_000011465.domtblout
  MIT9313/
    GCF_000011485.tblout
    GCF_000011485.domtblout
  summary.tsv
```

`module_matrix.tsv` is the file to open first when comparing genomes. One row per KEGG module,
one column per genome, holding the module confidence BLIMMP assigned:

```
module   description                                     MED4    MIT9301
M00001   Glycolysis (Embden-Meyerhof pathway)...         0.445   0.445
M00002   Glycolysis, core module involving three-...     0.61    0.61
```

`NA` means the module is missing from that genome's file rather than scored at zero.

### Sample IDs

The sample ID names the output subdirectory and prefixes every file inside it, so IDs must be unique. By default an ID is the name of the directory holding the input file, which gives `MED4` for `genomes/MED4/GCF_000011465.fna`.

If your genomes sit in one flat directory, every input has the same parent and would get the same ID. The pipeline stops with an error rather than merging them:

```
Duplicate sample IDs derived from parent directory names:
    'genomes' from 9 files
```

Two ways to resolve it. Either give each genome its own directory:

```
genomes/MED4/GCF_000011465.fna
genomes/MIT9313/GCF_000011485.fna
```

Or keep the flat layout and take IDs from file names with `--id_from_filename`:

```
genomes/GCF_000011465.fna    →  results/GCF_000011465/
genomes/GCF_000011485.fna    →  results/GCF_000011485/
```

Use whichever gives the more meaningful names. Directory names carry strain labels; file names usually carry accessions.

---

## Samplesheets

A samplesheet gives each genome its own row, so completeness and taxonomy can differ from one
to the next instead of applying to the whole run:

```csv
name,genome_files,sigma,taxonomy
MED4,genomes/MED4/MED4.fna,0.91,cyanobacteriota
SS120,genomes/SS120.fna,0.87,cyanobacteriota
bin_014,bins/bin_014/*.fa,0.62,
```

```bash
nextflow run main.nf --input samples.csv --hmms /path/to/kofam/
```

Blank `sigma` and `taxonomy` cells fall back to `--sigma` and `--taxonomy`, so anything shared
across rows can still be said once on the command line. A row naming several files is one
sample assembled from pieces, which is the metagenome bin case.

Switching an existing run onto a sheet does not re-run your searches.

Ten worked examples, the full column reference and the error messages are on their own page:
[samplesheets](samplesheets.md).

---

## Incremental runs

Add `-resume` to re-use work already done. Nextflow gives every task a fingerprint taken
from its input files, the command it ran, and the container it ran in. On a re-run it
compares fingerprints: anything that would run identically is copied from `work/` instead of
being computed again, and only the tasks whose fingerprint changed are run.

```bash
nextflow run main.nf \
    --genomes /path/to/genomes/ \
    --hmms /path/to/kofam/ \
    --busco_lineage bacteria_odb12 \
    --outdir results/ \
    -profile uahpc,singularity \
    -resume
```

This covers the common cases without any flags of its own:

Adding genomes runs the new ones and reads the rest from cache.

If you ran without `--busco_lineage` and want it after all, BUSCO and BLIMMP are computed while
the HMMER results are re-used.

Changing a parameter, a script, or a tool version re-runs whatever depended on it. Upgrading
the BLIMMP container re-runs BLIMMP and nothing else.

Masking is the expensive exception. Adding `--no_masking` rewrites the protein sequences, and the
search reads those, so every `HMMSEARCH` task is recomputed. Masking is on by default, so what
costs you here is turning it off partway through a large run.

An interrupted run picks up at the task it was interrupted on, not at the start of the
stage. A job that dies partway through the HMMER search keeps every chunk that had already
finished and repeats only the rest, so hitting a walltime limit costs you the tasks that
were in flight rather than the whole search.

To re-run everything from scratch, omit `-resume`, or delete `work/` to also reclaim the
disk.

### Which run does -resume continue?

`-resume` on its own continues whichever run you started most recently from the directory
you are standing in. Nextflow keeps that list per directory, not per pipeline, so a dry run
or a quick test you launched from the same place is what `-resume` will pick up, even though
you meant the real run.

Nothing warns you when this happens. The run you resumed by accident stored its files
somewhere else, so nothing matches, everything is computed again, and the run succeeds. What
you see is a pipeline that started over for no visible reason, which looks the same as a
cache that had nothing to offer.

Name the run when the directory holds more than one:

```bash
nextflow log                      # run names, session IDs, and the command each used
nextflow run main.nf ... -resume tender_magritte
```

Either a run name or a session ID works. A session ID is the safer choice in a script, since
it stays valid across later resumes of the same chain.

> **Note:** `-resume` reads the cache in `work/`, so keep that directory if you intend to
> re-run. Deleting `work/` and keeping `--outdir` leaves you with results you cannot extend
> incrementally.

If you used an older version of this pipeline, `--force` and the file-existence check it
overrode are both gone. `-resume` does that job properly.

---

## Dry run

`-stub-run` walks the whole pipeline without computing anything. Every process writes empty
placeholder outputs instead of real ones, and `publishDir` still fires, so every output
subdirectory is created. That is how you check the layout before committing to a full run.

```bash
nextflow run main.nf \
    --genomes 'genomes/*/*.fna' \
    --hmms /path/to/kofam/ \
    --outdir results/ \
    -stub-run \
    -profile standard
```

One dash, because this is a Nextflow option rather than one of ours.

Use `-profile standard` for stub runs. It runs locally with no containers, and there is no
reason to pull Singularity images to create directories.

After the run, inspect the output:

```bash
find results/ -type d
```

A stub run splits the HMM database into one dummy chunk rather than the hundreds a real run
uses, so you will see one `HMMSEARCH` task per genome instead of one per genome and chunk.
Trust the directory structure it shows you; ignore the job count.

> **Warning:** Do not point `-stub-run` at the same `--outdir` you intend to use for a real
> run. The placeholder files are published as if they were results and will overwrite real ones.

---

## Running on a cluster

Chunk size is chosen for you and printed at the top of the run, so the common case needs no flag
at all. Submitting the head job, overriding that choice, capping how many jobs reach the queue at
once, and
cleaning up afterwards are all in [running at scale](hpc.md).
