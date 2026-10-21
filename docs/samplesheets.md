# Samplesheets

`--genomes` and `--proteins` describe your inputs by where they sit on disk, and every other
setting applies to the whole run. That stops working as soon as your genomes differ from each
other. A draft assembly and a closed genome do not share a completeness, and a mixed-taxon set
does not share a set of priors.

`--input samples.csv` describes them one row at a time.

- [The columns](#the-columns)
- [Worked examples](#worked-examples)
  - [1. The smallest useful sheet](#1-the-smallest-useful-sheet)
  - [2. Each genome its own completeness](#2-each-genome-its-own-completeness)
  - [3. Several organisms, several sets of priors](#3-several-organisms-several-sets-of-priors)
  - [4. Metagenome bins spread across files](#4-metagenome-bins-spread-across-files)
  - [5. Starting from proteins](#5-starting-from-proteins)
  - [6. Genomes and proteins in one run](#6-genomes-and-proteins-in-one-run)
  - [7. Measuring with BUSCO and overriding two rows](#7-measuring-with-busco-and-overriding-two-rows)
  - [8. Annotating some rows without scoring them](#8-annotating-some-rows-without-scoring-them)
  - [9. Tab-separated instead of comma-separated](#9-tab-separated-instead-of-comma-separated)
  - [10. Moving an existing run onto a sheet](#10-moving-an-existing-run-onto-a-sheet)
- [Rules worth knowing](#rules-worth-knowing)
- [When something is wrong with the sheet](#when-something-is-wrong-with-the-sheet)

---

## The columns

| column | required | holds |
|---|---|---|
| `name` | yes | the sample ID: names the output subdirectory and prefixes every file in it |
| `genome_files` | one of the two | nucleotide FASTA, gene-called with pyrodigal |
| `protein_files` | one of the two | protein FASTA, skipping pyrodigal |
| `sigma` | no | this genome's completeness fraction |
| `taxonomy` | no | the priors to score this genome against |

Column order does not matter, and you can leave out any column you have nothing to put in.

---

## Worked examples

### 1. The smallest useful sheet

Two columns. Everything the genomes have in common stays on the command line.

```csv
name,genome_files
MED4,genomes/MED4.fna
MIT9313,genomes/MIT9313.fna
SS120,genomes/SS120.fna
```

```bash
nextflow run main.nf \
    --input samples.csv \
    --hmms /path/to/kofam/ \
    --sigma 0.95 \
    --taxonomy cyanobacteriota
```

All three genomes are scored at 95% complete against cyanobacterial priors. This is worth
doing even when nothing varies, because the sheet states which files went in and what they
were called, where a glob only records where you pointed on the day.

### 2. Each genome its own completeness

The usual reason to reach for a sheet. These are assemblies of different quality, and scoring
a 62% complete draft as though it were finished would report modules missing that are probably
just unsequenced.

```csv
name,genome_files,sigma
closed_genome,assemblies/closed.fna,1.0
good_draft,assemblies/draft_a.fna,0.94
rough_draft,assemblies/draft_b.fna,0.71
very_rough,assemblies/draft_c.fna,0.62
```

```bash
nextflow run main.nf --input samples.csv --hmms /path/to/kofam/
```

No `--sigma` on the command line, because every row carries its own.

### 3. Several organisms, several sets of priors

BLIMMP scores against background expectations of how often each KO turns up. One run can now
use a different background per genome.

```csv
name,genome_files,sigma,taxonomy
MED4,cyano/MED4.fna,0.91,cyanobacteriota
MIT9313,cyano/MIT9313.fna,0.95,cyanobacteriota
E_coli_K12,entero/K12.fna,1.0,pseudomonadota
B_subtilis,firm/subtilis.fna,1.0,bacillota
unknown_isolate,misc/isolate_7.fna,0.88,
```

```bash
nextflow run main.nf --input samples.csv --hmms /path/to/kofam/
```

`unknown_isolate` has an empty taxonomy, so it falls back to the all-bacteria priors. That is
the right choice when you do not know: a narrower prior is not automatically a better one, and
guessing wrong shifts your results. See [Taxonomic priors](overview.md#taxonomic-priors).

### 4. Metagenome bins spread across files

A bin usually arrives as a pile of contig files rather than one assembly. Point one row at all
of them with a glob and they are treated as one organism.

```csv
name,genome_files,sigma,taxonomy
bin_003,bins/bin_003/*.fa,0.78,
bin_014,bins/bin_014/*.fa,0.62,
bin_022,bins/bin_022/*.fa,0.91,cyanobacteriota
isolate_ref,refs/reference.fna,1.0,cyanobacteriota
```

```bash
nextflow run main.nf --input samples.csv --hmms /path/to/kofam/
```

The pieces of each bin are joined into one sequence and gene-called once. That is also the
better answer scientifically: pyrodigal trains its gene model on whatever you hand it, so
training once over a whole bin beats training separately on each contig file, and the gap
widens as the contigs get shorter.

Bins and ordinary assemblies can sit in the same sheet, as `isolate_ref` does here. A row with
one file is just a row whose glob matched once.

### 5. Starting from proteins

Use `protein_files` when you already have predicted proteins and want to skip gene calling.

```csv
name,protein_files,sigma,taxonomy
proteome_A,proteomes/A.faa,1.0,pseudomonadota
proteome_B,proteomes/B.faa,1.0,pseudomonadota
mag_proteins,proteomes/mag_*.faa,0.83,
```

```bash
nextflow run main.nf --input samples.csv --hmms /path/to/kofam/
```

Without `--busco_lineage`, a `sigma` value is the only way a protein row reaches BLIMMP — add
`--busco_lineage` and BUSCO measures these too, in `proteins` mode instead of `genome` mode (see
example 7). `mag_proteins` shows that a glob works here as well: several protein files get
joined into one proteome, which needs none of the gene-model care that joining nucleotides does.

### 6. Genomes and proteins in one run

One sheet can hold both. Fill in `genome_files` or `protein_files` per row, never both.

```csv
name,genome_files,protein_files,sigma,taxonomy
MED4,genomes/MED4.fna,,0.91,cyanobacteriota
MIT9313,genomes/MIT9313.fna,,0.95,cyanobacteriota
external_A,,proteomes/A.faa,1.0,cyanobacteriota
external_B,,proteomes/B.faa,1.0,cyanobacteriota
```

```bash
nextflow run main.nf --input samples.csv --hmms /path/to/kofam/
```

The genome rows go through pyrodigal; the protein rows skip it. Everything downstream is the
same for both.

### 7. Measuring with BUSCO and overriding two rows

`--busco_lineage` works alongside a sheet, on both `genome_files` and `protein_files` rows alike
— BUSCO measures a genome row in `genome` mode and a protein row in `proteins` mode, and either
way fills in the blank `sigma` cells. A cell you filled in beats the measurement. With bins from
several clades, name `auto_prok` instead of a dataset and BUSCO chooses one per row.

```csv
name,genome_files,protein_files,sigma,taxonomy
MED4,genomes/MED4.fna,,,cyanobacteriota
MIT9313,genomes/MIT9313.fna,,,cyanobacteriota
SS120,genomes/SS120.fna,,,cyanobacteriota
odd_assembly,genomes/odd.fna,,0.55,cyanobacteriota
known_complete,genomes/finished.fna,,1.0,cyanobacteriota
external_proteome,,proteomes/A.faa,,cyanobacteriota
```

```bash
nextflow run main.nf \
    --input samples.csv \
    --hmms /path/to/kofam/ \
    --busco_lineage bacteria_odb12
```

The first three take whatever BUSCO measures. `odd_assembly` and `known_complete` use the
numbers you typed, on the reasoning that you had a reason for typing them. `external_proteome`
has no `genome_files` to call genes from, so BUSCO measures its proteins directly instead — same
column, same fallback order, just a different BUSCO mode under the hood. BUSCO still runs on all
six rows, so its summary table stays complete.

The precedence is the cell, then BUSCO's measurement, then `--sigma`.

### 8. Annotating some rows without scoring them

A row with no completeness from any source is searched against KOfam and stops there.

```csv
name,genome_files,sigma
scored_1,genomes/a.fna,0.94
scored_2,genomes/b.fna,0.88
annotate_only,genomes/c.fna,
```

```bash
nextflow run main.nf --input samples.csv --hmms /path/to/kofam/
```

`annotate_only` gets its `.tblout` and `.domtblout` and no BLIMMP output. The pipeline says so
before it starts:

```
WARN: BLIMMP will not run for annotate_only: no completeness estimate. Fill in the
sigma column on those rows, or pass --sigma. They are annotated against KOfam and
stop there.
```

Add `--sigma 0.9` to that command and `annotate_only` would pick it up, since blank cells fall
back to the flag.

### 9. Tab-separated instead of comma-separated

Name the file `.tsv` and it is read as tab-separated. Useful when your paths or names contain
commas.

```tsv
name	genome_files	sigma	taxonomy
MED4	genomes/MED4.fna	0.91	cyanobacteriota
MIT9313	genomes/MIT9313.fna	0.95	cyanobacteriota
```

```bash
nextflow run main.nf --input samples.tsv --hmms /path/to/kofam/
```

In a `.csv`, a field containing a comma has to be quoted, which includes brace globs:

```csv
name,genome_files,sigma
MED4,"genomes/MED4/*.{fa,fna}",0.91
```

### 10. Moving an existing run onto a sheet

Switching does not re-run your searches. A sheet naming the files you already searched
produces the same task hashes they had under `--genomes`, so `-resume` reuses all of it.

If an earlier run used `--busco_lineage`, its `busco_summary.tsv` already holds the per-genome
completeness, and you can build the sheet from it:

```bash
{
  echo "name,genome_files,sigma,taxonomy"
  tail -n +2 results/busco_summary.tsv | while read -r id sigma; do
      echo "${id},/data/genomes/${id}/*.fna,${sigma},cyanobacteriota"
  done
} > samples.csv
```

```bash
nextflow run main.nf \
    --input samples.csv \
    --hmms /path/to/kofam/ \
    --outdir results_cyano/ \
    -resume
```

Measured on nine genomes: `PREDICT_PROTEINS` 9 of 9 cached, `HMMSEARCH` 2,457 of 2,457 cached,
`MERGE_RESULTS` 9 of 9 cached, with only BLIMMP re-running. With the same sigma values the
BLIMMP output was identical to the `--busco_lineage` run across all 36 files.

If you only want to re-score and have the `.domtblout` files already, `-entry blimmp_only` is
the shorter road. See [Re-scoring finished searches](overview.md#re-scoring-finished-searches).

---

## Rules worth knowing

**One row is one sample.** The `name` gives the output subdirectory and prefixes every file in
it, so two rows cannot share a name.

**A cell holds one path or a glob, never a bare directory.** A glob says which files it means,
where a directory leaves it open whether the other things in there were meant to come along.

**Files in a multi-file row are joined in sorted order.** Glob expansion order is not
guaranteed across filesystems, and the order reaches the task hash, so sorting keeps the same
sheet giving the same answer everywhere.

**Relative paths are read from the sheet's own directory**, not from wherever you launched
Nextflow, so a sheet and its data can be moved or copied together. This is the one place the
sheet deliberately differs from the command line, where a relative `--genomes` path is read
from the launch directory as usual.

**Blank cells fall back to the flags.** Anything shared across every row can still be said once
on the command line, and the sheet only has to carry what differs.

**`--input` cannot be combined with `--genomes` or `--proteins`.** The sheet already names every
input file.

**`-entry blimmp_only` does not take a sheet.** It starts from finished HMMER tables and uses
`--domtblout` and `--completeness`. Passing `--input` there is refused rather than accepted and
ignored.

---

## When something is wrong with the sheet

Every check happens before a single task is submitted, and every message names the row and the
line it is on. A few of the more common ones:

A directory where a glob belongs:

```
--input samples.csv: genome_files for 'MED4' is a directory: genomes/MED4
  Name the files with a glob instead, such as genomes/MED4/*.fna, so the sheet
  says which of the files in there belong to this sample.
```

A misspelled column, which stops the run rather than being ignored, so the typo costs you the
message instead of the column it was meant to fill:

```
--input samples.csv has columns this pipeline does not read: taxonmy
  Known columns are name, genome_files, protein_files, sigma, taxonomy.
```

A glob that matched nothing, which also tells you where relative paths are being read from:

```
--input samples.csv: genome_files for 'X' matched no files: genomes/X/*.fna
  Relative paths are read from the samplesheet's own directory, which is
  /data/project/sheets.
```

Both input columns filled on one row, neither filled, a repeated name, a missing `name` column,
and a `sigma` outside 0 to 1 are each caught the same way.
