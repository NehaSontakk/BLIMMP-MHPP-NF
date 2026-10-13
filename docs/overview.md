# What the pipeline does

For each genome, protein-coding genes are predicted with
[pyrodigal](https://github.com/althonos/pyrodigal), then searched against the full
[KOfam](https://www.genome.jp/tools/kofamkoala/) database with HMMER. Optionally, genome
completeness is estimated with BUSCO and passed to
[BLIMMP](https://github.com/TravisWheelerLab/BLIMMP) as a correction factor. BLIMMP then infers
the probability that each KEGG module is present, accounting for what the assembly is missing.

The work parallelises on two axes at once: every genome against every chunk of the HMM database.
Nine genomes and 270 chunks is 2,430 independent searches, which is why this is a Nextflow
pipeline rather than a shell script.

---

## The steps

Two lines of work start from the same genomes and run at the same time. One annotates the genes,
the other measures how complete each genome is. They meet only at the end, where BLIMMP needs a
result from both.

**Annotation**

1. `CONCAT_INPUTS`: only when a [samplesheet](samplesheets.md) row names several files; joins
   them into one sample
2. `PREDICT_PROTEINS`: pyrodigal turns nucleotide FASTA into protein FASTA. Skipped if you
   supply `--proteins`
3. `MASK_REPEATS`: tantan replaces low-complexity residues with `X`, which HMMER scores as
   unknown. On by default; `--no_masking` skips it
4. `GROUP_HMMS` / `SPLIT_HMMS`: divide the HMM database into chunks, as many as it takes to keep
   your machine or your queue busy. `SPLIT_HMMS` takes one concatenated `.hmm` file, `GROUP_HMMS`
   takes a directory of them. The run prints the size it chose and what to pass to change it; see
   [parallelism and chunk size](hpc.md#parallelism-and-chunk-size)
5. `HMMSEARCH`: every genome against every chunk, one job each
6. `MERGE_RESULTS`: concatenate the per-chunk results into one `.tblout` and one `.domtblout`
   per genome
7. `POSTPROCESS`: hit counts per genome, into `summary.tsv`

**Completeness**, only with `--busco_lineage`

8. `RUN_BUSCO`: assess each genome
9. `BUSCO_SUMMARY`: collect the values into `busco_summary.tsv`

**Module detection**

10. `RUN_BLIMMP`: per genome, from the merged `.domtblout` and a completeness value. Also writes
    a debug log and the BLIMMP version that produced the numbers
11. `CHECK_BLIMMP`: confirm BLIMMP scored something, and stop the run if it did not
12. `MODULE_MATRIX`: pivot the per-genome confidences into `module_matrix.tsv`, modules as rows
    and genomes as columns

---

## Why the proteins are masked

Low-complexity stretches, runs of one or two residues repeated, match many HMM profiles for
reasons that have nothing to do with homology. tantan replaces them with `X`, which HMMER scores
as an unknown residue, so those regions stop contributing score in either direction.

This is on by default. `--no_masking` turns it off, and changing it either way re-runs every
search, so it is a decision to make before a large run rather than partway through one.

The default was set by running nine *Prochlorococcus* strains both ways against the full KOfam
database. Masking removed 209 hits, and not one of them was above its KOfam threshold. It removed
hits scoring 7 against a threshold of 377 that were nonetheless being reported as present.

It also recovers hits, which is the less obvious half. KOfam scores about a third of its KOs on
the best single domain rather than on the full sequence. A low-complexity stretch sitting inside a
protein can break the alignment into fragments that each score below the threshold, while the
full-sequence score sits far above it. Masking the stretch lets HMMER align one long domain
instead. Threonine dehydratase went from a best domain of about 35 to about 485, against a
threshold of 261.5, in seven of the nine strains; aspartate kinase did the same in four.

The cost is that masking changes results. Thirteen module calls crossed 0.5 across those nine
strains, so results produced with masking are not directly comparable to results produced without
it. If you are comparing against numbers from an earlier run, check which setting produced them.

---

## How complete are your genomes?

BLIMMP has to know this, because a KEGG module missing a step means one thing in a finished
genome and something else in a draft that is missing a tenth of its genes. That fraction is
sigma, and you supply it one of two ways.

**If you are reasonably certain your genomes are complete, use `--sigma 1.0`.** Closed
assemblies, finished RefSeq genomes, anything you would describe as a genome rather than a
draft.

`--sigma 0.85` states one value for every input. Use it when you know the fraction, and when you
are starting from proteins, where BUSCO has nothing to measure and `--sigma` is the only route
to BLIMMP at all.

`--busco_lineage bacteria_odb12` measures sigma per genome with BUSCO. Reach for it when the
missing fraction is real and unknown: drafts, MAGs, anything assembled from short reads and left
in contigs.

`--busco_lineage auto_prok` measures it without your naming a clade. BUSCO scores each genome
against the domain root sets, places the markers it finds on a reference tree, and rescores
against the deepest dataset the placement supports. A directory of bins from a metagenome, each
from a different clade, is what this is for: each bin gets its own marker set, and
`busco_summary.tsv` records which one in its third column. A bin too fragmentary to place keeps the
domain root set, and the third column then says `bacteria_odb12.2` rather than a clade. On a 150
kb slice BUSCO placed 3 of its 47 markers and kept the root; on the whole genome it placed 35 and
chose the genus. BUSCO runs at least twice per genome
this way, once against each root set and once against the chosen dataset, with the placement in
between, so expect a few times the run time of a named dataset; `auto` searches both domains and
costs more than `auto_prok` or `auto_euk`. BUSCO may then measure two rows of that table against
different marker sets, which is fine when each sigma only has to be right for its own bin.

On a genome that is already complete, BUSCO does not return 1.0. It returns how much of one
marker set it found, and a marker set is drawn from a clade rather than from your organism. Nine
closed Prochlorococcus genomes score anywhere from 0.86 to 0.99 depending only on which of four
lineage ranks you name. Handing those numbers to BLIMMP asserts that some of each genome is
missing when none of it is, which raises every module presence probability, and it does so
unevenly: the genomes that diverge most from the marker set are penalised hardest, which is
usually the opposite of what you want.

Supply one or the other. With neither, the pipeline annotates against KOfam, stops there, and
warns you at the start that BLIMMP is not going to run. That is a perfectly good way to use it
if annotation is all you wanted.

---

## Taxonomic priors

Before it sees your data, BLIMMP starts from how often each KO turns up in genomes like yours. Those background frequencies are the priors, and BLIMMP ships several sets of them: one built from all bacteria, four at the rank of kingdom, and fifteen at the rank of phylum.

Without `--taxonomy` you get the all-bacteria set. Naming a group closer to your organisms gives BLIMMP a better starting expectation:

```bash
nextflow run main.nf ... --taxonomy cyanobacteriota
```

Short forms work when they are unambiguous, so `cyanobacteriota` finds the set whose full name is `cyanobacteriota_melainabacteria_group`. A name BLIMMP does not recognize stops the run and prints the names it does, grouped by rank. That check lives in BLIMMP rather than here, so the pipeline never carries a second copy of the list to fall out of date, and the failure arrives when the first BLIMMP job runs rather than at startup.

The choice is not free: on nine *Prochlorococcus* strains, moving from the all-bacteria priors to the phylum priors moved 31 module probabilities by more than 0.10, and around 10 per strain across the 0.5 mark. A narrower group is not automatically the better answer, since those sets are built from fewer genomes. The debug log records which set was used for each genome.

One setting applies to every input, so a run mixing organisms from different phyla should either use the bacteria-wide priors or be split into separate runs.

---

## One value per run, or one per genome

`--sigma` and `--taxonomy` each set one value for the whole run. When your genomes differ from
each other, describe them in a table and pass it with `--input`:

```csv
name,genome_files,sigma,taxonomy
MED4,genomes/MED4/MED4.fna,0.91,cyanobacteriota
SS120,genomes/SS120.fna,0.87,cyanobacteriota
bin_014,bins/bin_014/*.fa,0.62,
```

```bash
nextflow run main.nf --input samples.csv --hmms /path/to/kofam/
```

Each row names one sample and the file or files it comes from. `sigma` and `taxonomy` are
optional, and a blank cell falls back to whatever the command line gave. The third row is a
metagenome bin spread across several files, joined and gene-called once.

Moving an existing run onto a sheet does not re-run your searches. Ten worked examples are in
the [samplesheet guide](samplesheets.md).

---

## Re-scoring finished searches

The HMMER search takes hours; BLIMMP takes minutes. When you want to score the same genomes again under different taxonomic priors, or when you already have `.domtblout` files from somewhere else, there is a second way in that skips straight to BLIMMP:

```bash
nextflow run main.nf -entry blimmp_only \
    --domtblout 'results/*/*.domtblout' \
    --completeness results/busco_summary.tsv \
    --taxonomy cyanobacteriota \
    --outdir results_cyano/
```

`--domtblout` takes a directory of files or a glob, and sample IDs are derived the same way as everywhere else: from the parent directory, or from file names with `--id_from_filename`.

Sigma comes from one of two places. `--completeness` reads the `busco_summary.tsv` written by an earlier run, so every genome keeps the value BUSCO measured for it, and the result is identical to what the full pipeline would produce. `--sigma` sets one value for all of them instead. Give one or the other.

Flags belonging to the search stop the run rather than being ignored, since nothing here would use them: `--genomes`, `--proteins`, `--hmms`, `--busco_lineage`.

---

## Debug logs

BLIMMP keeps a running account of what it loaded and how it scored each KO, and the pipeline
saves that file next to the results for the genome it came from:

```
results/<name>/<name>_debug.log
```

It is the only file that records which completeness value and which taxonomic priors that
genome was scored with, and it is what tells you when BLIMMP ran but had nothing to work
with. An early container shipped without its KEGG data files and wrote results that looked
perfectly normal and were empty; the log said so on its first line and nowhere else did.

The cost is about 3.4 MB per genome. `--no_blimmp_log` turns it off — and, since BLIMMP_MHPP ties
the log to `-v`/`--verbose`, also skips writing the pass-1 comparison files alongside the
substituted ones the pipeline reads.

If BLIMMP fails partway through, the pipeline copies out whatever the log holds before the
run stops, so you keep the account of what happened up to the failure. That copy needs the
output directory to be reachable from inside the container, which is not true on every
system. When it isn't, the pipeline prints the working directory the log is sitting in so
you can go get it.

---
