# Test dataset

What `-profile test` runs. Small enough to finish in a couple of minutes, real enough that
finishing means something.

```bash
nextflow run main.nf -profile test,docker
```

## What is here

| | |
|---|---|
| `genomes/MED4/MED4.fna` | 150 kb of *Prochlorococcus marinus* MED4, GCF_000011465.1, from position 141,959 |
| `genomes/MIT9301/MIT9301.fna` | 150 kb of *Prochlorococcus marinus* MIT9301, GCF_000015965.1, from position 152,723 |
| `kofam/*.hmm` | 20 KOfam profiles out of the full 27,262 |

About 2.4 MB in total, most of it the HMM profiles.

## Why these regions and these profiles

Picked so the test produces real results. An empty output file exits 0 and looks like a working
run, and that has happened here before: BLIMMP 0.1.2 shipped without its KEGG data and wrote output
that was structurally perfect and entirely empty for four months.

The 20 KOs were chosen from BLIMMP's own confident calls on the full MED4 genome, filtered to
`KO_conf_final >= 0.9`, then narrowed to a 150 kb window holding as many of them as possible,
then to the 20 with the smallest profiles so the bundle stays under a few megabytes.

Both regions cover all 20. They are conserved core genes in a syntenic stretch, so the same
window in either strain carries the same set, which is why two genomes could be included
without the second one coming back empty and failing `CHECK_BLIMMP`.

## Checking the result

```bash
nextflow run main.nf -profile test,docker
assets/test/assert.sh results_test
```

`assert.sh` asserts numbers. Per genome it requires that at
least half the bundled KOfam profiles hit at E<1e-6, that all four BLIMMP output files exist and
are non-empty, and that the per-KO table has rows. It also checks that `summary.tsv` names every
sample.

Both test genomes hit 20 of 20, so the half threshold has room in it: one profile drifting below
threshold reports 19/20 and still passes, while a broken search, a mis-split HMM database or the
wrong genome slice fails.

It was checked against deliberately broken output as well as good: an emptied `.tblout` fails,
and a `dk.csv` cut down to its header alone fails with "BLIMMP scored nothing". That second case
is the BLIMMP 0.1.2 failure reproduced on purpose, because a check that cannot fail is worth
nothing.

## What it does and does not prove

Proves your setup works: containers pull and run, pyrodigal calls genes, the HMM database splits
into chunks, hmmsearch runs across the genome-by-chunk grid, results merge, BLIMMP scores them,
and the check confirms it scored something.

Does not prove anything biological. Twenty profiles out of 27,262 cannot fill a KEGG module,
and `--sigma 0.95` is a stand-in: a 150 kb slice is not 95% of a genome. Module confidences
from this dataset mean nothing.

## Rebuilding it

```
K00330 K00339 K00340 K00762 K00860 K00950 K01588 K01658 K01693 K02707 K02708 K03465 K05574 K05576 K05578 K05580 K05581 K05582 K05584 K12073
```

The selection came from `results_prochlorococcus/*/` in the development notebook, using each
strain's `_BLIMMP_dk.csv` for confident KO calls and the pyrodigal `.faa` headers for ORF
coordinates. Regenerating needs a full run against complete KOfam first.
