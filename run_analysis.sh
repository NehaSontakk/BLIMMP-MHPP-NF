#!/bin/bash
#
# Template: a real analysis, pinned so you can resume it.
#
# The other templates finish in one sitting. This one is for the run you come
# back to: a study's worth of genomes, completeness measured per genome, and a
# taxonomy chosen on purpose. Runs like that get interrupted, and resuming one
# correctly is most of what this file is about.
#
# Replace every YOUR_ and /path/to/ below, then: sbatch run_analysis.sh
#
#SBATCH --job-name=blimmp-analysis
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=nextflow-%j.out
#SBATCH --error=nextflow-%j.err

module load nextflow

export NXF_SINGULARITY_CACHEDIR=/path/to/shared/singularity_cache

# Resume by session id, never bare.
#
# `-resume` on its own continues whichever run is newest in .nextflow/history.
# A directory that has hosted a stub run or a quick test since your last real
# one will offer that instead, and because those used a different work
# directory, the resume matches nothing. Nothing warns you. The run recomputes
# everything and looks like a cache that simply had nothing to offer.
#
# Leave SESSION empty for the first run. Afterwards read the id out of
# `nextflow log` and paste it here, and every later run continues this one.
#
# A session is spent by anything that changes task hashes: a new container, a
# changed script block, a flag that alters a command line. When that happens the
# old id stops matching and the work is recomputed whatever you do, so record
# what spent it rather than resuming into confusion.
SESSION=

# --busco_downloads is a lineage you fetched beforehand, which keeps every BUSCO
# task off the network. Without it each task downloads its own copy, once per
# genome. Populate it once, from anywhere with network access:
#
#     busco --download YOUR_LINEAGE --download_path /path/to/busco_downloads
#
# `busco --list-datasets` prints the names. A narrower lineage measures
# completeness against more markers, which is a different measurement rather
# than a more accurate one, so pick the rank you can defend and say why.
#
# --taxonomy picks which priors BLIMMP scores against. It changes results, and
# a narrower group is not automatically better. Left off, it uses all bacteria.
#
# No --hmm_chunk_size here on purpose: the pipeline sizes chunks from
# --max_concurrent and the number of genomes, and logs what it chose. Set it
# only when you have a reason the log does not already account for.
nextflow run main.nf \
    --genomes '/path/to/genomes/*/*.fna' \
    --hmms /path/to/kofam/ \
    --busco_lineage YOUR_LINEAGE \
    --busco_downloads /path/to/busco_downloads \
    --taxonomy YOUR_TAXONOMY \
    --max_concurrent 50 \
    --outdir results/ \
    -profile slurm,singularity \
    -c conf/mysite.config \
    ${SESSION:+-resume "${SESSION}"}

# The genomes path is quoted because it is a glob the pipeline expands, not the
# shell. Unquoted, bash expands it first and Nextflow receives only the first
# match.
