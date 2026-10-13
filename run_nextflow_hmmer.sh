#!/bin/bash
#
# Template: a full production run.
#
# Nextflow runs here as a head job. It submits one job per genome and per HMM
# chunk, waits for them, and submits the next, so this allocation needs one core
# and very little memory however large the run is. --time is the number that
# matters: if this job hits the limit the pipeline stops, even with hundreds of
# child jobs still running.
#
# Replace every YOUR_ and /path/to/ below, then: sbatch run_nextflow_hmmer.sh
#
#SBATCH --job-name=nf-hmmsearch
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=nextflow-%j.out
#SBATCH --error=nextflow-%j.err

module load nextflow

# Pull each container image once, into somewhere shared and persistent, rather
# than into $HOME on every run.
export NXF_SINGULARITY_CACHEDIR=/path/to/shared/singularity_cache

# --max_concurrent caps how many jobs sit in the queue at once. Your cluster's
# per-user job limit is a reasonable starting value.
# -resume re-uses everything already in work/, so a run that stopped partway
# continues rather than starting over. It is harmless on a first run.
nextflow run main.nf \
    --genomes /path/to/genomes/ \
    --hmms /path/to/kofam/ \
    --busco_lineage bacteria_odb12 \
    --outdir results/ \
    --max_concurrent 200 \
    -profile slurm,singularity \
    -c conf/mysite.config \
    -resume

# If you would rather claim one large node and run everything inside it than
# dispatch to the scheduler, drop the #SBATCH allocation above to a whole node,
# and swap the two profile lines for `-profile server,singularity`. That keeps
# the work off the queue entirely, at the cost of the parallelism a cluster
# gives you.
