#!/bin/bash
#
# Template: one genome, one HMM chunk, verbose logging.
#
# For when a run fails on your cluster and you need to see what Nextflow is
# handing the scheduler. NXF_DEBUG=2 makes it print the command line
# for every task it submits, which is where a wrong account, a missing module or
# a memory flag your cluster rejects becomes visible.
#
# The other half of that evidence is in the failing task's work directory:
# .command.sh, .command.err and .command.run. docs/reporting-problems.md says
# which file holds what, and what to send us if it still makes no sense.
#
# Replace every YOUR_ and /path/to/ below, then: sbatch run_slurm_debug.sh
#
#SBATCH --job-name=nf-debug
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=nf-debug-%j.out
#SBATCH --error=nf-debug-%j.err

module load nextflow

NXF_DEBUG=2 nextflow run main.nf \
    --genomes /path/to/genomes/ \
    --hmms /path/to/kofam/ \
    --subsample_genomes 1 \
    --subsample_hmms 1 \
    -profile slurm,singularity \
    -c conf/mysite.config
