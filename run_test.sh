#!/bin/bash
#
# Template: a small run of your own data, on your own cluster.
#
# Different from `nextflow run main.nf -profile test,docker`: that one checks
# the installation against bundled data in 30 seconds and needs nothing from
# you, while this checks that your cluster, your account, your storage and your
# copy of KOfam work together on a handful of your own genomes.
#
# Replace every YOUR_ and /path/to/ below, then: sbatch run_test.sh
#
#SBATCH --job-name=nf-test
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --output=nextflow-%j.out
#SBATCH --error=nextflow-%j.err

# This job only submits and monitors, so one core is enough. The searches run
# in the jobs Nextflow submits from here.

# However your cluster provides Nextflow. `module avail nextflow` will say, and
# docs/configuration.md covers installing it yourself if it does not.
module load nextflow

nextflow run main.nf \
    --genomes /path/to/genomes/ \
    --hmms /path/to/kofam/ \
    --subsample_genomes 100 \
    --subsample_hmms 100 \
    --max_concurrent 20 \
    -profile slurm,singularity \
    -c conf/mysite.config
