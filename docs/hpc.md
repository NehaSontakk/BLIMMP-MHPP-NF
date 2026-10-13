# Running at scale

Everything here assumes a scheduler and a lot of genomes. For a first run, or a run on one
machine, start with [running the workflow](running.md).

The pipeline parallelises on two axes at once: every genome against every chunk of the HMM
database. That is what makes it fast, and it is also how it can swamp a shared queue, so chunk
size and the concurrency cap are worth understanding before a large run.

---

## Submitting the head job

The Nextflow process does no analysis. It submits jobs, waits for them, and submits the next
ones, so it wants one core and very little memory. Everything expensive happens in the jobs it
submits. What it does need is to outlive your login session, and clusters differ on where it is
allowed to run.

If your cluster lets you run Nextflow on the login node, launch it there and detach:

```bash
module load nextflow   # if needed

nohup nextflow run main.nf \
    --genomes /path/to/genomes/ \
    --hmms /path/to/kofam/ \
    --max_concurrent 200 \
    -profile slurm,singularity \
    -c conf/mysite.config \
    > nextflow.out 2>&1 &
```

A `screen` or `tmux` session does the same job without `nohup`.

If your cluster forbids that, or only offers Nextflow as a module on compute nodes, wrap the
same command in an sbatch script and submit it as an ordinary job. Four templates ship with the
repository, each with placeholders wherever a value has to be yours:

| Script | What it is |
|---|---|
| `run_test.sh` | A hundred of your genomes against a hundred profiles, to prove the cluster works |
| `run_nextflow_hmmer.sh` | A full production run, with `-resume` and a shared container cache |
| `run_analysis.sh` | A study you will resume: offline BUSCO, a taxonomy, and a pinned session |
| `run_slurm_debug.sh` | One genome, one chunk, and the command line for every job submitted |

Copy the closest one, fill it in, and submit:

```bash
cp run_test.sh run_myproject.sh
# replace every YOUR_ and /path/to/ in it
sbatch run_myproject.sh
```

Then watch it:

```bash
squeue --job <JOBID>
tail -f nextflow-<JOBID>.out
```

Give the head job a `--time` well beyond the expected total. It is idle almost all of that time,
but if it hits the limit the pipeline stops even though its child jobs are still running, and
you are left restarting from the cache rather than finishing.

---

## Three clusters, side by side

A Slurm cluster is not a standard. Two of them will disagree about how memory is asked for, what
the partitions are called, whether a QoS is required, how containers are provided, and what an
account looks like. The templates above carry placeholders in each of those spots because no
value would be right everywhere, and a script that runs here can fail on the next cluster over
for reasons that have nothing to do with this pipeline.

So here is the same run written for three of them. The first is the cluster this pipeline was
developed on and runs on daily. The second, CU Boulder Alpine, was ported by an agent working
through [docs/porting.md](porting.md): the bundled test dataset passes there, and so does the
nine-genome reference set against the full KOfam, in 7m18s over 102 tasks. The first attempt at
the config below, written straight from CU Boulder's docs and never run, was wrong in one place,
and the port turned up two things worth reading before a first run elsewhere. Both are noted
where they happened.

The third, TACC Stampede3, was ported the same way, and it is the one that shows how far apart
two Slurm clusters can be. The first two differ in what to call things. Stampede3 differs in
what a job is: it does not share nodes, so the one-job-per-task arrangement both other
clusters use would bill a whole 48-core node for a task asking for four. It is the first site
here whose config sets `executor = 'local'`, and the first where `--mem` is range-checked and
then discarded.

### University of Arizona HPC, which this was developed on

```bash
#!/bin/bash
#SBATCH --job-name=blimmp
#SBATCH --account=twheeler
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --output=nextflow-%j.out
#SBATCH --error=nextflow-%j.err

source /etc/profile.d/lmod.sh
module load nextflow/22.10.7.5853

export NXF_SINGULARITY_CACHEDIR=/xdisk/$USER/singularity_cache

# Pin the session rather than using bare -resume. -resume on its own continues
# whichever run started most recently from this directory, and a dry run or a
# quick test launched from the same place will be picked up instead of the one
# you meant. Nothing warns you; the run recomputes everything and looks fine.
SESSION=95512764-44bc-44b7-b2b6-d74239ec85bb

nextflow run main.nf \
    --genomes '/path/to/genomes/*/*.fna' \
    --hmms /path/to/kofam/ \
    --busco_lineage bacteria_odb12 \
    --hmm_chunk_size 100 \
    --max_concurrent 50 \
    --outdir results_prochlorococcus \
    -profile uahpc,singularity \
    -resume "${SESSION}"
```

`--mem-per-cpu` is not a style preference here. Slurm on this cluster ignores `#SBATCH --mem`
and silently allocates 1 MB, so asking for memory the ordinary way produces jobs killed for
exceeding a limit you believed you had raised. `-profile uahpc` handles the matching half inside
the pipeline, clearing Nextflow's `memory` directive so it does not emit `--mem` beside the
`--mem-per-cpu`, which sbatch would reject outright.

### CU Boulder Alpine, reached through [RMACC](https://www.colorado.edu/partnerships/rmacc/)

```bash
#!/bin/bash
#SBATCH --job-name=blimmp
#SBATCH --account=ucb###_asc1
#SBATCH --partition=acpu
#SBATCH --qos=cpu-long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=nextflow-%j.out
#SBATCH --error=nextflow-%j.err

module load slurm/alpine

# Apptainer is installed on every Alpine compute node with no module to load,
# and provides the `singularity` command that -profile singularity calls.
export NXF_SINGULARITY_CACHEDIR=/scratch/alpine/$USER/singularity_cache

nextflow run main.nf \
    --genomes '/path/to/genomes/*/*.fna' \
    --hmms /path/to/kofam/ \
    --busco_lineage bacteria_odb12 \
    --max_concurrent 200 \
    -profile slurm,singularity \
    -c conf/alpine.config
```

with `conf/alpine.config`:

```groovy
process {
    queue          = 'acpu'
    clusterOptions = '--account=ucb###_asc1 --qos=cpu-normal'

    withLabel: 'process_low'    { memory = { 4.GB  * task.attempt } }
    withLabel: 'process_medium' { memory = { 16.GB * task.attempt } }
}
```

Alpine honours `--mem`, so memory goes in Nextflow's own directive and `-profile slurm` needs
nothing special. `cpu-normal`, not `cpu-long`, is the default QoS to put in a site config:
Alpine's `cpu-long` rejects any job whose requested walltime is 24 hours or under ("Error 19"),
and `process_low`/`process_medium` ask for 1h/8h by default, so a config that defaults every task
to `cpu-long` fails at the first submission. `cpu-normal`'s own 24-hour cap still covers every
retry at the default settings (`maxRetries = 2`, so `process_medium` tops out at 8h × 3 = 24h);
a step that genuinely needs longer than a day should carry `--qos=cpu-long` on its own
clusterOptions rather than as the site-wide default. The head job above is unaffected, since its
own `--time=2-00:00:00` clears the floor.

Nextflow is a module on Alpine, but only where the job runs. A compute node's `MODULEPATH`
includes the software stack, so the `module load nextflow` in the script above finds
`nextflow/24.04.4` by default, with 22.10.6, 23.04 and 25.10.2 alongside it. That module brings
its own Java 21, so no `module load jdk` is needed beside it.

The login node is the exception, and it makes the cluster look at first as though it has no
Nextflow: its `MODULEPATH` holds only the Slurm modules, so `module avail nextflow` there returns
nothing. The stack is still on disk, and adding it by hand works, but the modulefile also reads
`$SCRATCHDIR`, which the login node does not set, and fails on a Lua error naming a line number
rather than the missing variable:

```bash
export SCRATCHDIR=/scratch/alpine          # set on compute nodes, not on the login node
module use /curc/sw/persistence-modules
module load nextflow
```

Since the head job is submitted with `sbatch` and runs on a compute node, none of that is needed
for a real run. It matters for a quick `-stub-run` on the login node, where
[installing Nextflow yourself](configuration.md#installing-it-yourself) against
`/curc/sw/jdk/18.0.1.1` is the shorter path: the login node's own `java` is 1.8, too old for
[any supported Nextflow](configuration.md#java-comes-first).

One thing the module changes quietly: it sets `NXF_HOME`, `NXF_TEMP` and `NXF_WORK` to
`$SCRATCHDIR/$USER`, so a run started under it puts its work directory on scratch rather than
in `./work`. That is where a work directory should be here, but it is not where `-resume` looks
if a later run does not load the same module.

### TACC Stampede3, reached through [ACCESS](https://allocations.access-ci.org/)

Stampede3 does not share nodes. TACC's user guide: *"TACC does not implement node-sharing on
any compute resource. Each Stampede3 compute node can be assigned to only one user at a time."*
Slurm is configured to match, `SelectType=select/linear` with `OverSubscribe=NO`. `sbatch
--test-only` costs no allocation and shows what that means for a job:

```
sbatch --test-only -N 1                  ->  48 processors on one node
sbatch --test-only -N 1 -c 4             ->  48 processors on one node
sbatch --test-only -N 1 --mem=8G         ->  48 processors on one node
sbatch --test-only -N 1 --mem=400G       ->  "Requested node configuration is not available"
```

Asking for four cores gets you forty-eight, and the bill is one node-hour either way. So the
arrangement every other cluster here uses, Nextflow's `slurm` executor submitting one job per
task, prices each task at a whole node. The nine-genome reference set is 102 tasks on Alpine;
run that way here, each task would bill a node for its wall time, and Stampede3's per-user cap
of 40 running jobs would serialise it besides.

The arrangement that fits is one batch allocation with Nextflow's `local` executor inside it,
scheduling against the node's own 48 cores. One node, held for the length of the run, 1 SU per
hour whatever the task count.

```bash
#!/bin/bash
#SBATCH --job-name=blimmp
#SBATCH --account=TG-BIO250056
#SBATCH --partition=skx
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=08:00:00
#SBATCH --output=nextflow-%j.out
#SBATCH --error=nextflow-%j.err

# No --mem and no --cpus-per-task. The node is yours entirely, and --mem is
# validated against its RealMemory while changing nothing about what you get.

module load tacc-apptainer          # compute nodes only; see below
export PATH=$WORK/nextflow-bin:$PATH

# Both caches on $WORK. $SCRATCH is purged, and a purged work directory takes
# -resume with it. The tacc-apptainer module points APPTAINER_CACHEDIR at
# $SCRATCH by default, so override it as well as the Nextflow one.
export NXF_SINGULARITY_CACHEDIR=$WORK/blimmp/singularity_cache
export APPTAINER_CACHEDIR=$WORK/blimmp/apptainer_cache

nextflow run main.nf \
    --genomes "$WORK/blimmp/genomes/*/*.fna" \
    --hmms    $WORK/blimmp/kofam \
    --sigma 1.0 \
    --hmm_chunks 24 \
    --outdir results_compare \
    -w $WORK/blimmp/work \
    -profile singularity \
    -c conf/stampede3.config
```

There is no `slurm` in the profile list and no `--max_concurrent`. On a local executor
`--max_concurrent` sets a queue size that nothing queues against, so it has no effect; the
pipeline sizes chunks from cores instead, and `--hmm_chunks` is the flag that does the work
here. Left to itself it cuts three chunks per genome, 27 tasks against 24 slots, which leaves
most of the node idle for the last wave. Twenty-four chunks is 216 tasks of about 1,160
profiles.

The nine-genome reference set reproduces `docs/reference_summary.tsv` exactly on one skx node,
in 9m26s over 264 tasks with nothing failed or retried. Every job behind this section started in
under a minute, on an allocation with no priority history, but that says more about the afternoon
than about `skx`: a later job sat pending behind 439 others. Size a walltime request for the run,
and expect to wait in the queue.

That first run held 16.3 of the node's 48 cores busy. Working out why led to a change to
`nextflow.config` that applies everywhere. `hmmsearch` parallelises over its target sequence
database, and the target is one genome's proteins. Timed on one chunk of 1159 profiles against
2565 proteins, alone on a node:

| `--cpu` | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| wall | 46.4s | 32.0s | 27.7s | 28.2s | 29.0s |
| total CPU | 48.8s | 49.2s | 49.8s | 50.8s | 53.3s |

Total CPU barely moves, so the added threads do no added work; the speedup saturates at 1.67x
and then reverses. `nextflow.config` now asks for two rather than four, which is 1.45x at 72%
efficiency instead of 1.67x at 42%, and that holds on every cluster. It shows up most here
because the local executor's slot count is cores divided by the search's cpus, so cores a
task does not use sit idle inside an allocation already paid for:

| | wall | CPU | cores busy | |
|---|---|---|---|---|
| `cpus = 4` | 9m26s | 2h40m | 16.3 of 48 | the old default |
| `cpus = 2` | 7m05s | 2h41m | 22.6 of 48 | the current one |
| `cpus = 1` | 4m54s | 2h49m | 34.6 of 48 | |

All three reproduce the reference table exactly, so the thread count does not move the answer.
Two is the pipeline default because it holds on clusters that hand out cores as well as on this
one. If you are holding a whole node anyway and want the run to finish sooner, `-c` a file
setting `withName: 'HMMSEARCH' { cpus = 1 }` over the top: nothing else on the node is waiting
for those cores. Every run behind this section cost under half an SU in total.

**Apptainer is a stub on the login nodes, and it exits 0.**

```
$ module load tacc-apptainer && apptainer pull docker://alpine:3.19
================================================================================
 Please do not run Apptainer on the login nodes! Instead, submit Apptainer
 job scripts to the queue with 'sbatch'.
================================================================================
$ echo $?
0
```

`/opt/apps/tacc-apptainer/*/bin/apptainer` on a login node is a shell script that prints that
notice and exits successfully, pulling nothing. Any pre-staging step that runs from a login node
and checks its exit status will report that it worked. Stage containers from inside a batch job
or an `idev` session instead. The compute nodes do reach the internet: `quay.io`, `ghcr.io`,
Docker Hub and `genome.jp` all answered from one, and the six images the pipeline needs pulled
on first use with no pre-staging.

Whole-node allocation has one more consequence. A library that finds
no thread limit sizes its pools from the machine's core count, and on a cluster that shares nodes
the cgroup makes that the right number by accident. Here it is 48. Nine concurrent BUSCO tasks each went
for 48 OpenBLAS threads, failed to allocate, and hung instead of exiting. All nine stuck, on an
idle node billing by the hour, with nothing in the Nextflow log to say so. TACC sets
`OMP_NUM_THREADS=1` in every batch job to prevent this, and `--cleanenv` discards it
before the container starts. `nextflow.config` now puts the limit on the container command line,
where `env -` cannot scrub it. With that, the nine-genome set with `--busco_lineage
bacteria_odb12` runs in 7m46s and BUSCO reports completeness between 0.922 and 0.966.

Nextflow has no module on Stampede3 and has to be installed, which is
[the ordinary path](configuration.md#installing-it-yourself); the system Java is already 17, so
nothing needs loading beside it. A login node will run `nextflow -version` and a `-stub-run`,
but only just: `ulimit -v` there is 8 GB against one visible core, and Nextflow's default thread
pool dies on it with `pthread_create failed (EAGAIN)` and a JVM that reports insufficient
memory. `NXF_OPTS='-Xms256m -Xmx1g -Xss512k -XX:ActiveProcessorCount=2'` gets a stub run through.
Inside a job, none of that applies.

### Reading a profile directory, and why it is threaded

Pointing `--hmms` at a directory makes `GROUP_HMMS` open every profile in it, and on a networked
filesystem almost all of that time is latency. The first run of the nine-genome set here spent
2h37m34s in that one task and 64 seconds of CPU, at 80 MB of memory, leaving the processor idle
for better than 99% of it. The task carries the `process_low` default of an hour, so it was
killed at the hour and retried with two.

`bin/group_hmms.py` now reads through a thread pool, which overlaps that latency. The same run
takes 2m08s, and the nine-genome set end to end went from 2h42m40s to 7m18s. Threads suffice
because the GIL is released during I/O, and `pool.map` preserves order, so a chunk holds exactly
what the serial version wrote, byte for byte.

Sixteen workers is the default because the curve flattens there. Measured on this filesystem
over disjoint cold slices of 2,000 profiles, one node, one moment:

| workers | 1 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|
| time | 97.1s | 19.6s | 11.3s | 9.3s | 8.3s |

Sixteen takes 89% of what 32 gets while asking half as much of a filesystem other people are
using. `BLIMMP_READ_WORKERS` overrides it, and 1 restores serial reading.

Treat the speedup as roughly tenfold rather than the seventyfold the two full runs suggest. A
serial read has one request outstanding and is fully exposed to how busy the filesystem is:
the same serial code measured 0.049 seconds per file during the sweep above and 0.32 during
that first run. Sixteen concurrent readers overlap those stalls and land near 0.0046 seconds
per file either way, which is why the threaded time is the reproducible half of the comparison.

A single concatenated `.hmm` is still worth keeping on shared storage for a database you will
search many times, since `SPLIT_HMMS` then streams one file and no run walks the directory at
all. With the reads threaded it is a convenience rather than a requirement:

```bash
find /path/to/kofam/profiles -name '*.hmm' -print0 | xargs -0 cat > /path/to/kofam/kofam.hmm
nextflow run main.nf --hmms /path/to/kofam/kofam.hmm ...
```

`xargs` rather than `cat profiles/*.hmm` because 27,864 paths is 2.2 MB of argument list against
a 2 MB `ARG_MAX`, and the shell reports that as `Argument list too long`. Do not reach for
`xargs -P` to speed this up: parallel `cat` processes writing to one redirected file interleave
their output, and the result differed from the correct file at byte 38 while still carrying the
right number of models, so the obvious check misses it.

### What differs, and what it costs to get wrong

| | UA HPC | Alpine | TACC Stampede3 |
|---|---|---|---|
| Partition | `standard` | `acpu` | `skx` (`skx-dev` for anything under 2h) |
| Account | `twheeler` | `ucb###_asc1` | `TG-BIO250056`, an ACCESS project code |
| Memory | `--mem-per-cpu` only; `--mem` is ignored | `--mem`, honoured per node | neither. Nodes are not shared, so the node's memory is already yours |
| QoS | not used | `cpu-normal` by default; `cpu-long` rejects anything 24h or under | not set by hand; one QoS per partition, applied automatically |
| Containers | Singularity, via the module system | Apptainer, already on every node | `module load tacc-apptainer`, **compute nodes only**; on a login node the binary is a stub that exits 0 |
| Nextflow | `module load nextflow/22.10.7.5853` | `module load nextflow` on a compute node; the login node sees no modules at all | no module at all; install it yourself. System Java is already 17 |
| Node sharing | yes | yes | **no.** One job holds one whole node, billed by the node-hour |
| Executor | `slurm`, one job per task | `slurm`, one job per task | **`local`, inside one batch allocation.** See below |
| Profile | `-profile uahpc,singularity` | `-profile slurm,singularity -c conf/alpine.config` | `-profile singularity -c conf/stampede3.config`, with no `slurm` |

Every row is a place a working script stops working. The memory row is the one that gives no
warning: the others produce an error at submission, while a memory flag the scheduler ignores
gives you jobs that die partway through a long run.

If your cluster matches neither, `conf/site_template.config` is the starting point and the rows
above are the questions to put to your own documentation. The Alpine values here come from
[CURC's job resources](https://curc.readthedocs.io/en/latest/running-jobs/job-resources.html),
[Alpine hardware](https://curc.readthedocs.io/en/latest/clusters/alpine/alpine-hardware.html)
and [containerization](https://curc.readthedocs.io/en/latest/software/containerization.html)
pages, which will outlast this table.

---

## BUSCO lineages, and why to download one first

`--busco_lineage` names a dataset BUSCO needs a copy of. Left to itself, BUSCO downloads that
copy when it runs, and every BUSCO task runs in its own work directory, so a nine-genome run
fetches the same dataset nine times. On a cluster whose compute nodes have no route off the
network, it fetches it zero times and the run stops.

`--busco_downloads` points at a directory you filled in beforehand, and the BUSCO step then runs
with `--offline`:

```bash
# Once, from a login node or anywhere with network access.
busco --download bacteria_odb12 --download_path /path/to/busco_downloads

nextflow run main.nf \
    --genomes '/path/to/genomes/*/*.fna' \
    --hmms /path/to/kofam/ \
    --busco_lineage bacteria_odb12 \
    --busco_downloads /path/to/busco_downloads \
    -profile uahpc,singularity
```

Two things to get right. Give an absolute path, because the directory is read from inside a
container that does not share your working directory. And check that the path is visible in
there: `singularity.autoMounts` covers `$HOME`, `/tmp` and the work directory, and anywhere else
needs a bind of your own.

`--busco_lineage auto_prok` (or `auto`, `auto_euk`) needs more in that directory than one
dataset. BUSCO first scores the genome against the domain root sets, then places it on a
reference tree, then scores it against whatever dataset the placement picks, and BUSCO only
learns which dataset that is partway through the run. So the directory has to hold the tree, the
roots, and every dataset it might pick, which means downloading the domain:

```bash
busco --download prokaryota --download_path /path/to/busco_downloads   # or eukaryota, or all
```

The pipeline checks for the tree and the root datasets before submitting anything, because
BUSCO 6.1.0 handles their absence badly. It looks for `lineages/bacteria_odb12.2` and its
siblings by exact name, because it also passes `--datasets_version odb12.2` to BUSCO. Left to its
default, that version follows the BUSCO release, so a container bump would change which datasets
auto-lineage chooses from without anything in the command line saying so. The pinned value is
`buscoDatasetsVersion()` in `main.nf`, and it moves with the container. With a root dataset missing it fails with an
`UnboundLocalError` from inside `AutoLineage.py` and never names a file. With the tree missing it
reports a failed download, in a mode where it was told not to download. The one thing the
pipeline cannot check is the dataset the placement will pick: if that one is missing, BUSCO does
say so, but only after the root runs and the placement, which on a 1.7 Mb Prochlorococcus genome
was most of the 81 seconds the whole task took. Download the domain and that case cannot arise.

```groovy
singularity.runOptions = '--cleanenv -B /path/to/busco_downloads'
```

`busco --list-datasets` prints the lineage names, which is worth running before you pick one.
The names carry the OrthoDB release: BUSCO 5 shipped `bacteria_odb10`, BUSCO 6 ships
`bacteria_odb12`, and a name from the wrong release is not found. Narrower lineages have far more
markers than broad ones, so they measure completeness more finely: `bacteria_odb12` has 116,
`cyanobacteriota_odb12` has 719, `prochlorococcus_odb12` has 1128.

## Parallelism and chunk size

Two axes multiply:

`PREDICT_PROTEINS` runs once per genome, and no genome needs anything from another, so they all
run at once.

The HMM database is split into chunks, and one `HMMSEARCH` job runs for every genome and chunk
pair. That product is where the job count comes from:

```
Total HMMSEARCH jobs = (number of genomes) × (number of HMM chunks)
```

**You do not have to pick this.** Left alone, the pipeline sizes the chunks itself and prints what
it chose:

```
HMM chunks: 6, about 4544 profiles each, from 27262 total (9 samples x 6 = 54 searches, filling 50 queue slots). Override with --hmm_chunk_size 4544 or --hmm_chunks 6.
```

It names the values it picked, so overriding is a matter of copying one of them and changing it
rather than working out what to pass. A run that sets either flag prints what that resolves to
instead, so a run always states its own chunking either way.

The rule is `chunks = slots / genomes`: enough tasks to keep the available slots busy, and no
more, because every task past that is another scheduler round trip buying no extra parallelism.
What counts as a slot is read from the resolved configuration rather than guessed, so it holds on
a site config of your own as much as on the profiles shipped here. On a scheduler it is the queue
size, which `--max_concurrent` sets. On a laptop nothing queues, so it is your cores divided by
what one search asks for, and the answer is usually one chunk.

Chunk size trades job count against scheduler overhead. Smaller chunks make more jobs that each
finish sooner, which spreads load across a large cluster. Larger chunks make fewer jobs and
waste less time on submission, which matters when the queue is slow to turn jobs around.

Pick it from your queue rather than from the database. On a busy cluster the wait dominates
everything else, and small chunks lose badly. Measured on nine genomes against all 27,262 KOfam
profiles: at `--hmm_chunk_size 100` each search ran for two to three seconds and then waited about
twenty minutes for a slot, and 2,457 jobs would have taken most of a day. The same work at
`--hmm_chunk_size 5000` is 54 jobs of a couple of minutes each and finishes in about ninety.
Total CPU is identical; what changes is how many queue waits you pay for it. `squeue` showing
`(Priority)` on everything is the signal to use larger chunks.

Either flag turns the automatic sizing off. Set profiles per chunk with `--hmm_chunk_size`, or fix
the number of chunks with `--hmm_chunks` when you want the job count to stay the same whatever the
database size. The two are mutually exclusive.

One floor applies even when the count is chosen for you: a chunk never holds fewer than 500
profiles. Nine genomes against all 27,262 profiles at 100 per chunk made 2,457 searches that each
computed for two to three seconds and then waited about twenty minutes for a slot on a busy
cluster. The same work in 54 larger tasks finished in about ninety minutes, for identical total
CPU. An explicit `--hmm_chunk_size` is yours and is not floored.

```bash
# Fixed chunk size of 100 profiles
nextflow run main.nf \
    --genomes data/genomes/ \
    --hmms data/kofam/profiles.hmm \
    --hmm_chunk_size 100

# Exactly 50 chunks regardless of HMM database size
nextflow run main.nf \
    --genomes data/genomes/ \
    --hmms data/kofam/profiles.hmm \
    --hmm_chunks 50
```

---

## Capping concurrency

Submitting a hundred thousand jobs at once will overwhelm a scheduler or trip a per-user limit.
`--max_concurrent` holds at most N jobs in the queue and submits more as they finish.

```bash
# Cap at 200 queued jobs
nextflow run main.nf \
    --genomes /path/to/genomes/ \
    --hmms /path/to/kofam/ \
    --max_concurrent 200 \
    -profile slurm,singularity \
    -c conf/mysite.config
```

Your cluster's per-user job limit is a good starting value; `sacctmgr show qos` or the local
documentation will have it. Too high trips the limit you were trying to avoid, too low leaves
cores idle.

---

## Testing at reduced scale

Before committing to a production run, put a handful of genomes through the whole pipeline on
the real cluster:

```bash
nextflow run main.nf \
    --genomes /path/to/genomes/ \
    --subsample_genomes 5 \
    --hmms /path/to/kofam/ \
    --subsample_hmms 50 \
    --max_concurrent 5 \
    -profile slurm,singularity \
    -c conf/mysite.config
```

`--subsample_genomes` takes the first N genomes from the input glob and `--subsample_hmms` the
first N profiles from the HMM database. Neither touches the other input.

This tests your cluster. To test the installation instead, `-profile test` runs bundled data in
about 30 seconds and asserts the result; see [before your first real
run](running.md#before-your-first-real-run).

---

## Reclaiming the space in work/

Everything `-resume` reuses lives in `work/`, so the two are one question. [Incremental
runs](running.md#incremental-runs) covers what gets reused.

Once a run has finished and published its outputs, `work/` can go, and `-resume` then has
nothing left to offer for that run. `nextflow clean -but <run-name> -f` is the narrower tool,
though it only removes bookkeeping for the runs you did not name, so it frees nothing when the
directories it points at were already deleted by hand.

`cleanup.sh`, in the repository root, removes the run artifacts from the current directory. It
lists what it found with sizes and waits for confirmation:

```bash
./cleanup.sh
```

| Path | What it is |
|---|---|
| `work/` | The task cache, and therefore the `-resume` cache |
| `.nextflow/` | Nextflow's own state: history and cache index |
| `.nextflow.log*` | Nextflow process logs |
| `nextflow-*.out/err` | Head job stdout and stderr |
| `results*/` | Every output directory |

It is worth running between test runs, and after a production run once you have the results
somewhere safe.

---

## Worked examples

```bash
# All genomes, all KOfam profiles, cap at 200 concurrent jobs
nextflow run main.nf \
    --genomes /path/to/all_genomes/ \
    --hmms /path/to/kofam/ \
    --max_concurrent 200 \
    -profile slurm,singularity \
    -c conf/mysite.config

# Same, but 500 chunks for finer parallelism:
# (num_genomes) x 500 searches, each of them very fast
nextflow run main.nf \
    --genomes /path/to/all_genomes/ \
    --hmms /path/to/kofam/ \
    --hmm_chunks 500 \
    --max_concurrent 500 \
    -profile slurm,singularity \
    -c conf/mysite.config

# Skip pyrodigal by supplying pre-computed protein FASTAs
nextflow run main.nf \
    --proteins /path/to/proteins/ \
    --hmms /path/to/kofam/ \
    --hmm_chunks 50 \
    --max_concurrent 100 \
    -profile slurm,singularity \
    -c conf/mysite.config
```
