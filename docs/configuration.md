# Configuration

Most runs need one thing, a profile:

```bash
nextflow run main.nf --genomes ... --hmms ... -profile laptop,docker
```

A profile name has two halves. The first says where the work runs, the second says where the
tools come from. Both halves live in `nextflow.config`, so there is nothing to copy or edit
before a first run.

**Where the work runs**

| | |
|---|---|
| `standard` | this machine; what you get if you name no profile at all |
| `laptop` | this machine, 4 GB for bookkeeping steps and 8 GB for the search |
| `server` | this machine, 8 GB and 32 GB, and 16 CPUs for BUSCO and BLIMMP |
| `uahpc` | University of Arizona HPC (Puma / El Gato) |
| `slurm` | any other Slurm cluster, with a small site config naming your account |

**Where the tools come from**

| | |
|---|---|
| `docker` | pull each tool's container image |
| `singularity` | the same images under Singularity or Apptainer, which is what HPC clusters have |
| `modules` | `module load` HMMER, no containers |
| omit it | the tools are already on your `PATH` |

Take one from each list: `-profile server,docker`, or `-profile uahpc,singularity`.

A site config (`-c mysite.config`) comes into it only for a cluster with no profile of its own,
and then only to supply the things nobody can guess, an account and a partition. It is applied
after the profile, so anything in it wins.

To see what a combination resolves to before running anything:

```bash
nextflow config -profile uahpc,singularity .
```

---

## Getting Nextflow

This pipeline needs Nextflow 22.10 or later, and is tested on 22.10.7, 24.04.4, 25.10.4 and
26.04.0.

```bash
nextflow -version
```

If that prints a version at or above 22.10, you are done and can skip the rest of this section.

### Java comes first

Nextflow runs on the JVM, and the version it needs changed partway through its history:

| Nextflow | Java needed |
|---|---|
| 22.10 through 24.10 | 11 or later |
| 24.11 and later | 17 or later |

```bash
java -version
```

The Nextflow installer does not install Java. This catches people out, because the install
itself is one command and the Java prerequisite is invisible until it fails.

On a cluster, look for a Java module before installing anything:

```bash
module avail openjdk    # or: module avail java
module load openjdk/17.0.2
```

Off a cluster, [SDKMAN](https://sdkman.io/) is what the Nextflow documentation recommends:

```bash
curl -s https://get.sdkman.io | bash
sdk install java 17.0.10-tem
```

Check the Java on the machine that actually runs Nextflow. On a cluster that is a compute node,
not the login node, and the two often differ: on the system this pipeline was developed against,
the login node reports Java 8 while the compute nodes have Java 18.

### On a cluster, try the module system first

```bash
module avail nextflow
module load nextflow
```

Modules are the least trouble when they exist, because someone else has already dealt with the
Java pairing. Check what version you actually got, since `module load nextflow` gives whichever
one the site made the default, and that may be newer or older than you expect:

```bash
nextflow -version
```

### Installing it yourself

No root needed. Nextflow is a single self-contained launcher:

```bash
curl -s https://get.nextflow.io | bash    # creates ./nextflow here
mkdir -p ~/bin && mv nextflow ~/bin/      # somewhere on your PATH
```

On a cluster, put it on shared storage rather than in a node-local directory, and point
Nextflow's own working files somewhere shared as well:

```bash
export NXF_HOME=/path/to/shared/storage/.nextflow
```

### Running a specific version

`NXF_VER` picks the version for a single run, downloading it on first use. One launcher can run
any of them, so there is no need to install several:

```bash
NXF_VER=24.04.4 nextflow run main.nf --help
```

That is also the quickest way to test whether a problem is version-related, and the way to pin a
version for a run you want to reproduce later.

### If your cluster has an older Nextflow than you want

Install your own with the steps above; it takes precedence over the module if it comes earlier
on your `PATH`. Remember the Java pairing: a self-installed 25.x or 26.x needs Java 17, which may
itself have to come from a module.

---

## Installed tools or containers?

If HMMER and pyrodigal are already on the machine, or loadable with `module load`, using them
directly is simpler and skips pulling images. Leave the second half of the profile off, or use
`modules`.

Containers pin exact versions and need nothing installed, which is what lets a run reproduce on
a machine you have never touched. Docker suits laptops and servers. Singularity or Apptainer is
the HPC answer, because Docker wants root and clusters do not hand that out.

`-profile modules` covers HMMER and nothing else. pyrodigal, BUSCO and BLIMMP have no module
form here, so it describes a run that stops at the HMMER results. tantan has no module either,
and masking is on by default, so a run without containers needs `tantan` on your `PATH` or
`--no_masking` to skip the step. Masking is worth keeping if you can: see
[what the pipeline does](overview.md).

You will not discover this the hard way. A run that names no container profile checks for the
tools it is about to need and stops before submitting anything, listing whatever is missing.
`-stub-run` skips the check, since a stub invokes no tool.

---

## A laptop or workstation

```bash
nextflow run main.nf \
    --genomes /path/to/genomes/ \
    --hmms /path/to/kofam/profiles.hmm \
    -profile laptop,docker
```

`laptop` asks for 4 GB for the bookkeeping steps and 8 GB for hmmsearch, which is where the
time goes. If your machine has more to give, override it without editing anything shipped:

```groovy
// conf/mymachine.config
process {
    withLabel: 'process_medium' { memory = '24 GB' }
    withName:  'HMMSEARCH'      { cpus = 4 }
}
```

```bash
nextflow run main.nf ... -profile laptop,docker -c conf/mymachine.config
```

Neither number does much for the search.
A search holds about 45 MB resident, so the memory figure is headroom rather than a target. And
hmmsearch parallelises over its target sequence database, which here is one genome's proteins:
measured on a 48-core node, one chunk goes 46.4s on one core to 27.7s on four and no further.
What does move is how many searches run at once, and the local executor already fills the machine
with them.

HMMSEARCH also takes its cpus from a `withName` selector in `nextflow.config`, which outranks any
`withLabel`, so raising `process_medium.cpus` moves BUSCO and BLIMMP and leaves the search where
it was. Name it directly if you want to change it. The measurement is in a comment in
`nextflow.config`.

Drop `,docker` if HMMER and pyrodigal are already installed. BUSCO and BLIMMP are not covered
that way, so a run using `--busco_lineage` or producing module calls wants the containers.

`--busco_lineage` also reaches the network, once per genome, unless you hand it a lineage you
downloaded beforehand with `--busco_downloads`. On a cluster whose compute nodes cannot reach
the network, that is the difference between a run and a stopped run: see
[Running on HPC](hpc.md#busco-lineages-and-why-to-download-one-first).

---

## A large shared server

```bash
nextflow run main.nf \
    --genomes /path/to/genomes/ \
    --hmms /path/to/kofam/ \
    -profile server,docker
```

Same as above with more of everything: 32 GB, and 16 CPUs for BUSCO and BLIMMP. The search
itself uses two whatever the profile says, because that is where it stops gaining. Nextflow's
local executor will not start more tasks than the machine can hold, so the searches queue
themselves.

---

## University of Arizona HPC

```bash
export NXF_SINGULARITY_CACHEDIR=/xdisk/$USER/singularity_cache

nextflow run main.nf \
    --genomes /path/to/genomes/ \
    --hmms /path/to/kofam/ \
    -profile uahpc,singularity
```

The `uahpc` profile carries the account, the partition, and one local peculiarity: Slurm here
ignores `#SBATCH --mem` and silently allocates 1 MB, so memory has to arrive as `--mem-per-cpu`.
Nextflow renders its own `memory` directive as `--mem`, and sbatch refuses to see both, so the
profile clears `memory` and puts the request in `clusterOptions` instead. Nothing warns you when
this is wrong; jobs are killed for exceeding a limit you thought you had raised.

`-profile uahpc,modules` swaps the containers for `module load hmmer/3.4`, which suits an
annotation-only run. UA HPC has no pyrodigal module, so install it with `pip install pyrodigal`
first; the preflight check stops the run otherwise.

The head job has to outlive your login session. See [Running on HPC](hpc.md#submitting-the-head-job).

---

## Another Slurm cluster

Two things cannot be guessed, so they go in a file of your own:

```bash
cp conf/site_template.config conf/mysite.config
# fill in queue and clusterOptions
nextflow run main.nf ... -profile slurm,singularity -c conf/mysite.config
```

`conf/site_template.config` is commented through, and the partition and the account are the
only lines that must change.

Check which memory flag your cluster honours before committing to a long run:

```bash
scontrol show config | grep -i -E 'DefMemPer|MaxMemPer'
```

If it wants `--mem-per-cpu`, follow the `uahpc` profile: clear `memory` and pass the request
through `clusterOptions`, because sbatch rejects `--mem` and `--mem-per-cpu` together.

If both come back `UNLIMITED`, the cluster may not share nodes at all, and the answer is to ask
for no memory and no cpus: the whole node is yours, RAM included. `conf/stampede3.config`
is the worked example, and it does not use the `slurm` profile either, because one sbatch job
per task would bill a whole node per task.

Set `NXF_SINGULARITY_CACHEDIR` to somewhere shared and persistent so images are pulled once
rather than per run. If your compute nodes have no outbound network, run
`nextflow run main.nf -profile test,singularity` once from a node that does, which fills the
cache with every image the pipeline uses.

### Why `--cleanenv`

The `singularity` profile sets
`runOptions = '--cleanenv --env OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1'`. The point of a container is that
everyone runs the same software, but Singularity hands your shell settings to the container as
it starts, and some of them redirect a tool to code outside the image. A working copy of BLIMMP
on `PYTHONPATH` will be found first, and your results then come from whatever is in that
directory today rather than from the version you asked for. Nothing says so. `--cleanenv` starts
the container without your settings.

`--cleanenv` also drops any thread limit the site set for the job, which is why the `--env` half
is there. Without it OpenBLAS inside the container sizes its pool from the whole machine, and
concurrent BUSCO tasks on a large node hang rather than fail: measured on a 48-core node, none of
nine finished without the limit and all nine finished in about 20 seconds with it.

If you set `singularity.runOptions` in your own config, your value replaces this one rather than
joining it, so carry **both** halves across; keeping `--cleanenv` without the `--env` is the
combination that hangs.

---

## AWS Batch

Never run, so treat what follows as a starting point rather than a recipe. Batch needs the work
directory on S3, which is a command-line option rather than a config setting:

```groovy
// conf/mysite.config
process.executor  = 'awsbatch'
process.queue     = 'your-job-queue'
aws.region        = 'us-east-1'
aws.batch.cliPath = '/home/ec2-user/miniconda/bin/aws'
```

```bash
nextflow run main.nf ... -profile docker -c conf/mysite.config \
    -work-dir s3://your-bucket/work
```

---

## Overriding a profile

Anything passed with `-c` is applied after the profile and wins, so changing one number never
means editing a file that shipped with the pipeline:

```groovy
// a longer wall clock for the heavy steps
process {
    withLabel: 'process_medium' { time = '24h' }
}
```

`withLabel` reaches BUSCO and BLIMMP but not the search, whose cpus are pinned by name in
`nextflow.config`. To change that one, say `withName: 'HMMSEARCH' { cpus = 4 }`.

```bash
nextflow run main.nf ... -profile uahpc,singularity -c conf/more_cores.config
```

`nextflow config -profile uahpc,singularity -c conf/more_cores.config .` prints the result, which
is quicker than discovering it from a job that has already been submitted.
