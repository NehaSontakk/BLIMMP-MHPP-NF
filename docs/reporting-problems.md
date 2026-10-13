# When something goes wrong

Nextflow runs this pipeline on laptops, shared servers, several schedulers and several container
engines. We have tested the combinations we have access to, and there are many we do not.
A cluster can differ in ways nobody can predict from the outside: how it wants memory requested,
whether compute nodes reach the network, which module names exist, how Singularity is built.

So if this does not work on your system, open an issue at
[github.com/NehaSontakk/BLIMMP-MHPP-NF/issues](https://github.com/NehaSontakk/BLIMMP-MHPP-NF/issues)
and we will fix it. That is more useful to us than a workaround you keep to yourself.

---

## First, two commands that usually locate the problem

```bash
nextflow run main.nf -profile test,docker      # or test,singularity on a cluster
assets/test/assert.sh results_test
```

This runs the whole pipeline on 2.4 MB of bundled data in about 30 seconds, and its outcome
narrows things down a long way:

| outcome | what it means |
|---|---|
| passes | The installation is fine, so the problem is in your data, your parameters, or the scale of your run. Say so in the report; it rules out a great deal. |
| fails | The problem is in the environment: Nextflow, the container engine, the scheduler, or the filesystem. This is the more useful failure, because it is small and fast to reproduce. |
| will not start | Usually Nextflow or Java. See [getting Nextflow](configuration.md#getting-nextflow). |

---

## What to put in the report

Paste the output of this. It is everything we would otherwise have to ask for, one message
instead of four:

```bash
{
  nextflow -v
  java -version 2>&1 | head -1
  for engine in docker singularity apptainer; do
    command -v $engine >/dev/null && $engine --version
  done
  for sched in sbatch qsub bsub; do
    command -v $sched >/dev/null && echo "$sched: $($sched --version 2>&1 | head -1)"
  done
  uname -sr
  git rev-parse --short HEAD 2>/dev/null
} 2>&1
```

Run it from your clone, so the last line records which commit you are on.

Then, in your own words:

**The command you ran.** The whole thing, including `-profile` and any `-c` file. If you used a
site config of your own, paste it. It is the single most common place the answer is hiding.

**What you expected, and what happened instead.**

**The error Nextflow printed.** Not just the last line. The block that starts
`Error executing process` names the work directory, which is where the evidence is:

```
Work dir:
  /path/to/work/3a/f19c8b...
```

**The files from that work directory.** They are small:

| file | what it holds |
|---|---|
| `.command.sh` | the command that actually ran, after all interpolation |
| `.command.err` | what the tool printed to stderr, usually the real message |
| `.command.out` | stdout |
| `.command.log` | both, interleaved, plus the container invocation |
| `.exitcode` | the exit status |

**`.nextflow.log`** from the directory you launched from. Attach it as a file rather than
pasting it; the last few hundred lines are the useful part but the header records your whole
resolved configuration.

**If BLIMMP ran but the numbers look wrong**, add `<outdir>/<sample>/<sample>_debug.log`. It
records the completeness value and the taxonomic priors that run used, which is the
usual explanation for results that look off rather than broken.

---

## Things worth checking before you write

Not to discourage the report, only because these four account for most of them.

**Memory that is requested but not granted.** Some clusters ignore `--mem` and want
`--mem-per-cpu`; a job then gets a default allocation and is killed for exceeding it. Nothing
warns you. `scontrol show config | grep -i -E 'DefMemPer|MaxMemPer'` says which yours honours,
and [configuration](configuration.md#another-slurm-cluster) covers the fix.

**Compute nodes with no route to the internet.** The first task tries to pull a container image
and hangs or fails. Run `-profile test,singularity` once from a node that can reach the network,
with `NXF_SINGULARITY_CACHEDIR` set to shared storage; that fills the cache for every image the
pipeline uses.

**A run that stops making progress, on an idle node.** Nine concurrent BUSCO tasks on a 48-core
Stampede3 node hung indefinitely: no error, no failed task, Nextflow still alive, and the node at
load 0.31 with nothing running. The cause is in the container's stderr, which
Nextflow does not surface for a task that never exits:

```
OpenBLAS error: Memory allocation still failed after 10 retries, giving up.
```

OpenBLAS sizes its thread pool from the machine's core count unless a thread limit is set, so
nine tasks on a whole node asked for 48 threads each, failed to allocate, and parked in
`futex_wait_queue` instead of exiting. Sites that hand out whole nodes usually set
`OMP_NUM_THREADS` to stop this, and `--cleanenv` throws it away. The pipeline now passes the
limit on the container command line, which is the only place it survives; if you meet this with
another tool, `--env OPENBLAS_NUM_THREADS=1` is the shape of the fix.

One task alone finishes normally either way, so this needs concurrency to reproduce. Check
`.command.err` in a stuck task's work directory before anything else.

**A staging step on a login node that succeeded and staged nothing.** TACC Stampede3 puts a
stub at `/opt/apps/tacc-apptainer/*/bin/apptainer` on its login nodes: it prints "Please do not
run Apptainer on the login nodes", pulls nothing, and **exits 0**. Short of looking for the file
itself, a script has no way to detect the failed pull. Where a cluster confines
container work to compute nodes, check a pull by listing what it produced:

```bash
apptainer pull docker://alpine:3.19 && ls -la "$APPTAINER_CACHEDIR"
```

**A download that restarted instead of resuming.** `fetch_reference_data.sh` pulls 1.6 GB of
KOfam from a server slow enough that the transfer takes most of an hour, which is long enough to
be reset partway through. It resumes with `curl -C -` now, but if you see one begin again from
zero, check that before you conclude the network is blocked. A firewall fails in the first
second. A connection reset fails after transferring a gigabyte.

**A `-resume` that recomputed everything.** Bare `-resume` continues the *last* run from that
directory, which may not be the one you meant. `nextflow log` lists them, and
`-resume <run-name>` names one. See [which run does -resume continue?](running.md#which-run-does--resume-continue).

**A parameter that was accepted and ignored.** Nextflow options take one dash and pipeline
parameters take two. This pipeline checks its own spelling and will stop, but a parameter
introduced by a config of your own is not checked.

---

## If you are reporting a wrong answer rather than a crash

A run that finishes and produces nothing is the failure mode this pipeline has been bitten by
before, so it is worth separating from a crash. `assets/test/assert.sh` exists because a
BLIMMP release once shipped without its KEGG data files and produced output with the right
shape, the right columns and no content, for four months.

If your results differ from an earlier run rather than looking broken, check whether the two used
the same masking setting. Masking is on by default and `--no_masking` turns it off; the two
produce different scores, and thirteen module calls changed across nine genomes in our own
comparison. `.nextflow.log` records the resolved configuration for each run.

If your own results look empty or uniform, check the per-KO table first:

```bash
wc -l <outdir>/<sample>/<sample>_BLIMMP_substituted_dk.csv
```

One line means the header alone, and nothing scored. Include that number and the debug log.
