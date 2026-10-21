#!/usr/bin/env nextflow

/*
 * hmmsearch_kofam: Parallel HMMER3 search of bacterial genomes against Kofam
 *
 * Usage:
 *   nextflow run main.nf --genomes /data/genomes/ --hmms /data/kofam/
 *   nextflow run main.nf --genomes 'data/genomes/*.fna' --hmms profiles.hmm
 *   nextflow run main.nf --proteins /data/proteins/ --hmms /data/kofam/
 *   nextflow run main.nf --genomes /data/genomes/ --hmms /data/kofam/ -profile slurm
 *   nextflow run main.nf --genomes /data/genomes/ --hmms /data/kofam/ --hmm_chunks 50
 *   nextflow run main.nf --genomes /data/genomes/ --hmms /data/kofam/ --max_concurrent 20
 *   nextflow run main.nf --genomes /data/genomes/ --hmms /data/kofam/ --busco_lineage bacteria_odb12
 *   nextflow run main.nf --genomes /data/bins/ --hmms /data/kofam/ --busco_lineage auto_prok
 *   nextflow run main.nf --proteins /data/proteins/ --hmms /data/kofam/ --sigma 0.95
 *
 * BLIMMP needs sigma, the fraction of each genome believed to be present. Either
 * --busco_lineage (measured per genome, needs --genomes) or --sigma (one value for
 * every input) provides it. With neither, the run annotates against Kofam, warns
 * that BLIMMP is being skipped, and stops there. --busco_lineage takes a dataset
 * name, or auto, auto_prok or auto_euk to let BUSCO choose one per genome.
 *
 * Re-running: add -resume to reuse completed tasks from work/. Nextflow decides what to
 * reuse by hashing each task's inputs, script, and container, so changing any of them
 * re-runs the affected tasks and leaves the rest cached.
 *
 * --genomes / --proteins: path to a directory of FASTA files, or a glob pattern
 * --hmms: path to a directory of .hmm files, or a single concatenated .hmm file
 *
 * Dry run (walks the whole graph, writing empty placeholder outputs):
 *   nextflow run main.nf --genomes 'data/genomes/*\/*.fna' --hmms /data/kofam/ -stub-run -profile standard
 *
 * Output is organized into one subdirectory per input directory:
 *   results/<strain>/<sample_id>.tblout
 *   results/<strain>/<sample_id>.domtblout
 *   results/summary.tsv
 */

nextflow.enable.dsl=2

// ---------------------------------------------------------------------------
// Processes
// ---------------------------------------------------------------------------

/*
 * Predict protein-coding genes with pyrodigal.
 * Skipped when the user supplies --proteins directly.
 */
process PREDICT_PROTEINS {
    tag "${strain_id}"
    label 'process_low'
    container 'quay.io/biocontainers/pyrodigal:3.7.1--py313hd72fa03_0'
    publishDir path: { "${params.outdir}/${strain_id}" }, mode: 'copy'

    input:
    tuple val(strain_id), path(genome)

    output:
    tuple val(strain_id), path("${genome.baseName}.faa")

    script:
    // pyrodigal has no -q: it is quiet already, and passing -q is an
    // "unrecognized arguments" error rather than a silent no-op.
    """
    pyrodigal -i ${genome} -a ${genome.baseName}.faa -o /dev/null
    """

    stub:
    """
    touch ${genome.baseName}.faa
    """
}

/*
 * Join the pieces of one sample into a single FASTA.
 *
 * Only runs for a samplesheet row naming more than one file, which is how a bin
 * or a draft assembly split across scaffold files gets treated as one sample.
 * Rows naming a single file bypass this entirely, so PREDICT_PROTEINS and
 * RUN_BUSCO still receive exactly one file and their task hashes are unchanged
 * for everyone who was here before samplesheets.
 *
 * Concatenating first is also the right call scientifically. pyrodigal trains its
 * gene model on whatever it is handed, so training once across a whole assembly
 * beats training separately on each contig file, and the gap widens as the
 * contigs get shorter.
 *
 * Borrows the pyrodigal container rather than pulling an image to run `cat`; the
 * output is on its way to pyrodigal regardless.
 */
process CONCAT_INPUTS {
    tag "${sample_id}"
    label 'process_low'
    container 'quay.io/biocontainers/pyrodigal:3.7.1--py313hd72fa03_0'

    input:
    tuple val(sample_id), path(pieces), val(extension)

    output:
    tuple val(sample_id), path("${sample_id}.${extension}")

    script:
    // Via a temporary name because a piece may already be called
    // <sample_id>.<extension>, and `cat x > x` truncates it before reading.
    """
    cat ${pieces} > concatenated.tmp
    mv concatenated.tmp ${sample_id}.${extension}
    """

    stub:
    """
    touch ${sample_id}.${extension}
    """
}

/*
 * Hard-mask low-complexity regions in the protein sequences with tantan.
 *
 * On by default; --no_masking turns it off. In its own process for the same
 * reason CONCAT_INPUTS is: masking rewrites every .faa, so every
 * HMMSEARCH task downstream gets a new hash, and a run made under the other
 * setting re-searches rather than resuming.
 *
 * The default was decided by measurement across nine Prochlorococcus strains on
 * 2026-08-20. Masking removed 209 hits, none of them above its KOfam threshold,
 * including hits scoring 7 against a threshold of 377 that BLIMMP was calling
 * present at confidence 1.0.
 *
 * It also recovers hits. KOfam scores many KOs on the best domain rather than
 * the full sequence, and a low-complexity stretch can fragment an alignment into
 * pieces that each score below the threshold while the full-sequence score sits
 * far above it. Masking the stretch lets HMMER align one long domain instead.
 * threonine dehydratase (K01754) went from a best domain of ~35 to ~485 against
 * a threshold of 261.5 in seven of the nine strains, and aspartate kinase
 * (K00928) from ~100 to ~460 against 398.1 in four.
 *
 * -p is protein mode and is not the default. Without it tantan reads a .faa as
 * DNA, warns on stderr, and writes a file anyway: lowercased, over different
 * regions. -x X substitutes a literal X, which HMMER scores as an unknown
 * residue.
 *
 * Expect the output rewrapped to 50 columns whatever the input wrapping was, and
 * masked spans running a few residues past the repeat itself, since tantan models
 * repeats probabilistically instead of matching character runs.
 */
process MASK_REPEATS {
    tag "${sample_id}"
    label 'process_low'
    container 'quay.io/biocontainers/tantan:51--h5ca1c30_1'

    input:
    tuple val(sample_id), path(proteins)

    output:
    tuple val(sample_id), path("${proteins.baseName}.masked.faa")

    script:
    """
    tantan -p -x X ${proteins} > ${proteins.baseName}.masked.faa
    """

    stub:
    """
    touch ${proteins.baseName}.masked.faa
    """
}

/*
 * Split a multi-model HMM file into chunks for parallel search.
 *
 * Two modes (mutually exclusive):
 *   --hmm_chunk_size N   → fixed number of HMM profiles per chunk (default 200)
 *   --hmm_chunks N       → fixed number of chunks (chunk size auto-calculated)
 */
process SPLIT_HMMS {
    label 'process_low'
    container 'python:3.11'

    input:
    path hmm_file
    val  auto_chunks

    output:
    path 'chunk_*.hmm'

    script:
    def chunks_arg     = params.hmm_chunks    ? "--chunks ${params.hmm_chunks}"         :
                         (auto_chunks         ? "--chunks ${auto_chunks}"               : '')
    def chunk_size_arg = params.hmm_chunk_size ? "--chunk-size ${params.hmm_chunk_size}" : ''
    def limit_arg      = params.subsample_hmms ? "--limit ${params.subsample_hmms}"      : ''
    """
    split_hmms.py ${chunks_arg} ${chunk_size_arg} ${limit_arg} ${hmm_file}
    """

    stub:
    """
    touch chunk_0000.hmm
    """
}

/*
 * Group individual HMM files from a directory into concatenated chunks.
 * Used when --hmms points to a directory of separate .hmm files.
 */
process GROUP_HMMS {
    label 'process_low'
    container 'python:3.11'

    input:
    val hmm_dir
    val auto_chunks

    output:
    path 'chunk_*.hmm'

    script:
    def chunks_arg     = params.hmm_chunks    ? "--chunks ${params.hmm_chunks}"         :
                         (auto_chunks         ? "--chunks ${auto_chunks}"               : '')
    def chunk_size_arg = params.hmm_chunk_size ? "--chunk-size ${params.hmm_chunk_size}" : ''
    def limit_arg      = params.subsample_hmms ? "--limit ${params.subsample_hmms}"      : ''
    """
    group_hmms.py ${chunks_arg} ${chunk_size_arg} ${limit_arg} ${hmm_dir}
    """

    stub:
    """
    touch chunk_0000.hmm
    """
}

/*
 * Search a chunk of HMM profiles against a set of protein sequences.
 * Produces per-sequence (tblout) and per-domain (domtblout) hit tables.
 */
process HMMSEARCH {
    tag "${sample_id}:${hmm_chunk.baseName}"
    label 'process_medium'
    container 'quay.io/biocontainers/hmmer:3.4--hb6cb901_4'

    input:
    tuple val(sample_id), path(proteins), path(hmm_chunk)

    output:
    tuple val(sample_id), path("${sample_id}.${hmm_chunk.baseName}.tblout"),    emit: tblout
    tuple val(sample_id), path("${sample_id}.${hmm_chunk.baseName}.domtblout"), emit: domtblout
    tuple val(sample_id), path("${sample_id}.${hmm_chunk.baseName}.out"),       emit: alignments

    script:
    // No --noali: that flag suppresses the per-domain alignment block in -o's
    // human-readable output, which is the whole point of keeping -o around.
    // tblout/domtblout stay the machine-readable record BLIMMP actually
    // parses; this file is for a human (or BLAST-style eyeballing) to read
    // the actual query/match/target alignment behind a given hit.
    """
    hmmsearch \
        --cpu ${task.cpus} \
        --tblout  ${sample_id}.${hmm_chunk.baseName}.tblout \
        --domtblout ${sample_id}.${hmm_chunk.baseName}.domtblout \
        -o ${sample_id}.${hmm_chunk.baseName}.out \
        ${hmm_chunk} ${proteins}
    """

    stub:
    """
    touch ${sample_id}.${hmm_chunk.baseName}.tblout
    touch ${sample_id}.${hmm_chunk.baseName}.domtblout
    touch ${sample_id}.${hmm_chunk.baseName}.out
    """
}

/*
 * Post-processing stub: run analyses across all merged results.
 * Receives all per-genome tblout files at once (collected channel).
 * Currently: count hits per genome at two E-value thresholds.
 * Add further analyses to postprocess.py as needed.
 */
process POSTPROCESS {
    label 'process_low'
    container 'python:3.11'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path tblouts   // all per-genome tblout files, staged into the work dir

    output:
    path 'summary.tsv'

    script:
    """
    postprocess.py *.tblout > summary.tsv
    """

    stub:
    """
    touch summary.tsv
    """
}

/*
 * Concatenate per-chunk results into a single tblout and domtblout per genome.
 * Preserves the three-line HMMER comment header from the first chunk file.
 */
process MERGE_RESULTS {
    tag "${sample_id}"
    label 'process_low'
    container 'ubuntu:22.04'
    publishDir path: { "${params.outdir}/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path('tblout/*'), path('domtblout/*'), path('alignments/*')

    output:
    tuple val(sample_id), path("${sample_id}.tblout"), path("${sample_id}.domtblout"), path("${sample_id}.out")

    script:
    """
    # --- tblout ---
    head -3 \$(ls tblout/* | sort | head -1) > ${sample_id}.tblout
    for f in tblout/*; do
        grep -v '^#' "\$f" >> ${sample_id}.tblout || true
    done

    # --- domtblout ---
    head -3 \$(ls domtblout/* | sort | head -1) > ${sample_id}.domtblout
    for f in domtblout/*; do
        grep -v '^#' "\$f" >> ${sample_id}.domtblout || true
    done

    # --- alignments (-o output) ---
    # Human-readable text, not tabular, so chunks are concatenated as-is
    # rather than header-stripped, each labelled by its source chunk.
    > ${sample_id}.out
    for f in alignments/*; do
        echo "==> \$(basename "\$f") <==" >> ${sample_id}.out
        cat "\$f" >> ${sample_id}.out
    done
    """

    stub:
    """
    touch ${sample_id}.tblout
    touch ${sample_id}.domtblout
    touch ${sample_id}.out
    """
}

/*
 * --busco_lineage names a dataset, or asks BUSCO to choose one. With auto,
 * BUSCO scores the genome against its domain-level root sets, places the
 * markers it found on a reference tree, and reruns against the deepest dataset
 * the placement supports. auto_prok and auto_euk skip the other domain's root
 * run. Returns the BUSCO flag for the auto values and null for a dataset name.
 */
def buscoAutoFlag(lineage) {
    ['auto': '--auto-lineage', 'auto_prok': '--auto-lineage-prok', 'auto_euk': '--auto-lineage-euk'][lineage]
}

/*
 * The dataset generation auto-lineage draws from. A named lineage carries its
 * own (bacteria_odb12); auto-lineage takes BUSCO's default, which follows the
 * BUSCO release (6.1.0 defaults to odb12.2). Pinning it here means a container
 * bump cannot silently change which datasets auto-lineage chooses from, and
 * lets checkBuscoAutoDownloads look for the exact directories BUSCO will open.
 * Bump it together with the busco container.
 */
def buscoDatasetsVersion() { 'odb12.2' }

/*
 * What an offline auto-lineage run needs in --busco_downloads: the placement
 * tree, the root dataset of every domain it may search, and whichever dataset
 * the placement then picks, which cannot be known ahead of the run. Checked
 * here because BUSCO 6.1.0 checks it badly: with a root dataset missing it dies
 * on an UnboundLocalError inside AutoLineage.py and never names a file, and
 * with the chosen dataset missing it does say so, but only after the root runs
 * and the placement, which is most of the work. This function cannot predict
 * the chosen dataset either, so the message says to download the whole domain,
 * which covers it.
 */
def checkBuscoAutoDownloads() {
    def dl = file(params.busco_downloads)
    def domain = ['auto': 'all', 'auto_prok': 'prokaryota', 'auto_euk': 'eukaryota'][params.busco_lineage]
    def roots = ['auto': ['archaea', 'bacteria', 'eukaryota'],
                 'auto_prok': ['archaea', 'bacteria'],
                 'auto_euk': ['eukaryota']][params.busco_lineage]

    // listFiles rather than files() with a glob, which matches files only.
    def placement = dl.resolve('placement_files')
    def lineages = dl.resolve('lineages')
    def missing = []
    if (!placement.isDirectory() || !placement.listFiles())
        missing << 'placement_files/'
    roots.each { root ->
        if (!lineages.resolve("${root}_${buscoDatasetsVersion()}").isDirectory())
            missing << "lineages/${root}_${buscoDatasetsVersion()}/"
    }
    if (!missing)
        return

    error """--busco_downloads ${params.busco_downloads} is missing what --busco_lineage ${params.busco_lineage} needs offline:
     |    ${missing.join('\n    ')}
     |  Auto-lineage places each genome on a reference tree and then scores it against the dataset
     |  the placement picks, so the directory has to hold the tree, the root datasets, and every
     |  dataset it might pick. Fill it on a machine with network access:
     |    busco --download ${domain} --download_path ${params.busco_downloads}
     |  The pipeline checks this because BUSCO reports a missing root dataset as a Python traceback.""".stripMargin()
}

/*
 * Assess completeness with BUSCO, on nucleotide genomes or on protein FASTAs.
 * Only runs when --busco_lineage is provided. busco_mode is 'genome' for
 * --genomes/genome_files rows and 'proteins' for --proteins/protein_files rows
 * (BUSCO's own -m vocabulary, used directly) - summary file naming and the
 * JSON schema BUSCO writes are identical between the two modes, verified
 * against a real protein FASTA, so the parsing below needs no branching.
 * Publishes the full BUSCO output directory and a per-genome TSV of
 * sample_id, completeness fraction, and the dataset that measured it.
 */
process RUN_BUSCO {
    tag "${strain_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/busco:6.1.0--pyhdfd78af_2'
    publishDir path: { "${params.outdir}/${strain_id}" }, mode: 'copy'

    input:
    tuple val(strain_id), path(genome), val(busco_mode)

    output:
    path "busco_${strain_id}/"
    tuple val(strain_id), path("${strain_id}.busco_completeness.tsv"), emit: completeness

    script:
    def lineage_arg = buscoAutoFlag(params.busco_lineage)
        ? "${buscoAutoFlag(params.busco_lineage)} --datasets_version ${buscoDatasetsVersion()}"
        : "-l ${params.busco_lineage}"
    // Without --busco_downloads, BUSCO fetches the lineage itself. Every task
    // does that separately, since each has its own work directory, so an N
    // genome run is N downloads of the same dataset and N chances for a compute
    // node with no outbound network to stop the run.
    def offline = params.busco_downloads
        ? "--offline --download_path ${file(params.busco_downloads).toAbsolutePath()}"
        : ''
    // BUSCO names the summary after the dataset it scored with, which under
    // auto-lineage is decided during the run, so the file is found rather than
    // named. Auto-lineage also leaves a short_summary.generic.* from the root
    // dataset beside it; the specific one is the result.
    """
    busco -i ${genome} -o busco_${strain_id} ${lineage_arg} -m ${busco_mode} -c ${task.cpus} ${offline}
    summary=\$(ls busco_${strain_id}/short_summary.specific.*.busco_${strain_id}.json)
    [ \$(printf '%s\\n' \${summary} | wc -l) -eq 1 ] || { echo "Expected one specific summary, found: \${summary}" >&2; exit 1; }
    python3 - "\${summary}" ${strain_id} > ${strain_id}.busco_completeness.tsv <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
print(sys.argv[2], round(d['results']['Complete percentage'] / 100, 4), d['lineage_dataset']['name'], sep='\\t')
PY
    """

    stub:
    """
    mkdir -p busco_${strain_id}
    printf "${strain_id}\t0.0\t${params.busco_lineage}\n" > ${strain_id}.busco_completeness.tsv
    """
}

/*
 * Collect per-genome BUSCO completeness values into a single summary table.
 * Columns: sample_id, completeness (0–1 fraction), busco_lineage (the dataset
 * that produced the number, which varies by row under --busco_lineage auto).
 */
process BUSCO_SUMMARY {
    label 'process_low'
    container 'ubuntu:22.04'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path completeness_files

    output:
    path 'busco_summary.tsv'

    script:
    """
    printf "sample_id\tcompleteness\tbusco_lineage\n" > busco_summary.tsv
    cat *.busco_completeness.tsv >> busco_summary.tsv
    """

    stub:
    """
    printf "sample_id\tcompleteness\tbusco_lineage\n" > busco_summary.tsv
    """
}

/*
 * Run BLIMMP module detection on a single genome.
 *
 * Takes the merged HMMER domtblout and sigma, the fraction of the genome
 * believed to be present. Sigma arrives as a number rather than as a BUSCO
 * file so that both sources look the same here: --busco_lineage measures it
 * per genome, --sigma states one value for every input.
 *
 * --taxonomy is passed through untouched. BLIMMP knows which prior sets it
 * ships and rejects a name it does not recognise with the valid ones listed
 * by rank, so repeating that list here would only give it somewhere to go
 * stale.
 */
process RUN_BLIMMP {
    tag "${sample_id}"
    label 'process_medium'
    // Built from github.com/NehaSontakk/BLIMMP_MHPP (Dockerfile in that repo).
    // Not yet published — build and push it before running this in anger.
    container 'ghcr.io/nehasontakk/blimmp_mhpp:26sep2026-mhpp1'
    publishDir path: { "${params.outdir}/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(domtblout), val(sigma), val(taxonomy)

    output:
    tuple val(sample_id), path("${sample_id}_BLIMMP_substituted_dk.csv"),      emit: dk
    path "${sample_id}_BLIMMP_substituted_modules.json",                      emit: modules
    path "${sample_id}_BLIMMP_substituted_module_steps.csv",                  emit: steps
    path "${sample_id}_BLIMMP_substituted_module_probabilities.csv",          emit: probabilities
    path "${sample_id}_BLIMMP_version.txt",                                  emit: version
    path "${sample_id}_debug.log", optional: true,                           emit: debug_log

    script:
    // BLIMMP (MHPP) runs its substitution pass by default and writes only the
    // *_BLIMMP_substituted_* files from it; the plain pass-1 files and the
    // debug log only appear with -v/--verbose, which writes both. So the old
    // -l/--no_blimmp_log split now rides on one flag: -v when the log is
    // wanted, nothing when it is not. Either way the substituted files -- the
    // ones this process reads -- get written.
    // The flag carries its own leading space, so that leaving it out gives
    // the same command string as before it existed rather than one with a
    // double space in it. The command string is part of the task hash, so a
    // stray space would re-run BLIMMP for everyone who never set the flag.
    def verbose_flag = params.no_blimmp_log ? '' : ' -v'
    // Arrives per sample so a samplesheet can score each row under its own
    // priors. Left to BLIMMP to validate, since BLIMMP holds the list of what it
    // ships and a copy here would only have somewhere to drift. Without it
    // BLIMMP uses domain-level priors.
    def tax_flag = taxonomy ? " -t ${taxonomy}" : ''
    def dest     = file("${params.outdir}/${sample_id}").toAbsolutePath()
    // Nextflow publishes nothing at all from a task that exits non-zero, so a
    // BLIMMP crash would leave its log in the work directory, which is where
    // nobody thinks to look. Copy it out on the way down instead. Best effort:
    // the copy needs ${params.outdir} to be visible from inside the container,
    // which held on this cluster but will not hold everywhere, so a failed copy
    // says so and leaves the real error to stand.
    def save_log = params.no_blimmp_log ? '' :
        "trap 'st=\$?; [ \$st -eq 0 ] || { mkdir -p ${dest} && cp -f ${sample_id}_debug.log ${dest}/ 2>/dev/null || echo \"[log] could not copy ${sample_id}_debug.log to ${dest}; it is in \$PWD\" >&2; }; exit \$st' EXIT"
    """
    ${save_log}
    BLIMMP ${domtblout} -f domtblout --sigma ${sigma}${tax_flag} --output ${sample_id}${verbose_flag}

    # BLIMMP (MHPP) isn't pip-packaged, so there is no dist metadata to ask for
    # a version; the image stamps one into /opt/blimmp/VERSION at build time.
    cat /opt/blimmp/VERSION > ${sample_id}_BLIMMP_version.txt 2>/dev/null \
        || echo unknown > ${sample_id}_BLIMMP_version.txt
    """

    stub:
    """
    touch ${sample_id}_BLIMMP_substituted_modules.json
    touch ${sample_id}_BLIMMP_substituted_dk.csv
    touch ${sample_id}_BLIMMP_substituted_module_steps.csv
    touch ${sample_id}_BLIMMP_substituted_module_probabilities.csv
    touch ${sample_id}_BLIMMP_version.txt
    """
}

/*
 * Decide whether a BLIMMP result is real.
 *
 * Separate from RUN_BLIMMP on purpose. A failing task publishes none of its
 * outputs, so making RUN_BLIMMP itself reject the result would strand the
 * debug log, the evidence, in the work directory. This way BLIMMP's own
 * success and our verdict on it are two different exit codes: RUN_BLIMMP
 * finishes and publishes, and this stops the run.
 */
process CHECK_BLIMMP {
    tag "${sample_id}"
    label 'process_low'
    // One failed comparison will fail again identically. Retrying wastes two
    // scheduler round trips and buries the message under NOTE lines.
    errorStrategy 'terminate'

    input:
    tuple val(sample_id), path(dk_csv)

    script:
    """
    # BLIMMP 0.1.2 shipped without its KEGG data files and produced output
    # that was structurally complete and entirely empty: correct columns,
    # 340 module rows, no annotation reaching the scoring step, exit 0.
    # Nine genomes of results looked fine for four months.
    #
    # The per-KO table is the direct read on whether anything scored. It had
    # 0 data rows when broken and 1980 when working. Do not test for a
    # non-zero module confidence instead: the broken output had 123 modules
    # above zero from the frequency prior alone.
    rows=\$(( \$(wc -l < ${dk_csv}) - 1 ))
    if [ "\${rows}" -lt 1 ]; then
        echo "[FATAL] ${sample_id}: BLIMMP wrote no per-KO rows to ${dk_csv}." >&2
        echo "        Every module would be reported absent regardless of the input." >&2
        echo "        Check that the BLIMMP container carries its KEGG data files." >&2
        echo "        The debug log published beside the results says what BLIMMP loaded." >&2
        exit 1
    fi
    """

    stub:
    """
    true
    """
}

/*
 * Pivot the per-genome BLIMMP module confidences into one table.
 *
 * Rows are KEGG modules, columns are genomes. This is the shape the results
 * are compared in: a module present in eight strains and absent in the ninth
 * is one row to read, where the per-genome files make it nine files to open.
 *
 * Runs once, after every genome has been scored, so it needs .collect().
 */
process MODULE_MATRIX {
    label 'process_low'
    container 'python:3.11'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path probability_files

    output:
    path 'module_matrix.tsv'

    script:
    """
    module_matrix.py *_BLIMMP_substituted_module_probabilities.csv > module_matrix.tsv
    """

    stub:
    """
    touch module_matrix.tsv
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/*
 * Pair each input file with its sample ID, failing if two inputs share one.
 *
 * A sample ID names the output subdirectory and prefixes every file in it, so
 * IDs must be unique. By default an ID is the name of the directory holding the
 * input, which suits the layout <genomes>/<strain>/<accession>.fna and yields
 * <strain>.
 *
 * A flat directory of genomes gives every input the same ID. That used to pass
 * silently: groupTuple() merged unrelated genomes into a single result and the
 * run reported one sample where there were many. Fail instead, and name the flag
 * that fixes it.
 */
def withSampleIds(List inputs, String kind) {
    def files = inputs
    if (!files)
        error "No ${kind} files matched. Check the path or glob you passed."

    def ids = files.collect { f -> params.id_from_filename ? f.baseName : f.parent.name }
    def duplicates = ids.countBy { it }.findAll { id, count -> count > 1 }

    if (duplicates) {
        def detail = duplicates.collect { id, count -> "    '${id}' from ${count} files" }.join('\n')
        def source = params.id_from_filename ? 'file names' : 'parent directory names'
        error """
        |Duplicate sample IDs derived from ${source}:
        |${detail}
        |
        |Sample IDs name the output subdirectory and prefix every file in it, so
        |results for these inputs would be merged into one.
        |
        |Either give each ${kind} its own directory, or ${params.id_from_filename
            ? 'rename the files so their names are unique'
            : 'add --id_from_filename to take IDs from file names instead'}.
        """.stripMargin()
    }

    return [files, ids].transpose().collect { f, id -> tuple(id, f) }
}

/*
 * Expand an input path into a file list.
 *
 * A bare directory is read as a flat directory of files with the given
 * extensions. Anything containing a '*' is used as written, which is how
 * nested layouts such as 'results/*\/*.domtblout' are given.
 *
 * files() rather than file(), which returns a bare Path for a spec holding no
 * glob and a list for one that does. Callers had to test the result's type
 * before they could count it. files() is a list in every case, and 26.04 warns
 * when file() is handed a glob.
 */
def expandInputs(String spec, String extensions) {
    def glob = (!spec.contains('*') && file(spec).isDirectory())
        ? "${spec}/*.${extensions}"
        : spec
    return files(glob)
}

/*
 * --help and the unknown-parameter check, called first thing by both entry
 * points.
 *
 * This lives in a function rather than at the top of the file because the strict
 * syntax parser, on by default from Nextflow 26.04, rejects statements that sit
 * outside a process, workflow, or function.
 */
/*
 * How many HMM chunks to cut the database into, when the user has not said.
 *
 * One HMMSEARCH task runs per genome and chunk, so the task count is
 * genomes x chunks. The number worth aiming for is however many tasks it takes
 * to keep the available slots busy, and no more: every task past that is
 * another scheduler round trip buying no extra parallelism.
 *
 *     chunks = slots / genomes
 *
 * What counts as a slot depends on where the run is happening, so it is read
 * from the resolved configuration rather than guessed from the profile name.
 * That way it holds for a site config of someone else's making as much as for
 * the profiles in this file. On a grid executor the cap is the queue size,
 * which --max_concurrent sets. On a local executor nothing queues and the limit
 * is cores divided by what one search asks for.
 *
 * Measured on this pipeline in August 2026: at 100 profiles a chunk, nine
 * genomes against the 27262 KOfam profiles made 2457 tasks that each computed
 * for two to three seconds and then waited about twenty minutes for a slot on a
 * busy cluster. The same work in 54 tasks finished in about ninety minutes.
 * Total CPU was identical. minProfilesPerChunk() exists so that a run with
 * many slots and one genome cannot rediscover that, since 27262 profiles across
 * 1000 slots would be 27 profiles a chunk.
 */
/*
 * Mirrors MIN_CHUNK_SIZE in bin/split_hmms.py and bin/group_hmms.py, so the
 * chunk size the log advertises is the one the scripts will use. A function
 * rather than a script-level constant: Nextflow 26.04's parser rejects a bare
 * assignment outside a process, workflow or function.
 */
def minProfilesPerChunk() { 500 }

def autoChunkCount(int genome_count) {
    def total = countProfiles()

    // The user has said. Report what that works out to anyway, so a run always
    // states its own chunking.
    if (params.hmm_chunk_size) {
        def n = total ? (Math.ceil(total / (params.hmm_chunk_size as int)) as int) : null
        log.info (n
            ? "HMM chunks: ${n}, ${params.hmm_chunk_size} profiles each, from ${total} total " +
              "(${genome_count} x ${n} = ${genome_count * n} searches). Set by --hmm_chunk_size."
            : "HMM chunks: ${params.hmm_chunk_size} profiles each, set by --hmm_chunk_size.")
        return 0
    }
    if (params.hmm_chunks) {
        def size = total ? (Math.ceil(total / (params.hmm_chunks as int)) as int) : null
        log.info (size
            ? "HMM chunks: ${params.hmm_chunks}, about ${size} profiles each, from ${total} " +
              "total (${genome_count} x ${params.hmm_chunks} = " +
              "${genome_count * (params.hmm_chunks as int)} searches). Set by --hmm_chunks."
            : "HMM chunks: ${params.hmm_chunks}, set by --hmm_chunks.")
        return 0
    }

    // workflow.session rather than a bare session: Nextflow 26.04's parser does
    // not define the latter. Checked on 22.10.7, 24.04.4, 25.10.4 and 26.04.0.
    def cfg      = workflow.session.config
    def executor = cfg.navigate('process.executor') ?: 'local'
    def slots
    if (executor == 'local') {
        // Nothing queues locally: Nextflow runs as many tasks as the cores
        // allow. So a slot is a core divided by what one search asks for, read
        // from whatever profile or -c file is in force rather than assumed.
        // HMMSEARCH's cpus is set by name in nextflow.config and a withName
        // selector outranks the label, so look there first.
        def per_task = (cfg.navigate('process.withName:HMMSEARCH.cpus')
                        ?: cfg.navigate('process.withLabel:process_medium.cpus')
                        ?: 4) as int
        slots = Math.max(1, Runtime.runtime.availableProcessors().intdiv(Math.max(1, per_task)))
    }
    else {
        slots = (cfg.navigate('executor.queueSize') ?: 100) as int
    }

    def chunks = Math.max(1, Math.ceil(slots / Math.max(1, genome_count)) as int)
    if (total) {
        // Mirrors the floor the chunking scripts apply, so the number reported
        // here is the number that will actually be used.
        def size = Math.ceil(total / chunks) as int
        if (size < minProfilesPerChunk()) {
            size   = Math.min(minProfilesPerChunk(), total)
            chunks = Math.ceil(total / size) as int
        }
        log.info "HMM chunks: ${chunks}, about ${size} profiles each, from ${total} total " +
                 "(${genome_count} ${genome_count == 1 ? 'sample' : 'samples'} x ${chunks} = " +
                 "${genome_count * chunks} searches, filling ${slots} " +
                 "${executor == 'local' ? 'local' : 'queue'} " +
                 "${slots == 1 ? 'slot' : 'slots'}). " +
                 "Override with --hmm_chunk_size ${size} or --hmm_chunks ${chunks}."
    }
    else {
        log.info "HMM chunks: ${chunks} (${genome_count} " +
                 "${genome_count == 1 ? 'sample' : 'samples'} x ${chunks} = " +
                 "${genome_count * chunks} searches, filling ${slots} " +
                 "${executor == 'local' ? 'local' : 'queue'} " +
                 "${slots == 1 ? 'slot' : 'slots'}). " +
                 "Override with --hmm_chunks ${chunks} or --hmm_chunk_size N."
    }
    return chunks
}

/*
 * How many profiles are on offer, when that is cheap to know.
 *
 * A directory of one profile per file is the KOfam layout and counting it is a
 * directory listing. A single concatenated .hmm has to be scanned to be
 * counted, and it can be a gigabyte, so that case reports nothing rather than
 * reading the file twice. The chunking scripts count it properly either way;
 * this is only so the log can name the size to override with.
 */
def countProfiles() {
    try {
        def input = file(params.hmms)
        if (!input.isDirectory())
            return 0
        def n = input.list().count { it.endsWith('.hmm') }
        return params.subsample_hmms ? Math.min(n, params.subsample_hmms as int) : n
    }
    catch (Exception e) {
        return 0
    }
}

def preflight() {
    if (params.help || params.help_all) {
        log.info paramsHelp(params.help_all as boolean)
        exit 0
    }
    if (params.validate_params)
        validateParamNames()
    checkCountParams()
    checkToolsOnPath()
}

/*
 * Numbers that have to be whole and positive.
 *
 * These arrive as strings, because a parameter declared null in nextflow.config
 * gives Nextflow no type to coerce to. A string is truthy whatever it says, so
 * `--hmm_chunk_size 0` sails past a plain `if (params.hmm_chunk_size)` and turns
 * into a division by zero further down, and `--hmm_chunks abc` into a
 * NumberFormatException. Both deserve a sentence rather than a stack trace.
 */
def checkCountParams() {
    ['hmm_chunk_size', 'hmm_chunks', 'max_concurrent'].each { name ->
        def value = params[name]
        if (value == null || value == '')
            return
        def text = value.toString()
        if (!text.isInteger() || (text as int) < 1)
            error "--${name} takes a whole number of 1 or more, not '${text}'."
    }
    ['subsample_genomes', 'subsample_hmms'].each { name ->
        def text = params[name].toString()
        if (!text.isInteger() || (text as int) < 0)
            error "--${name} takes a whole number of 0 or more, not '${text}'."
    }
}

/*
 * Refuse a run whose tools are not there, before any task is submitted.
 *
 * With -profile docker or singularity every tool arrives in its own image and
 * there is nothing to check. Without one, the tools have to be on PATH already,
 * and the failure otherwise is a "command not found" inside one task among
 * thousands, hours after submission.
 *
 * -profile modules covers hmmsearch through beforeScript, so it is not checked
 * under it. Nothing supplies pyrodigal, tantan, BUSCO or BLIMMP that way: no
 * module exists for any of them here.
 *
 * Only the tools a given run reaches are checked. Demanding one that the run
 * never invokes turns a working setup into a refused run, which is the failure
 * this function exists to prevent.
 */
def checkToolsOnPath() {
    // A stub run never invokes a tool, so requiring them would make -stub-run
    // useless on exactly the machine where you want to check the wiring first.
    if (workflow.containerEngine || workflow.stubRun)
        return

    def modules = workflow.profile.tokenize(',').contains('modules')
    def wanted = [:]
    // --proteins skips PREDICT_PROTEINS entirely, and a samplesheet with no
    // genome_files cell does the same, so neither needs a gene caller.
    // CONCAT_INPUTS borrows the pyrodigal image but only runs cat, which needs
    // nothing on PATH.
    if (!params.proteins && !(params.input && !sheetHasGenomeRows(params.input)))
        wanted['pyrodigal'] = 'gene calling'
    if (!modules)
        wanted['hmmsearch'] = 'the KOfam search'
    if (!params.no_masking)
        wanted['tantan'] = 'protein masking, which --no_masking turns off'
    if (params.busco_lineage)
        wanted['busco'] = '--busco_lineage'
    if (params.sigma != null || params.busco_lineage || params.input || params.completeness)
        wanted['BLIMMP'] = 'module detection'

    def missing = wanted.findAll { tool, _why ->
        ["bash", "-c", "command -v ${tool}"].execute().waitFor() != 0
    }
    if (!missing)
        return

    error """Not on PATH: ${missing.keySet().join(', ')}
     |
     |${missing.collect { tool, why -> "  ${tool.padRight(10)} needed for ${why}" }.join('\n')}
     |
     |  This run names no container profile, so every tool has to be installed
     |  already. Add -profile docker, or -profile singularity on a cluster, and
     |  they are pulled for you.""".stripMargin()
}

/*
 * nextflow_schema.json describes every parameter once. --help prints from it
 * and validateParamNames checks against it, so a parameter cannot be documented
 * in one place and validated against another, and adding one means editing a
 * single file.
 *
 * The file follows the nf-core schema format, so the nf-schema plugin can take
 * over this job unchanged if the Nextflow floor ever rises far enough to use it.
 * Reading it here instead keeps the pipeline working on a cluster with no
 * outbound network, where a plugin fetched at runtime would stop the run dead.
 */
def readSchema() {
    def schema_file = file("${projectDir}/nextflow_schema.json")
    if (!schema_file.exists())
        error "nextflow_schema.json not found in ${projectDir}. --help and the parameter check both read it."
    return new groovy.json.JsonSlurper().parseText(schema_file.text)
}

/* Wrap prose to a width and indent it, so help output stays readable at 80 columns. */
def wrapText(String text, int width, String indent) {
    def lines = []
    def line = new StringBuilder()
    text.split(/\s+/).each { word ->
        if (line.length() > 0 && line.length() + 1 + word.length() > width) {
            lines << line.toString()
            line = new StringBuilder()
        }
        if (line.length() > 0) line << ' '
        line << word
    }
    if (line.length() > 0) lines << line.toString()
    return lines.collect { "${indent}${it}" }.join('\n')
}

def paramsHelp(boolean full) {
    def schema = readSchema()
    def out = new StringBuilder()
    out << "\n${schema.title}\n\n"
    out << wrapText(schema.description, 78, '') << "\n"
    out << """
Usage
  nextflow run main.nf --input samples.csv --hmms /path/to/kofam/
  nextflow run main.nf --genomes 'genomes/*/*.fna' --hmms /path/to/kofam/ --sigma 0.95
  nextflow run main.nf --genomes /data/genomes/ --hmms /path/to/kofam/ --busco_lineage bacteria_odb12
  nextflow run main.nf -entry blimmp_only --domtblout 'results/*/*.domtblout' --sigma 1.0

Single-dash options belong to Nextflow: -profile, -resume, -c, -work-dir,
-entry. Everything below takes two dashes. Getting that backwards fails
without complaint, so check the dashes first when a flag seems to do nothing.
"""
    schema.allOf.each { ref ->
        def group = schema.definitions[ref['$ref'].tokenize('/').last()]
        out << "\n${group.title}\n"
        if (group.description)
            out << wrapText(group.description, 76, '  ') << "\n"
        group.properties.each { name, spec ->
            // The schema's `format` says what kind of string it is, and "<file>"
            // tells a reader more than "<string>" does.
            def label = [
                'file-path'     : 'file',
                'directory-path': 'directory',
                'path'          : 'path',
            ].get(spec.format, spec.type)
            def kind = spec.type == 'boolean' ? '' : " <${label}>"
            def has_default = spec.containsKey('default') &&
                !(spec.default instanceof Boolean && spec.default == false)
            out << "\n  --${name}${kind}${has_default ? "   [default: ${spec.default}]" : ''}\n"
            out << wrapText(spec.description, 72, '      ') << "\n"
            if (full && spec.help_text)
                out << wrapText(spec.help_text, 72, '      ') << "\n"
        }
    }
    out << (full ? "\n" : "\nAdd --help_all for the full explanation of each parameter.\n")
    return out.toString()
}

def knownParamNames() {
    def schema = readSchema()
    def names = [] as Set
    schema.allOf.each { ref ->
        names.addAll(schema.definitions[ref['$ref'].tokenize('/').last()].properties.keySet())
    }
    return names
}

/*
 * Which taxonomic priors a sample is scored against: its samplesheet cell if
 * there is a sheet, otherwise the run-level --taxonomy. Takes the lookup as an
 * argument because the strict syntax parser will not call a closure held in a
 * local variable.
 */
def taxonomyFor(String id, Map sheet_taxonomy) {
    return params.input ? sheet_taxonomy[id] : params.taxonomy
}

/*
 * camelCase to hyphenated. A real function rather than a closure in a local
 * variable, which the strict syntax parser will not let you call like one.
 */
def toKebab(String name) {
    return name.replaceAll(/([A-Z])/) { match -> "-${match[1].toLowerCase()}" }
}

/* Levenshtein distance, used only to guess what a misspelled parameter meant. */
def editDistance(String a, String b) {
    def previous = (0..b.length()).toList()
    a.toCharArray().eachWithIndex { ch, i ->
        def current = [i + 1]
        b.toCharArray().eachWithIndex { other, j ->
            current << [previous[j + 1] + 1,
                        current[j] + 1,
                        previous[j] + (ch == other ? 0 : 1)].min()
        }
        previous = current
    }
    return previous[b.length()]
}

/*
 * Stop on a parameter this pipeline does not read.
 *
 * Nextflow accepts every --name it is given, because a pipeline may invent any
 * parameter it likes and there is no list to check against. So a misspelling is
 * set, read by nothing, and the run carries on: --outdr writes results to
 * 'results' while you wait for them somewhere else, and --resume runs the whole
 * thing again from scratch. The schema gives us the list Nextflow lacks.
 */
def validateParamNames() {
    def known = knownParamNames()
    def unknown = params.keySet().findAll { !known.contains(it) }
    if (!unknown) return

    // Nextflow's own options, which are the likeliest thing to arrive here.
    // Writing one with two dashes is the most expensive typo available:
    // --resume runs the whole pipeline again from scratch, and --profile leaves
    // every task on the login node instead of going to the scheduler.
    def nextflow_options = ['resume', 'profile', 'config', 'work-dir', 'entry', 'params-file',
                            'revision', 'latest', 'offline', 'preview', 'stub-run', 'bg',
                            'log', 'name', 'with-report', 'with-trace', 'with-timeline',
                            'with-dag', 'dump-hashes', 'ansi-log'] as Set

    // Nextflow adds a camelCase twin for any parameter written with hyphens, so
    // --work-dir arrives as both 'work-dir' and 'workDir'. Report the mistake
    // once, in the hyphenated spelling.
    def reported = unknown.collect { name -> toKebab(name) }.unique().sort()

    def lines = reported.collect { name ->
        if (nextflow_options.contains(name))
            return "  --${name}   belongs to Nextflow, not to this pipeline. Write it with one dash: -${name}"
        def closest = known.min { editDistance(it, name) }
        return (closest != null && editDistance(closest, name) <= 3)
            ? "  --${name}   did you mean --${closest}?"
            : "  --${name}"
    }
    error """Unrecognised parameter${reported.size() > 1 ? 's' : ''}:
     |
     |${lines.join('\n')}
     |
     |Nextflow accepts any --name you type, so this would otherwise be set as a
     |parameter nothing reads, and the run would carry on without it.
     |
     |--help lists every parameter this pipeline takes. If a config of your own
     |declares extra parameters, --validate_params false turns this check off.""".stripMargin()
}

/*
 * Check one sigma value and return it as a number.
 *
 * Sigma reaches the pipeline three ways: the --sigma flag, a samplesheet column,
 * and the --completeness table. `source` names whichever one is being checked so
 * the message points at the thing the user has to edit.
 */
def checkSigma(value, String source) {
    if (!("${value}" ==~ /^[0-9]*\.?[0-9]+$/))
        error "${source} must be a number. Got: ${value}"
    def number = "${value}" as Double
    if (number <= 0 || number > 1)
        error """${source} must be greater than 0 and no greater than 1. Got: ${value}
         |  It is a fraction, so a genome estimated 85% complete is 0.85, and a
         |  sequence you believe is complete is 1.0.""".stripMargin()
    return number
}

/*
 * Check --sigma and return it as a number, or null when it was not given.
 * Shared by both entry points so the rule cannot drift between them.
 */
def validatedSigma() {
    return params.sigma == null ? null : checkSigma(params.sigma, '--sigma')
}

/*
 * Read a table of completeness values, keyed by sample ID.
 *
 * Written for the busco_summary.tsv this pipeline produces, which carries a
 * header, and equally happy with the headerless per-genome TSVs: any row whose
 * second field is not a number is skipped, which covers the header without
 * needing to know what it says. Columns past the second are ignored, so the
 * busco_lineage column and any file that predates it both read the same.
 */
def readCompleteness(String path) {
    def table_file = file(path)
    if (!table_file.exists())
        error "--completeness file not found: ${path}"

    def table = [:]
    table_file.eachLine { line ->
        def parts = line.trim().split(/\s+/)
        if (parts.size() < 2 || !(parts[1] ==~ /^[0-9]*\.?[0-9]+$/))
            return
        def value = parts[1] as Double
        if (value <= 0 || value > 1)
            error """--completeness ${path}: '${parts[0]}' has completeness ${parts[1]}.
             |  Values are fractions between 0 and 1, so 85% complete is 0.85.""".stripMargin()
        table[parts[0]] = value
    }

    if (!table)
        error """--completeness ${path} has no usable rows.
         |  Expected a sample ID and a completeness fraction as the first two columns,
         |  as written to busco_summary.tsv by a run using --busco_lineage.""".stripMargin()
    return table
}

/*
 * Turn one samplesheet cell into the list of files it names.
 *
 * A cell holds a single path or a glob, never a bare directory. A glob states
 * which files it means, where a directory leaves the reader guessing whether the
 * extra files in it were meant to be included.
 *
 * Relative paths resolve against the samplesheet's own directory rather than the
 * launch directory Nextflow would use, so a sheet and its data can be moved or
 * copied together and keep working.
 */
def resolveInputCell(String cell, base, String column, String name, String path) {
    def resolved = cell.startsWith('/') ? cell : "${base}/${cell}"

    if (!cell.contains('*') && file(resolved).isDirectory())
        error """--input ${path}: ${column} for '${name}' is a directory: ${cell}
         |  Name the files with a glob instead, such as ${cell}/*.fna, so the sheet
         |  says which of the files in there belong to this sample.""".stripMargin()

    def matched = files(resolved).findAll { it.exists() }
    if (!matched)
        error """--input ${path}: ${column} for '${name}' matched no files: ${cell}
         |  Relative paths are read from the samplesheet's own directory, which is
         |  ${base}.""".stripMargin()

    // Sorted because a glob's order decides concatenation order, which decides
    // the task hash. Unsorted, the same sheet could re-run work on a filesystem
    // that returned matches in a different order.
    return matched.sort()
}

/*
 * Read a samplesheet into one map per row.
 *
 * Columns are `name`, one of `genome_files` or `protein_files`, and the optional
 * `sigma` and `taxonomy`. A blank taxonomy falls back to --taxonomy here. A blank
 * sigma stays null and gets resolved in the workflow, where BUSCO's measurement
 * sits between the cell and the --sigma fallback; collapsing it here would let
 * --sigma outrank a measurement it should lose to.
 *
 * Everything is checked here, before a single task is submitted, and every
 * message names the line and the sample it came from.
 */
/*
 * Does any row of a samplesheet start from nucleotides?
 *
 * Only the preflight tool check asks, and only to decide whether this run
 * reaches a gene caller. It reads cells and nothing else: no glob resolution,
 * no column validation, no file existence checks. readSamplesheet does all of
 * that later and owns every error message about a sheet, so anything unreadable
 * here answers true and leaves the real complaint to the real parser. Guessing
 * true costs a run that names a gene caller it did not need; guessing false
 * would let one through to fail inside a task.
 */
def sheetHasGenomeRows(String path) {
    try {
        def sheet = file(path)
        if (!sheet.exists())
            return true
        def sep  = sheet.name.toLowerCase().endsWith('.tsv') ? '\t' : ','
        def rows = sheet.splitCsv(header: true, sep: sep, strip: true, quote: '"')
        if (!rows)
            return true
        if (!rows[0].keySet().contains('genome_files'))
            return false
        return rows.any { "${it['genome_files'] ?: ''}".trim() }
    }
    catch (Exception e) {
        return true
    }
}

def readSamplesheet(String path) {
    def sheet = file(path)
    if (!sheet.exists())
        error "--input samplesheet not found: ${path}"

    def base = sheet.toAbsolutePath().parent
    def sep  = sheet.name.toLowerCase().endsWith('.tsv') ? '\t' : ','
    // `quote` has to be named explicitly. Nextflow's splitCsv defaults to no
    // quote character at all, so without this a quoted field holding a comma is
    // cut in half and keeps its quote marks. A brace glob such as *.{fa,fna} is
    // exactly that case, and it arrived as '"..../*.{fa' before this was set.
    def rows = sheet.splitCsv(header: true, sep: sep, strip: true, quote: '"')
    if (!rows)
        error "--input ${path} has no rows below its header."

    def columns = rows[0].keySet()
    if (!columns.contains('name'))
        error """--input ${path} has no 'name' column. Found: ${columns.join(', ')}
         |  The first line is a header naming the columns, and 'name' is the one
         |  that gives each row its sample ID.""".stripMargin()

    def known = ['name', 'genome_files', 'protein_files', 'sigma', 'taxonomy']
    def unknown = columns.findAll { !(it in known) }
    if (unknown)
        error """--input ${path} has columns this pipeline does not read: ${unknown.join(', ')}
         |  Known columns are ${known.join(', ')}. A misspelled column would
         |  otherwise be ignored without comment, so it stops the run instead.""".stripMargin()

    def samples = []
    rows.eachWithIndex { row, index ->
        def line = index + 2   // the header is line 1
        def name = "${row['name'] ?: ''}".trim()
        if (!name)
            error "--input ${path} line ${line}: the name column is empty."

        def genome_cell  = "${row['genome_files'] ?: ''}".trim()
        def protein_cell = "${row['protein_files'] ?: ''}".trim()
        if (genome_cell && protein_cell)
            error """--input ${path} line ${line}: '${name}' fills in both genome_files and protein_files.
             |  A sample starts either from nucleotides, which get gene-called, or
             |  from proteins, which skip pyrodigal. Give one of the two.""".stripMargin()
        if (!genome_cell && !protein_cell)
            error """--input ${path} line ${line}: '${name}' names no input files.
             |  Fill in genome_files with nucleotide FASTAs, or protein_files with
             |  protein FASTAs to skip gene calling.""".stripMargin()

        def from_proteins = protein_cell as Boolean
        def column        = from_proteins ? 'protein_files' : 'genome_files'
        def sigma_cell    = "${row['sigma'] ?: ''}".trim()
        def taxonomy_cell = "${row['taxonomy'] ?: ''}".trim()

        samples << [
            name        : name,
            files       : resolveInputCell(genome_cell ?: protein_cell, base, column, name, path),
            from_proteins: from_proteins,
            sigma       : sigma_cell
                ? checkSigma(sigma_cell, "--input ${path} line ${line}: sigma for '${name}'")
                : null,
            taxonomy    : taxonomy_cell ?: params.taxonomy,
        ]
    }

    def duplicates = samples.countBy { it.name }.findAll { _name, count -> count > 1 }
    if (duplicates)
        error """--input ${path} uses the same name on more than one row: ${duplicates.keySet().sort().join(', ')}
         |  A name gives the output subdirectory and prefixes every file in it, so
         |  two rows sharing one would write over each other.""".stripMargin()

    return samples
}

workflow {
    preflight()

    // ---- Input validation ----
    if (params.input && (params.genomes || params.proteins))
        error """Provide --input or --genomes/--proteins, not both.
         |  A samplesheet already names every input file, one row per sample.""".stripMargin()
    if (!params.input && !params.genomes && !params.proteins)
        error "Provide --input (a samplesheet), --genomes (nucleotide FASTAs), or --proteins (protein FASTAs)"
    if (params.genomes && params.proteins)
        error "Provide --genomes or --proteins, not both"
    if (!params.hmms)
        error "Provide --hmms (path to Kofam profiles.hmm)"
    if (params.hmm_chunk_size && params.hmm_chunks)
        error "Provide --hmm_chunk_size or --hmm_chunks, not both"
    if (params.busco_downloads && !params.busco_lineage)
        error "--busco_downloads only does something with --busco_lineage, which is what runs BUSCO"
    if (params.busco_downloads && !file(params.busco_downloads).isDirectory()) {
        // A named lineage is one download. Auto-lineage needs the tree and
        // every dataset it might choose, which is the domain.
        def download = ['auto': 'all', 'auto_prok': 'prokaryota', 'auto_euk': 'eukaryota'][params.busco_lineage]
            ?: params.busco_lineage
        error file(params.busco_downloads).exists()
            ? "--busco_downloads is a file, and has to be a directory: ${params.busco_downloads}"
            : """--busco_downloads does not exist: ${params.busco_downloads}
             |  Populate it first, on a machine with network access:
             |    busco --download ${download} --download_path ${params.busco_downloads}
             |  Then every BUSCO task reads it instead of downloading its own copy.""".stripMargin()
    }
    if (params.busco_downloads && buscoAutoFlag(params.busco_lineage))
        checkBuscoAutoDownloads()

    // Sigma is the fraction of the genome believed to be present, and BLIMMP
    // needs one to weigh a missing KO against an incomplete assembly.
    //
    // On the command line both --sigma and --busco_lineage describe the whole
    // run, and nothing says which should win, so having both is refused. A
    // samplesheet does establish an order: a number typed on one genome's row is
    // a more specific claim than a number typed for the run, so the row wins,
    // then BUSCO's measurement, then --sigma filling in whatever is left blank.
    if (params.sigma != null && params.busco_lineage && !params.input)
        error """Provide --sigma or --busco_lineage, not both.
         |  BUSCO measures completeness separately for each genome; --sigma applies one
         |  value to every input. There is no sensible way to combine them.
         |
         |  With --input they do combine: the sigma column wins where it is filled in,
         |  BUSCO's measurement covers the rest, and --sigma is the last fallback.""".stripMargin()

    def sigma_value = validatedSigma()
    def samples     = params.input ? readSamplesheet(params.input) : null

    // ---- Protein sequences ----
    //
    // Three ways in, funnelling into two channels: nucleotides that need gene
    // calling, and proteins that do not. PREDICT_PROTEINS is invoked once below,
    // on a channel that is simply empty when nothing needs calling.
    genomes_ch        = Channel.empty()
    ready_proteins_ch = Channel.empty()
    def sample_count  = 0

    if (params.input) {
        if (params.subsample_genomes) {
            samples = samples.take(params.subsample_genomes as int)
            log.info "Subsampling to ${samples.size()} samples"
        }
        // A process can only be invoked once per workflow, so both kinds of
        // multi-file row go through CONCAT_INPUTS together and are told apart
        // afterwards by name. Single-file rows skip it, which is what keeps their
        // PREDICT_PROTEINS hashes identical to the form before samplesheets.
        // Written out as an if rather than a ternary: Nextflow's DSL2 compiler
        // rejects a process invocation used as an expression.
        def many = samples.findAll { it.files.size() > 1 }
        joined = Channel.empty()
        if (many)
            joined = CONCAT_INPUTS(Channel.fromList(
                many.collect { tuple(it.name, it.files, it.from_proteins ? 'faa' : 'fna') }))

        def protein_names = samples.findAll { it.from_proteins }.collect { it.name } as Set

        genomes_ch = Channel.fromList(
                samples.findAll { !it.from_proteins && it.files.size() == 1 }
                       .collect { tuple(it.name, it.files[0]) })
            .mix(joined.filter { !protein_names.contains(it[0]) })

        ready_proteins_ch = Channel.fromList(
                samples.findAll { it.from_proteins && it.files.size() == 1 }
                       .collect { tuple(it.name, it.files[0]) })
            .mix(joined.filter { protein_names.contains(it[0]) })

        sample_count = samples.size()
        log.info "Samplesheet ${params.input}: ${samples.size()} samples, " +
                 "${samples.count { !it.from_proteins }} from genomes, " +
                 "${samples.count { it.from_proteins }} from proteins"
    }
    else if (params.genomes) {
        def genome_files = expandInputs(params.genomes, '{fa,fna,fasta}')
        if (params.subsample_genomes) {
            genome_files = genome_files.take(params.subsample_genomes as int)
            log.info "Subsampling to ${genome_files.size()} genomes"
        }
        sample_count = genome_files.size()
        genomes_ch = Channel.fromList(withSampleIds(genome_files, 'genome'))
    }
    else {
        def protein_files = expandInputs(params.proteins, '{fa,faa,fasta}')
        if (params.subsample_genomes) {
            protein_files = protein_files.take(params.subsample_genomes as int)
            log.info "Subsampling to ${protein_files.size()} protein files"
        }
        sample_count = protein_files.size()
        ready_proteins_ch = Channel.fromList(withSampleIds(protein_files, 'protein FASTA'))
    }

    PREDICT_PROTEINS(genomes_ch)
    proteins_ch = PREDICT_PROTEINS.out.mix(ready_proteins_ch)

    // ---- Low-complexity masking (optional) ----
    //
    // After the mix, so one call site covers all three input routes.
    if (!params.no_masking) {
        MASK_REPEATS(proteins_ch)
        proteins_ch = MASK_REPEATS.out
    }

    // ---- BUSCO completeness (optional) ----
    //
    // Genome rows get genome mode; rows that arrived as proteins (--proteins,
    // or protein_files in a samplesheet) get protein mode instead of being
    // left out - BUSCO's own -m vocabulary names both, and the summary file
    // BUSCO writes is identical in shape either way, so RUN_BUSCO's parsing
    // needs no branching, only the mode tag threaded through the input tuple.
    if (params.busco_lineage) {
        busco_in = genomes_ch.map { id, f -> tuple(id, f, 'genome') }
            .mix(ready_proteins_ch.map { id, f -> tuple(id, f, 'proteins') })
        RUN_BUSCO(busco_in)
        BUSCO_SUMMARY(RUN_BUSCO.out.completeness.map { _id, tsv -> tsv }.collect())
        busco_completeness_ch = RUN_BUSCO.out.completeness  // tuple(id, tsv)
    }

    // ---- HMM chunks ----
    def hmm_input = file(params.hmms)
    auto_chunks = autoChunkCount(sample_count)
    if (hmm_input.isDirectory()) {
        GROUP_HMMS(params.hmms, auto_chunks)
        hmm_chunks = GROUP_HMMS.out.flatten()
    } else {
        SPLIT_HMMS(Channel.fromPath(params.hmms, checkIfExists: true), auto_chunks)
        hmm_chunks = SPLIT_HMMS.out.flatten()
    }

    // ---- Parallel hmmsearch: every protein file x every HMM chunk ----
    HMMSEARCH(proteins_ch.combine(hmm_chunks))

    // ---- Merge results per genome ----
    merge_in = HMMSEARCH.out.tblout
        .groupTuple()
        .join(HMMSEARCH.out.domtblout.groupTuple())
        .join(HMMSEARCH.out.alignments.groupTuple())
    MERGE_RESULTS(merge_in)

    // ---- Post-processing (runs once, after all merges complete) ----
    POSTPROCESS(MERGE_RESULTS.out.map { _id, tblout, _domtblout, _out -> tblout }.collect())

    // ---- BLIMMP module detection ----
    //
    // Runs whenever sigma can be established: measured per genome by BUSCO, or
    // stated once with --sigma. With neither, the run stops after annotation,
    // and says so rather than leaving you to notice the missing files.
    domtblout_ch = MERGE_RESULTS.out.map { id, _tblout, domtblout, _out -> tuple(id, domtblout) }

    // Per-sample sigma and taxonomy, looked up by name. Both empty without a
    // samplesheet, where the run-level flags supply them instead.
    def sheet_sigma    = samples ? samples.collectEntries { [(it.name): it.sigma] }    : [:]
    def sheet_taxonomy = samples ? samples.collectEntries { [(it.name): it.taxonomy] } : [:]

    // Rows BLIMMP cannot score, named now rather than left to be discovered as
    // missing files afterwards. Under --busco_lineage every row, genome or
    // protein, always gets a measurement, so only a run with neither
    // --busco_lineage nor any sigma anywhere lands rows here.
    def unscored = (samples ?: [])
        .findAll { s -> s.sigma == null && sigma_value == null && !params.busco_lineage }
        .collect { it.name }
    if (unscored)
        log.warn "BLIMMP will not run for ${unscored.join(', ')}: no completeness estimate. Fill in the sigma column on those rows, or pass --sigma. They are annotated against KOfam and stop there."

    def blimmp_in = null
    if (params.busco_lineage) {
        // BUSCO now measures every row, genome or protein (see the BUSCO
        // completeness block above). A filled-in sigma cell outranks the
        // measurement either way.
        blimmp_in = domtblout_ch
            .join(busco_completeness_ch.map { id, tsv -> tuple(id, tsv.text.tokenize()[1]) })
            .map { id, domtblout, measured ->
                tuple(id, domtblout, sheet_sigma[id] ?: measured, taxonomyFor(id, sheet_taxonomy)) }
    }
    else if (sigma_value != null || sheet_sigma.any { _name, value -> value != null }) {
        blimmp_in = domtblout_ch
            .filter { !unscored.contains(it[0]) }
            .map { id, domtblout ->
                tuple(id, domtblout, sheet_sigma[id] ?: sigma_value, taxonomyFor(id, sheet_taxonomy)) }
    }

    if (blimmp_in != null) {
        RUN_BLIMMP(blimmp_in)
        CHECK_BLIMMP(RUN_BLIMMP.out.dk)
        MODULE_MATRIX(RUN_BLIMMP.out.probabilities.collect())
    }
    else if (params.proteins) {
        // Reached without --busco_lineage and without any sigma anywhere.
        // Annotation on its own is a legitimate use of this pipeline, so this
        // is a warning rather than a refusal.
        log.warn "BLIMMP will not run: it needs --sigma or --busco_lineage to establish completeness. Use --sigma 1.0 if the proteomes are complete. Continuing with HMMER annotation only."
    }
    else {
        log.warn "BLIMMP will not run: no completeness estimate. Add --busco_lineage to measure it, or --sigma to state it. Continuing with HMMER annotation only."
    }
}

/*
 * Entry point for domtblout files that already exist:
 *
 *   nextflow run main.nf -entry blimmp_only \
 *       --domtblout 'results/*\/*.domtblout' \
 *       --completeness results/busco_summary.tsv \
 *       --taxonomy cyanobacteriota
 *
 * The search is the expensive part of this pipeline and BLIMMP takes minutes,
 * so re-scoring finished searches under different priors should not mean
 * running hmmsearch again. Sigma comes either from --completeness, which reads
 * the busco_summary.tsv a full run writes and so reproduces per-genome values
 * exactly, or from --sigma for one value across the board.
 */
workflow blimmp_only {
    preflight()

    if (!params.domtblout)
        error """-entry blimmp_only needs --domtblout.
         |  A directory of .domtblout files, or a glob such as
         |  'results/*/*.domtblout' for one genome per subdirectory.""".stripMargin()

    // Nothing else in the pipeline runs here, so flags belonging to the search
    // are refused rather than accepted and ignored.
    ['genomes', 'proteins', 'hmms', 'busco_lineage', 'input'].each { name ->
        if (params[name])
            error """--${name} does nothing with -entry blimmp_only.
             |  This entry point starts from finished HMMER tables. Drop the flag, or
             |  drop -entry to run the whole pipeline.""".stripMargin()
    }

    def sigma_value = validatedSigma()
    if (sigma_value != null && params.completeness)
        error """Provide --sigma or --completeness, not both.
         |  --completeness gives each sample its own value; --sigma applies one to
         |  all of them.""".stripMargin()
    if (sigma_value == null && !params.completeness)
        error """-entry blimmp_only needs --completeness or --sigma.
         |  --completeness results/busco_summary.tsv re-uses the values BUSCO
         |  measured. --sigma 1.0 states one value for every input.""".stripMargin()

    def pairs = withSampleIds(expandInputs(params.domtblout, 'domtblout'), 'domtblout')

    def with_sigma
    if (params.completeness) {
        def table = readCompleteness(params.completeness)
        def missing = pairs.collect { id, _f -> id }.findAll { !table.containsKey(it) }
        if (missing)
            error """No completeness value for ${missing.join(', ')} in ${params.completeness}.
             |  Sample IDs come from ${params.id_from_filename ? 'file names' : 'parent directory names'};
             |  the table has ${table.keySet().sort().join(', ')}.""".stripMargin()
        with_sigma = pairs.collect { id, f -> tuple(id, f, table[id], params.taxonomy) }
    }
    else {
        with_sigma = pairs.collect { id, f -> tuple(id, f, sigma_value, params.taxonomy) }
    }

    log.info "BLIMMP over ${with_sigma.size()} pre-computed domtblout files"

    RUN_BLIMMP(Channel.fromList(with_sigma))
    CHECK_BLIMMP(RUN_BLIMMP.out.dk)
    MODULE_MATRIX(RUN_BLIMMP.out.probabilities.collect())
}
