#!/usr/bin/env bash
#
# Check that a -profile test run actually did something.
#
#   nextflow run main.nf -profile test,docker
#   assets/test/assert.sh results_test
#
# Exiting 0 is not evidence. BLIMMP 0.1.2 shipped without its KEGG data files
# and produced output that was structurally perfect and entirely empty: correct
# columns, 340 module rows, every task green, for four months. Every check here
# asserts a number.

set -uo pipefail

outdir="${1:-results_test}"
hmm_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/kofam"
samples=(MED4 MIT9301)

# At least this fraction of the bundled KOfam profiles must hit each genome.
# The bundled set was chosen so that all of them hit both test regions, so the
# observed value is 20/20. Half is a regression guard with room to spare: it
# catches a broken search, a mis-split HMM database, or the wrong genome slice,
# without failing over one profile drifting below threshold.
min_fraction=50
evalue_cutoff=1e-6

fail=0
note () { printf '  %-34s %s\n' "$1" "$2"; }
bad  () { printf '  %-34s %s\n' "$1" "FAIL: $2"; fail=1; }

total_hmms=$(find "${hmm_dir}" -name '*.hmm' | wc -l)
if [ "${total_hmms}" -eq 0 ]; then
    echo "No HMM profiles found in ${hmm_dir}" >&2
    exit 1
fi
echo "Bundled KOfam profiles: ${total_hmms}"
echo "Output directory:       ${outdir}"
echo

for s in "${samples[@]}"; do
    echo "${s}"
    d="${outdir}/${s}"
    if [ ! -d "${d}" ]; then
        bad "output directory" "${d} does not exist"
        continue
    fi

    # 1. Did the search find the profiles we bundled?
    tblout="${d}/${s}.tblout"
    if [ ! -s "${tblout}" ]; then
        bad "hmmsearch results" "${tblout} missing or empty"
    else
        hits=$(awk -v c="${evalue_cutoff}" '!/^#/ && $5+0 < c+0 {print $3}' "${tblout}" | sort -u | wc -l)
        pct=$(( 100 * hits / total_hmms ))
        if [ "${pct}" -lt "${min_fraction}" ]; then
            bad "KOs hit at E<${evalue_cutoff}" "${hits}/${total_hmms} (${pct}%), need ${min_fraction}%"
        else
            note "KOs hit at E<${evalue_cutoff}" "${hits}/${total_hmms} (${pct}%)"
        fi
    fi

    # 2. Did BLIMMP write everything it should?
    for f in _BLIMMP_substituted_dk.csv _BLIMMP_substituted_modules.json _BLIMMP_substituted_module_steps.csv _BLIMMP_substituted_module_probabilities.csv _BLIMMP_version.txt; do
        [ -s "${d}/${s}${f}" ] || bad "output file" "${s}${f} missing or empty"
    done

    # The version file is only worth having if it names a version.
    ver="${d}/${s}_BLIMMP_version.txt"
    if [ -s "${ver}" ]; then
        v=$(head -1 "${ver}")
        [ "${v}" != "unknown" ] && note "BLIMMP version" "${v}" \
                                || bad  "BLIMMP version" "recorded as unknown"
    fi

    # 3. Did it score anything? A well-formed empty table passes every other check.
    dk="${d}/${s}_BLIMMP_substituted_dk.csv"
    if [ -s "${dk}" ]; then
        rows=$(( $(wc -l < "${dk}") - 1 ))
        [ "${rows}" -ge 1 ] && note "per-KO rows" "${rows}" \
                            || bad  "per-KO rows" "${rows}, BLIMMP scored nothing"
    fi
    echo
done

# 4. The cross-genome summary should name every sample.
summary="${outdir}/summary.tsv"
echo "summary.tsv"
if [ ! -s "${summary}" ]; then
    bad "summary.tsv" "missing or empty"
else
    for s in "${samples[@]}"; do
        grep -q "^${s}\b" "${summary}" && note "row for ${s}" "present" \
                                       || bad  "row for ${s}" "absent from summary.tsv"
    done
fi

# 5. The module matrix should carry one column per sample and every module.
matrix="${outdir}/module_matrix.tsv"
echo
echo "module_matrix.tsv"
if [ ! -s "${matrix}" ]; then
    bad "module_matrix.tsv" "missing or empty"
else
    cols=$(head -1 "${matrix}" | awk -F'\t' '{print NF}')
    want=$(( ${#samples[@]} + 2 ))
    [ "${cols}" -eq "${want}" ] && note "columns" "${cols}" \
                                || bad  "columns" "${cols}, expected ${want}"
    rows=$(( $(wc -l < "${matrix}") - 1 ))
    [ "${rows}" -ge 300 ] && note "module rows" "${rows}" \
                          || bad  "module rows" "${rows}, expected the full KEGG module set"
    # A matrix of nothing but NA would pass both checks above.
    scored=$(awk -F'\t' 'NR>1 { for (i=3; i<=NF; i++) if ($i != "NA" && $i+0 > 0) { c++; break } } END { print c+0 }' "${matrix}")
    [ "${scored}" -ge 1 ] && note "modules with a score" "${scored}" \
                          || bad  "modules with a score" "0, every cell is NA or zero"
fi

echo
if [ "${fail}" -eq 0 ]; then
    echo "PASS"
else
    echo "FAIL"
fi
exit "${fail}"
