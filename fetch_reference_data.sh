#!/usr/bin/env bash
#
# Fetch the reference data a real run needs, from the places it comes from.
#
#   ./fetch_reference_data.sh --kofam /path/to/kofam
#   ./fetch_reference_data.sh --genomes /path/to/genomes
#   ./fetch_reference_data.sh --busco prochlorococcus_odb12 /path/to/busco_downloads
#
# Any combination, in one call or several. Each part checks for what it would
# download and skips it, so re-running after an interruption resumes rather than
# starts over.
#
# KOfam comes from a pinned monthly snapshot. KOFAM_RELEASE=current takes the
# newest instead; see the note above that variable for what that costs.
#
# Nothing here is transferred from a person. All of it is public and versioned
# upstream, which is what makes a port reproducible: the next cluster runs the
# same command rather than waiting for someone to send a copy.

set -uo pipefail

# KOfam is published as dated monthly snapshots under archives/, alongside a
# profiles.tar.gz at the top of the tree that is whatever is newest. Pin a
# snapshot, because the rolling copy changes under you: a download on
# 2026-08-28 held 27,864 profiles where one a few weeks earlier held 27,262,
# and a KO count that differs between two clusters should say something about
# the clusters rather than about the day each of them downloaded.
#
# Set KOFAM_RELEASE=current to take the newest instead, and give up that
# guarantee. Available dates are listed at
# https://www.genome.jp/ftp/db/kofam/archives/
KOFAM_RELEASE="${KOFAM_RELEASE:-2026-07-02}"

if [ "${KOFAM_RELEASE}" = current ]; then
    KOFAM_BASE=https://www.genome.jp/ftp/db/kofam
else
    KOFAM_BASE="https://www.genome.jp/ftp/db/kofam/archives/${KOFAM_RELEASE}"
fi
KOFAM_URL="${KOFAM_BASE}/profiles.tar.gz"
KOFAM_LIST_URL="${KOFAM_BASE}/ko_list.gz"
EFETCH=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi

# The nine Prochlorococcus genomes the pipeline was developed against, as
# name:accession. All are closed single-contig RefSeq assemblies, which is why
# docs/overview.md tells you to run them with --sigma 1.0 rather than measuring
# completeness. Replace this list with your own; the format is the layout the
# pipeline expects, one genome per subdirectory.
GENOMES=(
    MED4:NC_005072.1
    MIT0801:NZ_CP007754.1
    MIT9211:NC_009976.1
    MIT9301:NC_009091.1
    MIT9303:NC_008820.1
    MIT9312:NC_007577.1
    MIT9313:NC_005071.1
    NATL2A:NC_007335.2
    SS120:NC_005042.1
)

die () { echo "error: $*" >&2; exit 1; }
say () { printf '%s\n' "$*"; }

usage () {
    sed -n '3,18p' "$0" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}

fetch_kofam () {
    local dir="$1"
    mkdir -p "${dir}" || die "cannot create ${dir}"

    local have
    have=$(find "${dir}" -maxdepth 1 -name '*.hmm' 2>/dev/null | head -30000 | wc -l)
    if [ "${have}" -gt 20000 ]; then
        say "KOfam: ${have} profiles already in ${dir}, skipping"
        return 0
    fi

    say "KOfam: downloading profiles.tar.gz, 1.6 GB compressed and about 7.2 GB extracted"
    say "  from ${KOFAM_URL}"

    # Resume rather than restart. This is 1.6 GB from a server that hands it
    # over at about 400 kB/s, so the transfer runs for the better part of an
    # hour and is long enough to be reset partway through: on Stampede3 it died
    # at 1.32 GB with "OpenSSL SSL_read: Connection reset by peer". curl's own
    # --retry does not help, because without -C - each retry reopens the same
    # -o file at byte zero.
    #
    # The skip check at the top of this function counts .hmm files, so it does
    # not see a part-downloaded tarball either, and the header of this script
    # says an interrupted download resumes, which was not true of this part.
    # This fixes both: the size check below skips a tarball that is already
    # whole, and -C - continues one that is not.
    local want have_bytes
    want=$(curl -fsIL "${KOFAM_URL}" 2>/dev/null \
           | tr -d '\r' | awk 'BEGIN{IGNORECASE=1} /^content-length:/ {n=$2} END {print n}')
    have_bytes=$(wc -c < "${dir}/profiles.tar.gz" 2>/dev/null || echo 0)

    if [ -n "${want}" ] && [ "${have_bytes}" = "${want}" ]; then
        say "KOfam: profiles.tar.gz already complete, ${want} bytes, skipping download"
    else
        [ "${have_bytes}" -gt 0 ] && say "KOfam: resuming from ${have_bytes} of ${want:-?} bytes"
        curl -fL --retry 3 --retry-delay 5 -C - -o "${dir}/profiles.tar.gz" "${KOFAM_URL}" \
            || die "download failed part-way. ${dir}/profiles.tar.gz holds what arrived; re-run this script to resume from there rather than deleting it. If outbound HTTP is blocked here, fetch it elsewhere and copy it in"
    fi

    say "KOfam: extracting"
    tar xzf "${dir}/profiles.tar.gz" -C "${dir}" || die "extract failed"

    # The tarball unpacks into a profiles/ subdirectory. The pipeline wants a
    # directory of .hmm files, so flatten it.
    #
    # Moved with find rather than `mv profiles/*.hmm`: 27,864 paths is 2.2 MB
    # of argument list against a 2 MB ARG_MAX, so the glob form dies with
    # "Argument list too long". With its stderr discarded that failure was
    # invisible, and the count check below then reported a tarball holding no
    # profiles at all, which points at the download rather than at the move.
    if [ -d "${dir}/profiles" ]; then
        find "${dir}/profiles" -maxdepth 1 -name '*.hmm' -exec mv -t "${dir}/" {} + \
            || die "could not move profiles out of ${dir}/profiles into ${dir}"
        rmdir "${dir}/profiles" 2>/dev/null || say "  note: ${dir}/profiles still holds non-.hmm files"
    fi
    rm -f "${dir}/profiles.tar.gz"

    have=$(find "${dir}" -maxdepth 1 -name '*.hmm' | wc -l)
    [ "${have}" -gt 20000 ] || die "only ${have} profiles after extraction, expected over 20000"

    # Leave the release on disk. A directory of 27,864 anonymous .hmm files
    # otherwise cannot be told apart from a directory of 27,262 of them, and
    # that difference is the first thing to check when two sites disagree.
    printf '%s\n' "${KOFAM_RELEASE}" > "${dir}/KOFAM_RELEASE"
    say "KOfam: ${have} profiles in ${dir}, release ${KOFAM_RELEASE}"

    # ko_list is not read by this pipeline. It is small, and it is what you need
    # to look up what a KO actually is, so it comes along.
    curl -fsL -o "${dir}/ko_list.gz" "${KOFAM_LIST_URL}" && gunzip -f "${dir}/ko_list.gz" \
        && say "KOfam: ko_list alongside it, for looking up KO names"
}

fetch_genomes () {
    local dir="$1"
    mkdir -p "${dir}" || die "cannot create ${dir}"
    local n=0 skipped=0

    for entry in "${GENOMES[@]}"; do
        local name="${entry%%:*}" acc="${entry##*:}"
        local out="${dir}/${name}/${name}.fna"
        if [ -s "${out}" ]; then
            skipped=$((skipped + 1)); continue
        fi
        mkdir -p "${dir}/${name}"
        curl -fsL --retry 3 --retry-delay 3 \
            "${EFETCH}?db=nuccore&id=${acc}&rettype=fasta&retmode=text" -o "${out}" \
            || { rm -f "${out}"; die "failed to fetch ${name} (${acc})"; }
        # An efetch error comes back as a 200 with an HTML or empty body.
        head -1 "${out}" | grep -q '^>' || { rm -f "${out}"; die "${name} (${acc}) did not return FASTA"; }
        n=$((n + 1))
        say "  ${name}  ${acc}  $(grep -v '^>' "${out}" | tr -d '\n' | wc -c | awk '{printf "%.2f Mb", $1/1e6}')"
        sleep 1   # NCBI asks for no more than three requests a second without a key
    done
    say "genomes: ${n} downloaded, ${skipped} already present, in ${dir}"
    say "  pass them as --genomes '${dir}/*/*.fna', quoted so the pipeline expands the glob"
}

fetch_busco () {
    local lineage="$1" dir="$2"
    command -v busco >/dev/null || die "busco is not on PATH. Run this inside the BUSCO container, or use the pipeline's own --busco_downloads workflow"
    if [ -d "${dir}/lineages/${lineage}" ]; then
        say "BUSCO: ${lineage} already in ${dir}, skipping"
        return 0
    fi
    mkdir -p "${dir}"
    say "BUSCO: downloading ${lineage}"
    busco --download "${lineage}" --download_path "${dir}" || die "busco --download failed"
    say "  pass it as --busco_lineage ${lineage} --busco_downloads ${dir}"
    say "  and see docs/overview.md first: on a closed genome, --sigma 1.0 is the better answer"
}

[ $# -eq 0 ] && usage 0
command -v curl >/dev/null || die "curl is required"

while [ $# -gt 0 ]; do
    case "$1" in
        --kofam)   [ $# -ge 2 ] || die "--kofam needs a directory";   fetch_kofam "$2";   shift 2 ;;
        --genomes) [ $# -ge 2 ] || die "--genomes needs a directory"; fetch_genomes "$2"; shift 2 ;;
        --busco)   [ $# -ge 3 ] || die "--busco needs a lineage and a directory"
                   fetch_busco "$2" "$3"; shift 3 ;;
        -h|--help) usage 0 ;;
        *)         die "unknown argument: $1" ;;
    esac
done
