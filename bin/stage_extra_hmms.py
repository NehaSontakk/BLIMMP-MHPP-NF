#!/usr/bin/env python3
"""
stage_extra_hmms.py - Rename custom HMM profile(s) to fit alongside a KOfam
directory, without copying or modifying that directory itself.

hmmsearch's domtblout reports a profile's own internal NAME line as the query
name, not its filename, so a custom profile has to be renamed on that line --
not just given a KNNNNN.hmm filename -- to be treated identically to a real
KOfam entry by the rest of the pipeline. Renamed ids come from a reserved
block, K99002 and up (K99001 is BLIMMP_MHPP's own addition -- see
assets/K99001.hmm in that project), which real KEGG KO numbering has not
reached, so a custom profile can never collide with a real one.

Writes two files: a single merged .hmm file ready to hand to hmmsearch (so the
rest of the pipeline sees it as one more chunk), and a mapping.tsv of
synthetic_id -> original_name, which main.nf's workflow.onComplete reads
afterward to put the original name back in that run's own result files.

Usage:
    stage_extra_hmms.py --base-hmms /path/to/kofam extra.hmm
    stage_extra_hmms.py --base-hmms /path/to/kofam /path/to/extra_profiles/
"""

import argparse
import glob
import os
import re
import sys

RESERVED_START = 99002  # see module docstring

NAME_RE = re.compile(r'^NAME\s+(\S+)', re.MULTILINE)


def iter_hmm_records(text):
    """Yield each profile's raw text from a file that may hold one or many,
    concatenated the way a KOfam profiles.hmm or a multi-profile custom file
    would be. Each record runs from its own header line to its own '//'."""
    record_lines = []
    for line in text.splitlines(keepends=True):
        record_lines.append(line)
        if line.strip() == '//':
            yield ''.join(record_lines)
            record_lines = []
    if any(l.strip() for l in record_lines):
        raise ValueError("file does not end on a '//' terminator -- truncated or not an HMM file?")


def existing_synthetic_ids(base_hmms):
    """Scan the base --hmms path (a directory of .hmm files, or one big
    concatenated file) for NAME lines already in the reserved block, without
    loading a multi-gigabyte KOfam file into memory at once."""
    ids = set()
    if os.path.isdir(base_hmms):
        for path in glob.glob(os.path.join(base_hmms, 'K99[0-9][0-9][0-9].hmm')):
            ids.add(os.path.splitext(os.path.basename(path))[0])
        return ids
    with open(base_hmms) as fh:
        for line in fh:
            if line.startswith('NAME'):
                m = re.match(r'^NAME\s+(K99\d{3})\b', line)
                if m:
                    ids.add(m.group(1))
    return ids


def next_free_id(used):
    n = RESERVED_START
    while f'K{n}' in used:
        n += 1
    return f'K{n}'


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        '--base-hmms', required=True,
        help="the --hmms path already in use, scanned only to avoid id collisions, never copied")
    ap.add_argument('--out-hmm', default='extra_hmms_chunk.hmm')
    ap.add_argument('--out-mapping', default='extra_hmms_mapping.tsv')
    ap.add_argument('extra_hmms', help='a single .hmm file or a directory of them')
    args = ap.parse_args()

    if os.path.isdir(args.extra_hmms):
        sources = sorted(glob.glob(os.path.join(args.extra_hmms, '*.hmm')))
        if not sources:
            sys.exit(f"error: no .hmm files found in {args.extra_hmms}")
    elif os.path.isfile(args.extra_hmms):
        sources = [args.extra_hmms]
    else:
        sys.exit(f"error: {args.extra_hmms} is neither a file nor a directory")

    used = existing_synthetic_ids(args.base_hmms)
    mapping = []

    with open(args.out_hmm, 'w') as out:
        for src in sources:
            with open(src) as fh:
                text = fh.read()
            if not text.lstrip().startswith('HMMER3/'):
                sys.exit(f"error: {src} does not look like an HMMER3 profile (no HMMER3/ header)")
            for record in iter_hmm_records(text):
                m = NAME_RE.search(record)
                original_name = m.group(1) if m else os.path.splitext(os.path.basename(src))[0]
                synthetic_id = next_free_id(used)
                used.add(synthetic_id)
                renamed = NAME_RE.sub(f'NAME  {synthetic_id}', record, count=1)
                out.write(renamed)
                mapping.append((synthetic_id, original_name, os.path.abspath(src)))
                print(f"extra-hmm: '{original_name}' -> {synthetic_id} (from {src})", file=sys.stderr)

    with open(args.out_mapping, 'w') as f:
        f.write('synthetic_id\toriginal_name\toriginal_source\n')
        for row in mapping:
            f.write('\t'.join(row) + '\n')


if __name__ == '__main__':
    main()
