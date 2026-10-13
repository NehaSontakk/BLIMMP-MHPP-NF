#!/usr/bin/env python3
"""
group_hmms.py - Concatenate individual HMM files into chunks for parallel search.

Reads all *.hmm files from a directory, sorts them, and writes groups of files
into chunk_0000.hmm, chunk_0001.hmm, ...  Each chunk is a concatenated HMM
file suitable for hmmsearch.

Chunk size is determined by (in order of precedence):
  --chunks N      divide files evenly into N chunks (ceiling division)
  --chunk-size N  fixed number of files per chunk
  (neither)       default of 200 files per chunk

Usage:
    group_hmms.py /path/to/kofam/
    group_hmms.py --chunks 50 /path/to/kofam/
    group_hmms.py --chunk-size 100 /path/to/kofam/
    group_hmms.py --chunk-size 100 --limit 500 /path/to/kofam/
"""

import argparse
import glob
import math
import os
import sys
from concurrent.futures import ThreadPoolExecutor

MAX_CHUNKS = 9999   # chunk_{:04d} supports 0000–9999
MIN_CHUNK_SIZE = 500  # floor when --chunks is given, see below

# Concatenating a KOfam-sized profile directory is bound by per-file latency,
# not by any work: on a networked filesystem the open costs far more than the
# read, and the processor is idle throughout. Reading through a thread pool
# hides that latency, and threads are enough because the GIL is released
# during I/O. Measured on a GPFS scratch mount, 1500 profiles: 84 seconds
# serially against 12 with 16 workers, for byte-identical output. On a local
# disk there is little latency to hide and this neither helps nor hurts much.
#
# Raise BLIMMP_READ_WORKERS on a filesystem with worse latency, or set it to 1
# to read serially.
READ_WORKERS = max(1, int(os.environ.get('BLIMMP_READ_WORKERS', '16')))

# Files held in memory at once. Mapping the pool over a whole chunk would read
# it entirely into memory, which for KOfam is over a gigabyte per chunk.
READ_BATCH = 64


def read_hmm(path):
    """Read one profile whole. Separate so a thread pool can map over it."""
    with open(path) as fh:
        return fh.read()


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        '--chunks', type=int, metavar='N',
        help='total number of output chunks (chunk size is auto-calculated)',
    )
    mode.add_argument(
        '--chunk-size', type=int, default=200, metavar='N',
        help='number of HMM files per chunk (default: 200)',
    )
    parser.add_argument(
        '--limit', type=int, default=0, metavar='N',
        help='use only the first N HMM files; 0 means no limit (default: 0)',
    )
    parser.add_argument('hmm_dir', help='directory of individual .hmm files')
    args = parser.parse_args()

    # Collect and sort all .hmm files in the directory
    all_files = sorted(glob.glob(os.path.join(args.hmm_dir, '*.hmm')))
    if not all_files:
        sys.exit(f"Error: no .hmm files found in {args.hmm_dir}")

    # Apply subsample limit before calculating chunk size
    files = all_files[:args.limit] if args.limit else all_files
    total = len(files)

    if args.chunks:
        if args.chunks > MAX_CHUNKS:
            sys.exit(
                f"Error: --chunks {args.chunks} exceeds the maximum of {MAX_CHUNKS}.\n"
                f"Use --chunks {MAX_CHUNKS} or fewer."
            )
        chunk_size = math.ceil(total / args.chunks)
        # A chunk small enough to finish in a couple of seconds costs more in
        # scheduler round trips than it saves in parallelism. Nine genomes
        # against 27262 profiles at 100 per chunk made 2457 searches of two to
        # three seconds each, and spent the rest of the day queueing. The floor
        # only ever reduces the chunk count, so an explicit --chunk-size is
        # untouched by it.
        if chunk_size < MIN_CHUNK_SIZE:
            chunk_size = min(MIN_CHUNK_SIZE, total) or 1
    else:
        chunk_size = args.chunk_size

    expected_chunks = math.ceil(total / chunk_size)
    if expected_chunks > MAX_CHUNKS:
        min_chunk_size = math.ceil(total / MAX_CHUNKS)
        sys.exit(
            f"Error: {total} HMM files with --chunk-size {chunk_size} would produce "
            f"{expected_chunks} chunks, exceeding the maximum of {MAX_CHUNKS}.\n"
            f"Use --chunk-size {min_chunk_size} or larger."
        )

    # Write each chunk as a concatenated HMM file. pool.map preserves the order
    # of its input, so a chunk holds exactly what the serial version would.
    with ThreadPoolExecutor(max_workers=READ_WORKERS) as pool:
        for chunk_idx, start in enumerate(range(0, total, chunk_size)):
            chunk_files = files[start:start + chunk_size]
            out_path = f'chunk_{chunk_idx:04d}.hmm'
            with open(out_path, 'w') as out:
                for i in range(0, len(chunk_files), READ_BATCH):
                    for data in pool.map(read_hmm, chunk_files[i:i + READ_BATCH]):
                        out.write(data)


if __name__ == '__main__':
    main()
