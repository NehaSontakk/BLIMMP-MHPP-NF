#!/usr/bin/env python3
"""
split_hmms.py - Split a concatenated HMM file into chunks for parallel search.

In a concatenated HMM file, each model ends with a '//' line.  This script
splits the file into chunk_0000.hmm, chunk_0001.hmm, ... such that each chunk
contains at most CHUNK_SIZE models.

Chunk size is determined by (in order of precedence):
  --chunks N      divide models evenly into N chunks (ceiling division)
  --chunk-size N  fixed number of models per chunk
  (neither)       default of 200 models per chunk

Usage:
    split_hmms.py input.hmm
    split_hmms.py --chunks 50 input.hmm
    split_hmms.py --chunk-size 100 input.hmm
    split_hmms.py --chunk-size 100 --limit 500 input.hmm
"""

import argparse
import math
import os
import sys

MAX_CHUNKS = 9999   # chunk_{:04d} supports 0000–9999
MIN_CHUNK_SIZE = 500  # floor when --chunks is given, see below


def count_models(path):
    """Count HMM models by counting '//' delimiter lines."""
    count = 0
    with open(path) as fh:
        for line in fh:
            if line.rstrip() == '//':
                count += 1
    return count


def split(input_path, chunk_size, limit):
    """Write chunk_NNNN.hmm files, each containing at most chunk_size models."""
    chunk_idx = 0
    models_written = 0
    out = open(f'chunk_{chunk_idx:04d}.hmm', 'w')

    with open(input_path) as fh:
        for line in fh:
            # Stop early if a subsample limit was requested
            if limit and models_written >= limit:
                break

            out.write(line)

            if line.rstrip() == '//':
                models_written += 1
                # Roll over to a new chunk file when this one is full
                if models_written % chunk_size == 0:
                    out.close()
                    chunk_idx += 1
                    out = open(f'chunk_{chunk_idx:04d}.hmm', 'w')

    out.close()

    # The last file is empty when total models is an exact multiple of chunk_size
    last = f'chunk_{chunk_idx:04d}.hmm'
    if os.path.getsize(last) == 0:
        os.remove(last)


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
        help='number of models per chunk (default: 200)',
    )
    parser.add_argument(
        '--limit', type=int, default=0, metavar='N',
        help='stop after the first N models; 0 means no limit (default: 0)',
    )
    parser.add_argument('input', help='concatenated HMM file to split')
    args = parser.parse_args()

    # Always count models so we can validate the chunk count before writing.
    total = count_models(args.input)
    if args.limit:
        total = min(total, args.limit)

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
            f"Error: {total} models with --chunk-size {chunk_size} would produce "
            f"{expected_chunks} chunks, exceeding the maximum of {MAX_CHUNKS}.\n"
            f"Use --chunk-size {min_chunk_size} or larger."
        )

    split(args.input, chunk_size, args.limit)


if __name__ == '__main__':
    main()
