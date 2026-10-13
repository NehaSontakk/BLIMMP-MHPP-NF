#!/usr/bin/env python3
"""
postprocess.py - Summarize HMMER tblout results across genomes.

For each genome, counts hits from the tblout (per-sequence) output:
  - matches with full-sequence E-value < 1
  - matches with full-sequence E-value < 1e-6

Writes a tab-separated table to stdout:
  genome_name  count<1  count<1e-6

Usage:
    postprocess.py results/*.tblout > summary.tsv
"""

import sys
import os


def count_hits(tblout_path):
    """Return (count_lt1, count_lt1e6) for hits in a tblout file."""
    count_lt1 = 0
    count_lt1e6 = 0
    with open(tblout_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 5:
                continue
            try:
                evalue = float(fields[4])   # column 5: full-sequence E-value
            except ValueError:
                continue
            if evalue < 1:
                count_lt1 += 1
            if evalue < 1e-6:
                count_lt1e6 += 1
    return count_lt1, count_lt1e6


def main(paths):
    print("genome_name\tcount<1\tcount<1e-6")
    for path in sorted(paths):
        genome_name = os.path.basename(path)
        if genome_name.endswith('.tblout'):
            genome_name = genome_name[:-7]
        count_lt1, count_lt1e6 = count_hits(path)
        print(f"{genome_name}\t{count_lt1}\t{count_lt1e6}")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("Usage: postprocess.py <file.tblout> [file2.tblout ...]")
    main(sys.argv[1:])
