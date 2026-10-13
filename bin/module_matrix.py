#!/usr/bin/env python3
"""
module_matrix.py - Build a module x genome table of BLIMMP module confidences.

Reads the per-genome <name>_BLIMMP_module_probabilities.csv files that BLIMMP
writes and pivots them into one table: one row per KEGG module, one column per
genome, holding module_confidence. Modules absent from a genome's file are NA.

Writes tab-separated output to stdout:
  module  description  <genome1>  <genome2>  ...

Usage:
    module_matrix.py A_BLIMMP_module_probabilities.csv B_...csv > module_matrix.tsv
"""

import csv
import os
import sys

SUFFIX = '_BLIMMP_module_probabilities.csv'


def sample_name(path):
    """Recover the sample ID BLIMMP was given from the file it wrote."""
    base = os.path.basename(path)
    if base.endswith(SUFFIX):
        return base[:-len(SUFFIX)]
    return os.path.splitext(base)[0]


def read_one(path):
    """Return {module: (description, confidence)} for one genome."""
    out = {}
    with open(path, newline='') as fh:
        for row in csv.DictReader(fh):
            module = row.get('module')
            if not module:
                continue
            out[module] = (row.get('module_description', ''),
                           row.get('module_confidence', ''))
    return out


def main(paths):
    per_sample = {}
    for path in paths:
        name = sample_name(path)
        # Two files claiming the same sample would silently overwrite each
        # other, and the result would look complete.
        if name in per_sample:
            sys.exit(f"module_matrix.py: two files both name the sample {name}")
        per_sample[name] = read_one(path)

    samples = sorted(per_sample)
    descriptions = {}
    for name in samples:
        for module, (description, _conf) in per_sample[name].items():
            descriptions.setdefault(module, description)

    writer = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
    writer.writerow(['module', 'description'] + samples)
    for module in sorted(descriptions):
        row = [module, descriptions[module]]
        for name in samples:
            hit = per_sample[name].get(module)
            row.append(hit[1] if hit else 'NA')
        writer.writerow(row)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("Usage: module_matrix.py <file_BLIMMP_module_probabilities.csv> [...]")
    main(sys.argv[1:])
