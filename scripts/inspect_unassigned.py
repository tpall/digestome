#!/usr/bin/env python3
"""Inspect genomes that carry mcrA but received no methanogenesis route.

These are the informative failures. A genome with no mcrA and no route is simply
incomplete; a genome with mcrA and no route means the panel can see a methanogen
and cannot say what kind, which is either a missing detector or real biology the
gates do not cover. The obligate methylotrophs among Methanofastidiosia were
found exactly this way.

Usage:
  inspect_unassigned.py --work <gtdb work dir> [--stratum methanogen]
"""
import argparse
import collections
import csv
import pathlib
import re
import sys

GATE = re.compile(r'\[([x ?])\]\s+(methanogenesis:\s*\S+)')
MARKER = re.compile(r'^\s{2,}(\S+)\s+([+-])\s*$', re.M)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--work', required=True)
    ap.add_argument('--stratum', default='methanogen')
    ap.add_argument('--require', default='mcrA',
                    help='only report genomes where this marker is present (default mcrA)')
    args = ap.parse_args()

    W = pathlib.Path(args.work)
    meta = {r['accession']: r for r in csv.DictReader(open(W / 'gtdb_sample.tsv'), delimiter='\t')}

    hits = []
    for f in sorted((W / 'scored').glob('*.summary.txt')):
        acc = f.name[:-len('.summary.txt')]
        if meta.get(acc, {}).get('stratum') != args.stratum:
            continue
        txt = f.read_text()
        if any(mark == 'x' for mark, _ in GATE.findall(txt)):
            continue
        marks = dict(MARKER.findall(txt))
        if marks.get(args.require) != '+':
            continue
        hits.append((acc, marks, txt))

    if not hits:
        sys.exit(f'no {args.stratum} genomes with {args.require} and no route')

    print(f'{len(hits)} {args.stratum} genomes carry {args.require} but received no route\n')
    missing = collections.Counter()
    for acc, marks, txt in hits:
        m = meta[acc]
        ranks = m['lineage'].split(';')
        genus = ranks[5] if len(ranks) > 5 else '?'
        present = sorted(g for g, s in marks.items() if s == '+')
        absent = sorted(g for g, s in marks.items() if s == '-')
        for g in absent:
            missing[g] += 1
        print(f'  {acc}  {m["order"]} / {genus}')
        print(f'      present: {" ".join(present) or "(none)"}')
        print(f'      absent : {" ".join(absent) or "(none)"}')
        na = [line.strip() for line in txt.splitlines() if 'NOT ASSESSABLE' in line]
        for line in na[:2]:
            print(f'      {line}')
        print()

    print('markers absent across this group (a shared gap points at a missing detector):')
    for g, c in missing.most_common():
        print(f'  {g}: {c} of {len(hits)}')

    print('\norders represented:')
    for o, c in collections.Counter(meta[a]['order'] for a, _, _ in hits).most_common():
        print(f'  {o}: {c}')


if __name__ == '__main__':
    main()
