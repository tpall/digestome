#!/usr/bin/env python3
"""Summarise the GTDB representative benchmark from scored per-genome output.

Reads the `[x] / [ ] / [?]` branch-gate lines that panel_scored.py prints. A gate
is only counted as assigned when it is marked `[x]`: `[?]` means not assessable,
which is a distinct state from absent and must not be collapsed into either.

Usage:
  summarise_gtdb_benchmark.py --work <gtdb work dir> [--scored DIR] [--pf05369 pf05369_hits.txt]

Also prints route calls by GTDB lineage (acetate by genus, methyl and H2/CO2 by order), which is
the check that a methanogen got the right route, not only that a non-methanogen got none.
"""
import argparse
import collections
import csv
import pathlib
import re
import sys

GATE = re.compile(r'\[([x ?])\]\s+(methanogenesis:\s*\S+)')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--work', required=True)
    ap.add_argument('--scored', default='scored', help='scored subdirectory of --work (default: scored)')
    ap.add_argument('--pf05369', default=None,
                    help='optional accession<TAB>count file for the MtmB family')
    args = ap.parse_args()

    W = pathlib.Path(args.work)
    meta = {r['accession']: r for r in csv.DictReader(open(W / 'gtdb_sample.tsv'), delimiter='\t')}

    routes, notassessable = collections.defaultdict(set), collections.defaultdict(set)
    SC = W / args.scored
    for f in SC.glob('*.summary.txt'):
        acc = f.name[:-len('.summary.txt')]
        for mark, label in GATE.findall(f.read_text()):
            branch = label.split(':')[1].strip()
            if mark == 'x':
                routes[acc].add(branch)
            elif mark == '?':
                notassessable[acc].add(branch)

    scored = {f.name[:-len('.modules.tsv')] for f in SC.glob('*.modules.tsv')}
    if not scored:
        sys.exit(f'no scored output under {SC}')

    tally = collections.Counter()
    for a in scored:
        tally[(meta.get(a, {}).get('stratum', '?'), bool(routes.get(a)))] += 1

    print('============ GTDB representative benchmark ============')
    print(f'genomes scored: {len(scored)}')
    for st in ('background', 'methanogen'):
        print(f'  {st:11s} route assigned: {tally[(st, True)]:5d}   no route: {tally[(st, False)]:5d}')

    fp = sorted(a for a in scored if meta.get(a, {}).get('stratum') == 'background' and routes.get(a))
    print(f'\nFALSE POSITIVES (background genomes given a route): {len(fp)}')
    for a in fp[:25]:
        m = meta[a]
        print(f'  {a}  {m["phylum"]} / {m["order"]} -> {sorted(routes[a])}')

    fn = sorted(a for a in scored if meta.get(a, {}).get('stratum') == 'methanogen' and not routes.get(a))
    print(f'\nMETHANOGENS WITH NO ROUTE: {len(fn)}')
    for o, c in collections.Counter(meta[a]['order'] for a in fn).most_common(12):
        print(f'  {o}: {c}')

    print('\nroute assignments among methanogens:')
    mix = collections.Counter(r for a in scored
                              if meta.get(a, {}).get('stratum') == 'methanogen'
                              for r in routes.get(a, ()))
    for r, c in mix.most_common():
        print(f'  {r}: {c}')

    lineage = {}
    tax = W / 'gtdbtk_sample.tsv'
    if tax.exists():
        for line in tax.read_text().splitlines():
            f = line.split('\t')
            if len(f) > 1 and f[0].startswith('GC'):
                lineage[f[0]] = f[1]
    def rank(a, r):
        for tok in lineage.get(a, '').split(';'):
            if tok.startswith(r + '__'):
                return re.sub(r'_[A-Z]+$', '', tok[3:]) or '?'
        return '?'
    if lineage:
        print('\nroute calls by GTDB lineage (all genomes):')
        for route, r in (('acetoclastic', 'g'), ('methylotrophic', 'o'), ('hydrogenotrophic', 'o')):
            c = collections.Counter(rank(a, r) for a in scored if route in routes.get(a, ()))
            level = {'g': 'genus', 'o': 'order'}[r]
            print(f'  {route} by {level} ({sum(c.values())}): ' +
                  ', '.join(f'{k} {v}' for k, v in c.most_common()))

    if args.pf05369:
        p = pathlib.Path(args.pf05369)
        if p.exists():
            hits = {}
            for line in p.read_text().splitlines():
                if line.strip():
                    a, n = line.split('\t')
                    hits[a] = int(n)
            byst = collections.Counter(meta.get(a, {}).get('stratum', '?') for a in hits)
            print(f'\nPF05369 (MtmB), not a panel marker: {len(hits)} genomes with a hit')
            print(f'  in methanogens: {byst["methanogen"]} of '
                  f'{sum(1 for a in scored if meta.get(a,{}).get("stratum")=="methanogen")}')
            print(f'  in background:  {byst["background"]} of '
                  f'{sum(1 for a in scored if meta.get(a,{}).get("stratum")=="background")}'
                  f'   <- these are the false positives')
            bg = sorted(a for a in hits if meta.get(a, {}).get('stratum') == 'background')
            for a in bg[:25]:
                m = meta[a]
                print(f'    {a}  {m["phylum"]} / {m["order"]}  ({hits[a]} protein(s))')
            print('  hits by methanogen order:')
            for o, c in collections.Counter(
                    meta[a]['order'] for a in hits
                    if meta.get(a, {}).get('stratum') == 'methanogen').most_common(12):
                print(f'    {o}: {c}')


if __name__ == '__main__':
    main()
