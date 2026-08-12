#!/usr/bin/env python3
"""Quantify what the contig-level rule costs, by comparing two scored directories.

--contig-level requires every term of a gate on one contig. That buys
single-organism evidence and pays for it in operons split across contigs. On
binned genomes the payment is pure loss, because a bin is already one organism;
the rule exists for whole-assembly input, where genome-level pooling has no
defensible error model. This measures the size of the payment.

The loss should scale with assembly fragmentation, so the report breaks it down
by contig count as well as in total.

Usage:
  compare_contig_mode.py --work <gtdb work dir> --genome-dir scored --contig-dir scored_contig
"""
import argparse
import collections
import csv
import pathlib
import re
import sys

GATE = re.compile(r'\[([x ?])\]\s+((?:methanogenesis|acetogenesis|syntrophy|substrate):\s*[^\n—]+)')


def read(d):
    """-> {accession: {gate_label: mark}} where mark is x, space or ?"""
    out = {}
    for f in d.glob('*.summary.txt'):
        acc = f.name[:-len('.summary.txt')]
        out[acc] = {lab.strip(): mark for mark, lab in GATE.findall(f.read_text())}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--work', required=True)
    ap.add_argument('--genome-dir', default='scored')
    ap.add_argument('--contig-dir', default='scored_contig')
    args = ap.parse_args()

    W = pathlib.Path(args.work)
    gen, ctg = read(W / args.genome_dir), read(W / args.contig_dir)
    common = sorted(set(gen) & set(ctg))
    if not common:
        sys.exit('no genomes scored in both directories')
    meta = {}
    sample = W / 'gtdb_sample.tsv'
    if sample.exists():
        meta = {r['accession']: r for r in csv.DictReader(open(sample), delimiter='\t')}

    # contig count per genome, from the proteome, so the loss can be read against
    # how fragmented the assembly is
    contigs = {}
    pdir = W / 'proteomes_prodigal'
    if pdir.exists():
        for acc in common:
            f = pdir / f'{acc}.faa'
            if not f.exists():
                continue
            seen = set()
            for line in f.read_text().splitlines():
                if line.startswith('>'):
                    pid = line[1:].split()[0]
                    m = re.match(r'^(.*)_\d+$', pid)
                    seen.add(m.group(1) if m else pid)
            contigs[acc] = len(seen)

    print(f'genomes compared: {len(common)}\n')

    labels = sorted({l for d in (gen, ctg) for v in d.values() for l in v})
    print(f'{"gate":48s} {"genome":>7s} {"contig":>7s} {"lost":>6s} {"gained":>7s}')
    tot_lost = 0
    for lab in labels:
        g = sum(1 for a in common if gen[a].get(lab) == 'x')
        c = sum(1 for a in common if ctg[a].get(lab) == 'x')
        lost = sum(1 for a in common if gen[a].get(lab) == 'x' and ctg[a].get(lab) != 'x')
        gained = sum(1 for a in common if ctg[a].get(lab) == 'x' and gen[a].get(lab) != 'x')
        tot_lost += lost
        print(f'{lab[:48]:48s} {g:7d} {c:7d} {lost:6d} {gained:7d}')

    def routed(d, a):
        return any(v == 'x' for k, v in d[a].items() if k.startswith('methanogenesis'))

    print('\nany methanogenesis route, by stratum:')
    for st in ('methanogen', 'background'):
        sub = [a for a in common if meta.get(a, {}).get('stratum') == st] if meta else common
        if not sub:
            continue
        g = sum(1 for a in sub if routed(gen, a))
        c = sum(1 for a in sub if routed(ctg, a))
        print(f'  {st:11s} genome {g:5d}   contig {c:5d}   lost {g - c:4d}'
              f'   ({100 * (g - c) / g:.0f}% of calls)' if g else f'  {st}: none')

    if contigs:
        print('\nloss against assembly fragmentation:')
        bins = [(1, 1), (2, 10), (11, 50), (51, 200), (201, 10 ** 9)]
        for lo, hi in bins:
            sub = [a for a in common if lo <= contigs.get(a, 0) <= hi]
            if not sub:
                continue
            g = sum(1 for a in sub if routed(gen, a))
            c = sum(1 for a in sub if routed(ctg, a))
            name = f'{lo}' if lo == hi else (f'{lo}-{hi}' if hi < 10 ** 9 else f'{lo}+')
            pct = f'{100 * (g - c) / g:.0f}%' if g else 'n/a'
            print(f'  {name:>8s} contigs  n={len(sub):5d}  routed {g:5d} -> {c:5d}  lost {pct}')

    print(f'\ntotal gate calls lost across all gates: {tot_lost}')


if __name__ == '__main__':
    main()
