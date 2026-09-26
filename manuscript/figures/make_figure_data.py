#!/usr/bin/env python3
"""Reduce the scored benchmark outputs to the small tables the manuscript figures are drawn from.

    bash manuscript/figures/fetch_benchmarks.sh      # scored outputs off the cluster, gitignored
    python3 manuscript/figures/make_figure_data.py   # -> manuscript/figures/data/*.tsv, committed
    Rscript manuscript/figures/plot_figures.R        # -> manuscript/figures/fig-*.{pdf,png}

No number in a figure is typed. The same scored directories feed the tables in Results, so a figure
and the table beside it cannot disagree without one of these scripts changing.

Inputs, under results/benchmarks/ (see fetch_benchmarks.sh for where each comes from):
  gtdb/final               the GTDB sample under the lineage policy (Results, tbl-gtdb, tbl-routes)
  gtdb/family_level_aug11  the same genomes under the family-level check (tbl-routes, "before"; d2daee1)
  gtdb/contig_final        the final scorer with --contig-level, paired with gtdb/final
  gtdb/contig_counts.tsv   contigs per genome, counted from the Prodigal proteomes
  catalogue/final          the 1,401 digester MAGs (tbl-secretion)
"""
import collections
import csv
import os
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT / 'bin'))
from panel_scored import load_acetate_lineages, lineage_in_role  # noqa: E402

B = ROOT / 'results' / 'benchmarks'
GTDB, CAT = B / 'gtdb', B / 'catalogue'
AFTER, BEFORE, CONTIG = 'final', 'family_level_aug11', 'contig_final'
CATALOGUE = 'final'
OUT = HERE / 'data'

GATE = re.compile(r'\[([x ?])\]\s+methanogenesis:\s*(\S+)')
ROUTES = ('acetoclastic', 'hydrogenotrophic', 'methylotrophic')


def need(p):
    if not p.exists():
        sys.exit(f'!! {p} is missing: run manuscript/figures/fetch_benchmarks.sh first')
    return p


def routes(d):
    """-> {accession: set of methanogenesis routes marked [x]}; [?] is not assessable, never a call"""
    out = {}
    for f in need(d).glob('*.summary.txt'):
        out[f.name[:-len('.summary.txt')]] = {r for m, r in GATE.findall(f.read_text()) if m == 'x'}
    if not out:
        sys.exit(f'!! no scored genomes under {d}')
    return out


def write(name, header, rows):
    OUT.mkdir(exist_ok=True)
    with open(OUT / name, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)
    print(f'  data/{name}  ({len(rows)} rows)')


meta = {r['accession']: r for r in csv.DictReader(open(need(GTDB / 'gtdb_sample.tsv')), delimiter='\t')}
lineage = {}
for line in need(GTDB / 'gtdbtk_sample.tsv').read_text().splitlines():
    f = line.split('\t')
    if len(f) > 1 and f[0].startswith('GC'):
        lineage[f[0]] = f[1]
policy = load_acetate_lineages()


def rank(a, r):
    for tok in lineage.get(a, '').split(';'):
        if tok.startswith(r + '__'):
            return re.sub(r'_[A-Z]+$', '', tok[3:]) or '?'
    return '?'


# --- Figure: route calls by lineage, before and after lineage confirmation -----------------------
# Lineage is the level each gate is confirmed at: genus for the acetoclastic route, order for the
# other two. Anaerobic methane and alkane oxidisers sit inside methanogen orders, so they are pulled
# out as their own group: that is the error the policy removes, and hiding it inside
# Methanosarcinales would hide the result.
before, after = routes(GTDB / BEFORE), routes(GTDB / AFTER)
common = sorted(set(before) & set(after))


RELABEL_AS_METHANOGEN = {'o__Methanonatronarchaeales'}   # as in scripts/summarise_gtdb_benchmark.py


def klass(a):
    """The benchmark's three classes: oxidisers by lineage policy, one sourced relabel, else the stratum."""
    if lineage.get(a) and lineage_in_role(lineage[a], 'methane_oxidiser', policy):
        return 'methane_oxidiser'
    st = meta.get(a, {}).get('stratum', '?')
    if st == 'background' and meta.get(a, {}).get('order') in RELABEL_AS_METHANOGEN:
        return 'methanogen'
    return st


def group(a, route):
    if lineage.get(a) and lineage_in_role(lineage[a], 'methane_oxidiser', policy):
        return 'Anaerobic methane / alkane oxidisers'
    if route == 'acetoclastic':
        # GTDB placeholder genera (DQIP01, Fen-7) are pooled by family: the policy treats the two
        # families differently, and one row per placeholder hides that.
        g = rank(a, 'g')
        return f'unnamed {rank(a, "f")} genera' if re.search(r'\d', g) else g
    return rank(a, 'o')


rows = []
for route in ROUTES:
    tally = collections.defaultdict(lambda: [0, 0])
    for a in common:
        if route in before[a]:
            tally[group(a, route)][0] += 1
        if route in after[a]:
            tally[group(a, route)][1] += 1
    for g, (b, k) in sorted(tally.items(), key=lambda kv: -kv[1][0]):
        rows.append((route, g, b, k))
write('route_calls_by_lineage.tsv', ('route', 'lineage', 'calls_family_check', 'calls_lineage_policy'), rows)

# --- Figure: what --contig-level costs, against assembly fragmentation ---------------------------
# Paired with the final call set, and counted over the same methanogen class as tbl-gtdb, so the
# genome-level totals here are the paper's methanogens with a route.
contig = routes(GTDB / CONTIG)
ncontig = {}
for line in need(GTDB / 'contig_counts.tsv').read_text().splitlines():
    a, n = line.split('\t')
    ncontig[a] = int(n)
BINS = [('1', 1, 1), ('2–10', 2, 10), ('11–50', 11, 50), ('51–200', 51, 200), ('201+', 201, 10 ** 9)]
methanogens = [a for a in set(after) & set(contig) if klass(a) == 'methanogen']
rows = []
for name, lo, hi in BINS:
    sub = [a for a in methanogens if lo <= ncontig.get(a, 0) <= hi]
    g = sum(1 for a in sub if after[a])
    c = sum(1 for a in sub if contig[a])
    rows.append((name, len(sub), g, c))
write('contig_level_loss.tsv', ('contigs', 'genomes', 'routed_genome_level', 'routed_contig_level'), rows)

# --- Figure: carries the family versus exports it, 1,401 digester MAGs ---------------------------
carry, export = collections.Counter(), collections.Counter()
stage = {}
for f in need(CAT / CATALOGUE).glob('*.modules.tsv'):
    for r in csv.DictReader(open(f), delimiter='\t'):
        if r['core_found'] not in ('', 'NA') and int(r['core_found']) > 0:
            carry[r['module']] += 1
        if r['core_secreted'] not in ('', 'NA') and int(r['core_secreted']) > 0:
            export[r['module']] += 1
for r in csv.DictReader(open(ROOT / 'assets' / 'AD_methanogenesis_panel.tsv'), delimiter='\t'):
    stage.setdefault(r['module'], r['stage'])
rows = [(m, stage.get(m, '?'), carry[m], export[m]) for m in sorted(carry, key=lambda m: -carry[m])]
write('secretion_by_module.tsv', ('module', 'stage', 'genomes_carrying', 'genomes_exporting'), rows)

print('figure data written from the scored outputs, nothing retyped')
