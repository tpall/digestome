#!/usr/bin/env python3
"""Build a name -> order/family table from a GTDB taxonomy file, for coarse lineages.

Catalogue taxonomies derived from NCBI organism names carry a genus only, and often an order or
family name in the genus slot (g__Thermoplasmatales). panel_scored.py --genus-orders uses this table
to place such lineages at order level, so the lineage policy (assets/acetate_lineages.tsv) applies to
the catalogue benchmark as it does to GTDB-Tk output. Genus suffixes (_A, _B) are merged; a genus
that maps to more than one order is left out (ambiguous). Only genus names are mapped: names above
genus mean different things in NCBI and GTDB (NCBI "Thermoplasmatales" held what GTDB calls
Methanomassiliicoccales), so a coarse lineage with a higher-rank name stays unplaced.

    build_genus_orders.py gtdb_taxonomy.tsv > genus_orders.tsv
"""
import collections, re, sys

orders = collections.defaultdict(set)
families = collections.defaultdict(set)
for line in open(sys.argv[1]):
    ranks = dict(t.split('__', 1) for t in line.rstrip('\n').split('\t')[-1].split(';') if '__' in t)
    o, f = ranks.get('o', ''), ranks.get('f', '')
    strip = lambda x: re.sub(r'_[A-Z]+$', '', x)
    g = ranks.get('g', '')
    if g:
        orders[strip(g)].add(strip(o)); families[strip(g)].add(strip(f))
print('#\tbuilt from', sys.argv[1])
print('name\torder\tfamily')
for name in sorted(orders):
    if len(orders[name]) == 1:
        fam = next(iter(families[name])) if len(families[name]) == 1 else ''
        o = next(iter(orders[name]))
        print(f'{name}\t{o}\t{fam}')
