#!/usr/bin/env python3
"""Which lineages confirm an acetoclastic call. No databases needed.

    python3 tests/test_acetoclastic_lineage.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from panel_scored import acetoclastic_lineage, load_acetate_lineages

A = 'd__Archaea;p__Halobacteriota;c__Methanosarcinia;'
CASES = [
    # confirmed: the acetate-using genera, under old and new names, with GTDB suffixes
    (A + 'o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina barkeri', True),
    (A + 'o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina_A;s__Methanosarcina_A sp1', True),
    (A + 'o__Methanotrichales;f__Methanotrichaceae;g__Methanothrix;s__Methanothrix soehngenii', True),
    ('d__Archaea;g__Methanosaeta;s__Methanosaeta concilii', True),          # NCBI-style, genus only
    (A + 'o__Methanotrichales;f__Methanotrichaceae;g__Methanothrix_B;s__Methanothrix_B sp1', True),
    (A + 'o__Methanotrichales;f__Methanotrichaceae;g__Methanocrinis;s__Methanocrinis harundinaceus', True),
    # not confirmed: an unnamed genus is left for review even in Methanotrichaceae (fail closed)
    (A + 'o__Methanotrichales;f__Methanotrichaceae;g__JAAYUN01;s__JAAYUN01 sp1', False),
    # not confirmed: methylotrophic Methanosarcinaceae genera
    (A + 'o__Methanosarcinales;f__Methanosarcinaceae;g__Methanolobus;s__Methanolobus tindarius', False),
    (A + 'o__Methanosarcinales;f__Methanosarcinaceae;g__Methanococcoides;s__Methanococcoides burtonii', False),
    (A + 'o__Methanosarcinales;f__Methanosarcinaceae;g__Methanohalophilus;s__Methanohalophilus mahii', False),
    # not confirmed: an unnamed Methanosarcinaceae genus is left for review
    (A + 'o__Methanosarcinales;f__Methanosarcinaceae;g__UBA123;s__UBA123 sp1', False),
    # not confirmed: hydrogenotrophs that carry ACDS for carbon fixation
    ('d__Archaea;p__Halobacteriota;c__Methanomicrobia;o__Methanomicrobiales;f__Methanoculleaceae;'
     'g__Methanoculleus;s__Methanoculleus bourgensis', False),
    ('d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanothermobacteraceae;'
     'g__Methanothermobacter;s__Methanothermobacter thermautotrophicus', False),
    # a genus name that only starts like an acetoclast is not one
    (A + 'o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcinales_X;s__x', False),
]

fails = [(l, want) for l, want in CASES if acetoclastic_lineage(l) != want]
for l, want in fails:
    print(f'FAIL expected {want}: {l}')

# Every name in the policy table must exist verbatim in the GTDB taxonomy in use,
# or an exact match silently fails (e.g. an order split into Methanomicrobiales_A).
# Set GTDB_TAXONOMY to gtdb_taxonomy.tsv (or an ar53 summary) to run this part.
tax = os.environ.get('GTDB_TAXONOMY')
if tax:
    seen = set()
    with open(tax) as fh:
        for line in fh:
            for tok in line.rstrip('\n').split('\t')[-1].split(';'):
                if '__' in tok:
                    seen.add((tok[0], tok.split('__', 1)[1]))
    roles = load_acetate_lineages()
    for role, names in roles.items():
        for rank, name in sorted(names):
            if (rank, name) not in seen and name != 'Methanosaeta':   # Methanosaeta: NCBI-only name
                fails.append((f'{rank}__{name}', f'in {os.path.basename(tax)}'))
                print(f'FAIL {rank}__{name} ({role}) not found verbatim in {tax}')
    print(f'policy names checked against {tax}')
print(f'{len(CASES) - len([f for f in fails if f[1] in (True, False)])}/{len(CASES)} lineage cases passed')
sys.exit(1 if fails else 0)
