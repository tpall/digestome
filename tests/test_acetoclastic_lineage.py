#!/usr/bin/env python3
"""Which lineages confirm an acetoclastic call. No databases needed.

    python3 tests/test_acetoclastic_lineage.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from panel_scored import acetoclastic_lineage, not_hydrogenotrophic_lineage, lineage_in_role, load_acetate_lineages

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

# hydrogenotrophic deny list: oxidative C1 users are excluded, true hydrogenotrophs and Methanosarcina are not
H2_CASES = [
    (A + 'o__Methanotrichales;f__Methanotrichaceae;g__Methanothrix;s__Methanothrix soehngenii', True),
    (A + 'o__Methanotrichales;f__Methanotrichaceae;g__Methanocrinis;s__Methanocrinis harundinaceus', True),
    (A + 'o__Methanosarcinales;f__Methanosarcinaceae;g__Methanolobus;s__Methanolobus tindarius', True),
    ('d__Archaea;g__Methanosaeta;s__Methanosaeta concilii', True),
    (A + 'o__Methanotrichales;f__Methanotrichaceae;g__UBA204;s__UBA204 sp1', True),
    (A + 'o__Methanosarcinales;f__Methanoperedenaceae;g__Methanoperedens;s__Methanoperedens nitroreducens', True),
    (A + 'o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina barkeri', False),
    ('d__Archaea;p__Halobacteriota;c__Methanomicrobia;o__Methanomicrobiales;f__Methanoculleaceae;'
     'g__Methanoculleus;s__Methanoculleus bourgensis', False),
    ('d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanothermobacteraceae;'
     'g__Methanothermobacter;s__Methanothermobacter thermautotrophicus', False),
]
h2_fails = [(l, want) for l, want in H2_CASES if not_hydrogenotrophic_lineage(l) != want]
for l, want in h2_fails:
    print(f'FAIL hydrogenotrophic deny list, expected {want}: {l}')
fails += h2_fails

# comMT-only methyl calls: accepted only in listed methylotroph lineages
M_CASES = [
    ('d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;'
     'g__Methanosphaera;s__Methanosphaera stadtmanae', True),
    ('d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;f__Methanomethylophilaceae;'
     'g__Methanoprimaticola;s__x', True),
    ('d__Archaea;p__Methanobacteriota_A;c__Methanofastidiosia;o__Methanofastidiosales;f__Methanofastidiosaceae;'
     'g__Methanofastidiosum;s__x', True),
    ('d__Archaea;p__Methanobacteriota;c__Methanococci;o__Methanococcales;f__Methanococcaceae;'
     'g__Methanococcus;s__Methanococcus maripaludis', False),
    ('d__Archaea;p__Halobacteriota;c__Methanomicrobia;o__Methanomicrobiales;f__Methanoculleaceae;'
     'g__Methanoculleus;s__Methanoculleus bourgensis', False),
    ('d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;'
     'g__Methanobacterium;s__x', False),
]
m_fails = [(l, want) for l, want in M_CASES if lineage_in_role(l, 'methyl_via_comMT') != want]
for l, want in m_fails:
    print(f'FAIL methyl_via_comMT, expected {want}: {l}')
fails += m_fails
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
N = len(CASES) + len(H2_CASES) + len(M_CASES)
print(f'{N - len([f for f in fails if f[1] in (True, False)])}/{N} lineage cases passed')
sys.exit(1 if fails else 0)
