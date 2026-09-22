#!/usr/bin/env python3
"""Which lineages confirm an acetoclastic call. No databases needed.

    python3 tests/test_acetoclastic_lineage.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bin'))
from panel_scored import acetoclastic_lineage

A = 'd__Archaea;p__Halobacteriota;c__Methanosarcinia;'
CASES = [
    # confirmed: the acetate-using genera, under old and new names, with GTDB suffixes
    (A + 'o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina barkeri', True),
    (A + 'o__Methanosarcinales;f__Methanosarcinaceae;g__Methanosarcina_A;s__Methanosarcina_A sp1', True),
    (A + 'o__Methanotrichales;f__Methanotrichaceae;g__Methanothrix;s__Methanothrix soehngenii', True),
    ('d__Archaea;g__Methanosaeta;s__Methanosaeta concilii', True),          # NCBI-style, genus only
    # confirmed: every described Methanotrichaceae is an obligate acetoclast
    (A + 'o__Methanotrichales;f__Methanotrichaceae;g__JAAYUN01;s__JAAYUN01 sp1', True),
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
print(f'{len(CASES) - len(fails)}/{len(CASES)} passed')
sys.exit(1 if fails else 0)
