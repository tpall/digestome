#!/usr/bin/env python3
"""Draw a phylum-stratified sample of GTDB species representatives.

The digester MAG benchmark establishes specificity against genomes assembled and
binned by one group with one pipeline. This sampler builds the complementary set:
isolate-derived reference genomes, which are not binned at all, so a specificity
result that holds across both cannot be an artefact of anyone's binning.

Two strata:

  methanogens   every representative in a methanogen order, taken whole. These
                measure sensitivity, so sampling them would only add noise.
  background    everything else, allocated across phyla in proportion to the
                square root of phylum size. Proportional allocation would spend
                the whole budget on Pseudomonadota; equal allocation would give
                a 36,000-genome phylum the same weight as a singleton. The
                square root is the usual compromise, and it matters here because
                false positives are likeliest in the large, well-sampled phyla.

Selection within a phylum is by even stride over the sorted accession list, not
by a random draw: the sample is then a pure function of the input files and
reproduces exactly without carrying a seed.

Usage:
  sample_gtdb_reps.py --taxonomy-dir <release/taxonomy> --n-background 2000 \
      --out gtdb_sample.tsv
"""
import argparse
import pathlib
import sys

# GTDB orders whose members are methanogens. Methanomassiliicoccales and
# Methanofastidiosales are the obligate methylotrophs that the digester benchmark
# showed carry none of the C1 carriers; they are the reason the terminal
# methyltransferase is in the panel, so they must be in the sensitivity set.
METHANOGEN_ORDERS = {
    'o__Methanobacteriales', 'o__Methanococcales', 'o__Methanomicrobiales',
    'o__Methanosarcinales', 'o__Methanotrichales', 'o__Methanocellales',
    'o__Methanopyrales', 'o__Methanomassiliicoccales', 'o__Methanofastidiosales',
    'o__Methanomethylicales',
    # Methyl-reducing haloalkaliphilic methanogens. Omitted from the first draw,
    # which then counted a correct call on GCF_029854155.1 as a false positive.
    # The panel found it via mcrA before the label list did.
    'o__Methanonatronarchaeales',
}


def read_taxonomy(paths):
    """-> [(ncbi_accession, phylum, order, full_lineage)]. Drops the GTDB RS_/GB_ prefix."""
    rows = []
    for p in paths:
        for line in pathlib.Path(p).read_text().splitlines():
            if not line.strip():
                continue
            acc, lineage = line.split('\t')[:2]
            for prefix in ('RS_', 'GB_'):
                if acc.startswith(prefix):
                    acc = acc[len(prefix):]
            ranks = lineage.split(';')
            phylum = ranks[1] if len(ranks) > 1 else 'p__'
            order = ranks[3] if len(ranks) > 3 else 'o__'
            rows.append((acc, phylum, order, lineage))
    return rows


def stride_sample(items, k):
    """k items spread evenly over a sorted list. Deterministic, no seed."""
    n = len(items)
    if k >= n:
        return list(items)
    return [items[(i * n) // k] for i in range(k)]


def allocate(sizes, budget):
    """Allocate `budget` across strata in proportion to sqrt(size)."""
    weights = {k: v ** 0.5 for k, v in sizes.items()}
    total = sum(weights.values())
    alloc = {k: max(1, int(budget * w / total)) for k, w in weights.items()}
    # Hand back what over-allocation to small strata borrowed, largest first.
    for k in sorted(alloc, key=lambda k: -sizes[k]):
        if sum(alloc.values()) <= budget:
            break
        alloc[k] = max(1, alloc[k] - (sum(alloc.values()) - budget))
    return {k: min(v, sizes[k]) for k, v in alloc.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--taxonomy-dir', required=True,
                    help='GTDB release taxonomy/ holding *_taxonomy_*_reps.tsv')
    ap.add_argument('--n-background', type=int, default=2000,
                    help='non-methanogen genomes to draw (default 2000)')
    ap.add_argument('--out', required=True, help='TSV: accession, phylum, order, lineage, stratum')
    args = ap.parse_args()

    tdir = pathlib.Path(args.taxonomy_dir)
    files = sorted(tdir.glob('*_taxonomy_*_reps.tsv'))
    if not files:
        sys.exit(f'no *_taxonomy_*_reps.tsv under {tdir}')
    rows = read_taxonomy(files)
    print(f'  {len(rows)} representatives from {len(files)} file(s)', file=sys.stderr)

    methanogens = sorted([r for r in rows if r[2] in METHANOGEN_ORDERS])
    background = [r for r in rows if r[2] not in METHANOGEN_ORDERS]

    by_phylum = {}
    for r in background:
        by_phylum.setdefault(r[1], []).append(r)
    for v in by_phylum.values():
        v.sort()

    sizes = {k: len(v) for k, v in by_phylum.items()}
    alloc = allocate(sizes, args.n_background)
    picked = []
    for phylum, items in sorted(by_phylum.items()):
        picked += stride_sample(items, alloc[phylum])

    # GTDB representatives are one per species, but guard against an accession
    # appearing in both domain files.
    seen, out = set(), []
    for r, stratum in [(m, 'methanogen') for m in methanogens] + [(b, 'background') for b in picked]:
        if r[0] in seen:
            continue
        seen.add(r[0])
        out.append((r[0], r[1], r[2], r[3], stratum))

    with open(args.out, 'w') as fh:
        fh.write('accession\tphylum\torder\tlineage\tstratum\n')
        for row in out:
            fh.write('\t'.join(row) + '\n')

    print(f'  methanogens (all):      {len(methanogens)}', file=sys.stderr)
    print(f'  background (sampled):   {len(picked)} from {len(by_phylum)} phyla', file=sys.stderr)
    print(f'  total written:          {len(out)} -> {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
