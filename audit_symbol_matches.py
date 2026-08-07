#!/usr/bin/env python3
"""Audit non-curated panel matches for gene-symbol collisions.

A gene symbol is not a unique key. `fmdA` is the panel's formylmethanofuran
dehydrogenase subunit A AND a legitimate NCBIfam symbol for formamidase, an
unrelated enzyme that was being scored as a methanogenesis marker. This finds
the rest of that class: matches where the enzyme the panel is asking for and the
product NCBIfam gives the model have nothing in common.

Heuristic, not a verdict. It scores lexical overlap between the panel's `enzyme`
text and NCBIfam's `product_name`, then ranks the weakest first for review.
Zero overlap means "look at this", not "this is wrong" -- synonyms and
abbreviations that escape the alias table will land here too.

Usage:
  audit_symbol_matches.py --panel AD_methanogenesis_panel.tsv \
      --map db/ad_panel_map.tsv --ncbifam _work/hmm_PGAP.tsv [--basis gene_symbol]

Stdlib only.
"""
import argparse, csv, re, sys

SKIP = {'multi','core','use','dbcan','see','pep','gh','lip','est'}
def gene_tokens(field):
    field = re.sub(r'\(.*?\)', '', field)
    out = []
    for p in re.split(r'[\/\+\s]+', field):
        for q in p.split('-'):
            q = q.strip()
            if len(q) >= 3 and re.search(r'[a-z]', q) and q.lower() not in SKIP:
                out.append(q)
    return out

# Words that carry no discriminating signal between two enzyme descriptions.
STOP = {'the','and','of','a','an','for','with','type','family','protein','subunit',
        'putative','related','domain','containing','system','chain','component',
        'alpha','beta','gamma','delta','epsilon','large','small','group'}

# Abbreviations the panel uses that NCBIfam spells out (and vice versa).
ALIAS = {
    'h4mpt': 'tetrahydromethanopterin', 'thf': 'tetrahydrofolate',
    'mfr': 'methanofuran', 'dep': 'dependent', 'dehydrog': 'dehydrogenase',
    'com': 'coenzyme', 'coa': 'coenzyme', 'red': 'reductase',
    'meth': 'methyl', 'synth': 'synthase', 'transf': 'transferase',
}

def terms(text):
    out = set()
    for w in re.split(r'[^A-Za-z0-9]+', (text or '').lower()):
        if not w or w in STOP or len(w) < 3:
            continue
        out.add(ALIAS.get(w, w))
    return out

def overlap(a, b):
    """Shared terms, allowing prefix matches for longer words (dehydrogenase/ases)."""
    hit = set()
    for x in a:
        for y in b:
            if x == y or (len(x) >= 5 and len(y) >= 5 and (x.startswith(y) or y.startswith(x))):
                hit.add(x)
    return hit

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--panel', required=True)
    ap.add_argument('--map', required=True)
    ap.add_argument('--ncbifam', required=True, help='hmm_PGAP.tsv')
    ap.add_argument('--basis', default='gene_symbol',
                    help="match_basis to audit: gene_symbol, ec, or all")
    a = ap.parse_args()

    # gene token -> enzyme text the panel is asking for
    want = {}
    for r in csv.DictReader(open(a.panel), delimiter='\t'):
        for g in gene_tokens(r['gene']):
            want.setdefault(g.lower(), set()).add(
                f"{r.get('enzyme','')} {r.get('branch','')}")

    prod, tax = {}, {}
    for r in csv.DictReader(open(a.ncbifam), delimiter='\t'):
        acc = (r.get('#ncbi_accession') or r.get('ncbi_accession') or '').strip()
        if acc:
            prod[acc] = (r.get('product_name') or '').strip()
            tax[acc]  = (r.get('taxonomic_range_name') or '').strip()

    rows = []
    for r in csv.DictReader(open(a.map), delimiter='\t'):
        basis = (r.get('match_basis') or '').strip()
        keep = (basis == a.basis if a.basis in ('gene_symbol','curated')
                else basis.startswith('ec:') if a.basis == 'ec'
                else basis != 'curated')
        if not keep:
            continue
        acc, gene = r['#model_accession'], (r.get('gene') or '')
        p = prod.get(acc, '')
        shared = overlap(terms(' '.join(want.get(gene.lower(), set()))), terms(p))
        rows.append((len(shared), gene, acc, basis, p, tax.get(acc, ''), sorted(shared)))

    rows.sort(key=lambda t: (t[0], t[1]))
    sus = [r for r in rows if r[0] == 0]
    weak = [r for r in rows if r[0] == 1]

    print(f"audited {len(rows)} '{a.basis}' matches")
    print(f"  {len(sus)} with NO shared term  <-- review these")
    print(f"  {len(weak)} with a single shared term")
    print()
    for label, group in (('NO OVERLAP', sus), ('WEAK (1 shared term)', weak)):
        if not group:
            continue
        print(f"=== {label} ===")
        for n, gene, acc, basis, p, tx, shared in group:
            print(f"  {gene:<10} {acc:<13} {p}")
            print(f"  {'':<10} {'':<13} panel wants: "
                  f"{' | '.join(sorted(want.get(gene.lower(), {'?'})))}")
            if tx:
                print(f"  {'':<10} {'':<13} ncbifam range: {tx}")
            if shared:
                print(f"  {'':<10} {'':<13} shared: {', '.join(shared)}")
            print()

if __name__ == '__main__':
    main()
