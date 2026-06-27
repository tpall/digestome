#!/usr/bin/env python3
"""Score AD/methanogenesis pathway completeness from hmmsearch (+ optional dbCAN)
output against the KEGG-less marker-gene panel. Part of DRAM feature/kegg-less-ad.

Inputs:
  --tblout   hmmsearch --tblout produced with db/ad_panel.hmm (run hmmsearch with
             --cut_nc for NCBIfam trusted cutoffs, or pass --evalue here).
  --map      db/ad_panel_map.tsv from build_ad_hmm_db.sh (model -> gene/module).
  --panel    AD_methanogenesis_panel.tsv (module / branch / tier definitions).
  --name     genome / MAG identifier for the output.
Outputs:
  --out      per-module completeness TSV.
  stdout     branch-gate + diagnostic-marker summary.

Stdlib only. SCAFFOLD — review thresholds and gate logic before production.
Run once per genome/MAG; loop externally for a community profile.
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

def load_panel(path):
    modules = {}
    for r in csv.DictReader(open(path), delimiter='\t'):
        m = r['module']
        d = modules.setdefault(m, {'branch': r['branch'], 'core': set(), 'accessory': set()})
        tier = r['tier'] if r['tier'] in ('core', 'accessory') else 'accessory'
        for g in gene_tokens(r['gene']):
            d[tier].add(g)
    return modules

def load_map(path):
    acc2gene, name2gene = {}, {}
    for r in csv.reader(open(path), delimiter='\t'):
        if not r or r[0].startswith('#'):
            continue
        acc, name, gene = (r + ['', '', ''])[:3]
        if acc:  acc2gene[acc]   = gene
        if name: name2gene[name] = gene
    return acc2gene, name2gene

def parse_tblout(path, evalue):
    hits = set()
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f = line.split()
        if len(f) < 6:
            continue
        qname, qacc = f[2], f[3]            # --tblout: query = the HMM model
        try:
            ev = float(f[4])
        except ValueError:
            ev = 0.0
        if evalue is not None and ev > evalue:
            continue
        hits.add((qacc, qname))
    return hits

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tblout', required=True)
    ap.add_argument('--map', required=True)
    ap.add_argument('--panel', required=True)
    ap.add_argument('--name', default='genome')
    ap.add_argument('--out', required=True)
    ap.add_argument('--dbcan', help='optional run_dbcan overview.txt for HYDROL-CARB')
    ap.add_argument('--evalue', type=float, default=None,
                    help='E-value cutoff; omit if hmmsearch used --cut_nc')
    a = ap.parse_args()

    modules = load_panel(a.panel)
    acc2gene, name2gene = load_map(a.map)
    hits = parse_tblout(a.tblout, a.evalue)

    found = set()
    for qacc, qname in hits:
        g = (acc2gene.get(qacc) or name2gene.get(qname)
             or name2gene.get(qacc) or acc2gene.get(qname))
        if g:
            found.add(g.lower())

    # optional CAZyme presence -> satisfies carbohydrate hydrolysis
    cazyme_families = set()
    if a.dbcan:
        for line in open(a.dbcan):
            if not line.strip() or line.startswith('Gene'):
                continue
            cazyme_families.update(re.findall(r'\b(?:GH|CE|PL|GT|CBM)\d+\b', line))
        if cazyme_families:
            found.add('cazyme')

    def has(*gs):
        return any(g.lower() in found for g in gs)

    # ---- per-module completeness ----
    rows = []
    for m, d in sorted(modules.items()):
        if m == 'HYDROL-CARB':
            fc, exp = (1 if 'cazyme' in found else 0), 1
        else:
            exp = len(d['core'])
            fc = len({g for g in d['core'] if g.lower() in found})
        pct = round(100 * fc / exp, 1) if exp else 0.0
        rows.append((m, d['branch'], exp, fc, pct))

    # ---- branch gates (own logic, replacing KEGG modules) ----
    c1 = [g for g in ('fwdB','fmdB','ftr','mch','mtd','hmd','mer','mtrA') if g.lower() in found]
    gates = {
        'methanogenesis: hydrogenotrophic': has('mcrA') and len(c1) >= 4,
        'methanogenesis: acetoclastic':     has('mcrA') and has('cdhA') and (has('acs') or (has('ackA') and has('pta'))),
        'methanogenesis: methylotrophic':   has('mcrA') and has('mtaB','mttB','mtbB','mtmB','mtsA'),
        'acetogenesis: Wood-Ljungdahl':     has('fhs') and (has('acsB') or has('cdhC')) and has('cooS','acsA'),
        'syntrophy: butyrate oxidation':    has('bcd') and has('crt') and has('hbd') and has('thlA','atoB'),
    }
    genus_hint = ('Methanothrix/Methanosaeta (acs)' if has('acs')
                  else 'Methanosarcina (ackA+pta)' if (has('ackA') and has('pta'))
                  else '-')

    # ---- write module table ----
    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['genome','module','branch','core_expected','core_found','pct_complete'])
        for m, b, exp, fc, pct in rows:
            w.writerow([a.name, m, b, exp, fc, pct])

    # ---- stdout summary ----
    print(f"# {a.name} — AD/methanogenesis KEGG-less summary")
    print(f"mcrA (master methanogen marker): {'PRESENT' if has('mcrA') else 'absent'}")
    print(f"acetoclastic genus hint: {genus_hint}")
    print("branch gates:")
    for k, v in gates.items():
        print(f"  [{'x' if v else ' '}] {k}")
    print("diagnostic markers:")
    for mk in ('mcrA','fhs','fwdB','mtrA','cdhA','acs','mttB','mtaB','hydA','frhA'):
        print(f"  {mk:6} {'+' if has(mk) else '-'}")
    if cazyme_families:
        print(f"CAZyme families (hydrolysis): {', '.join(sorted(cazyme_families))}")
    print(f"\nmodule table -> {a.out}")

if __name__ == '__main__':
    main()
