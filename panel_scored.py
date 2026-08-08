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

Identity note: everything joins on the panel's `marker_id` (e.g. MEG-CORE:mcrA),
never on the gene name. Gene names are 3-4 letter labels that are unique neither
within the panel nor in NCBIfam, and using one as a key hid real errors -- EC
matches collapsed every model sharing an EC onto one arbitrary name, and an
accession serving two panel rows could record only one of them. Gene tokens
survive only as display text and as an input to `has()`, which resolves them
through the panel so gate expressions stay readable.

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
    """module -> {branch, core, accessory, scoring}.

    `scoring` comes from the panel's own column so detection policy lives with
    the data, not in module names hardcoded here:
      scored   normal HMM scoring off the core tier (default)
      assumed  universal central metabolism -- present by assumption, no marker
      dbcan    satisfied by any CAZyme call from run_dbcan
      pending  no license-clean detector wired yet; reported, never scored
    """
    modules, by_gene = {}, {}
    for r in csv.DictReader(open(path), delimiter='\t'):
        m = r['module']
        mid = (r.get('marker_id') or '').strip() or \
              f"{m}:{(gene_tokens(r['gene']) or [r['gene']])[0]}"
        d = modules.setdefault(m, {'branch': r['branch'], 'core': set(),
                                   'accessory': set(), 'scoring': None})
        tier = r['tier'] if r['tier'] in ('core', 'accessory') else 'accessory'
        # Tiers hold marker_ids, not gene names -- see the module docstring.
        d[tier].add(mid)
        for g in gene_tokens(r['gene']):
            by_gene.setdefault(g.lower(), set()).add(mid)
        s = (r.get('scoring') or '').strip().lower()
        if s:
            if d['scoring'] and d['scoring'] != s:
                print(f"!! {m}: conflicting scoring modes "
                      f"'{d['scoring']}' vs '{s}' -- keeping the first",
                      file=sys.stderr)
            elif s not in ('scored', 'assumed', 'dbcan', 'pending'):
                print(f"!! {m}: unknown scoring mode '{s}' -- treating as 'scored'",
                      file=sys.stderr)
                s = 'scored'
            d['scoring'] = d['scoring'] or s
    for d in modules.values():
        d['scoring'] = d['scoring'] or 'scored'   # panels without the column
    return modules, by_gene

def load_map(path):
    """accession/model-name -> set(marker_id).

    Sets, not single values: one model can legitimately serve two panel rows
    (the same methylmalonyl-CoA mutase marks both ACID-PROP and ACET-SYN-PROP).
    Storing one gene per accession meant the second row scored 0 with its model
    sitting in the DB.
    """
    acc2mid, name2mid = {}, {}
    for r in csv.reader(open(path), delimiter='\t'):
        if not r or r[0].startswith('#'):
            continue
        acc, name, mid = (r + ['', '', ''])[:3]
        if not mid:
            continue
        if acc:  acc2mid.setdefault(acc, set()).add(mid)
        if name: name2mid.setdefault(name, set()).add(mid)
    return acc2mid, name2mid

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

    modules, by_gene = load_panel(a.panel)
    acc2mid, name2mid = load_map(a.map)
    hits = parse_tblout(a.tblout, a.evalue)

    found = set()                                  # marker_ids, not gene names
    for qacc, qname in hits:
        for src in (acc2mid.get(qacc), name2mid.get(qname),
                    name2mid.get(qacc), acc2mid.get(qname)):
            if src:
                found |= src
                break

    # Markers that COULD be hit, i.e. have at least one model in the pressed DB.
    # Distinguishing "looked for and absent" from "never looked for" is the whole
    # point -- a marker with no detector must never read as a negative result.
    detectable = set()
    for d in (acc2mid, name2mid):
        for mids in d.values():
            detectable |= mids

    # optional CAZyme presence -> satisfies carbohydrate hydrolysis
    cazyme_families = set()
    if a.dbcan:
        for line in open(a.dbcan):
            if not line.strip() or line.startswith('Gene'):
                continue
            cazyme_families.update(re.findall(r'\b(?:GH|CE|PL|GT|CBM)\d+\b', line))

    def has(*names):
        """True if any named marker was found.

        Accepts a marker_id ('MEG-CORE:mcrA') or a bare gene token ('mcrA').
        Tokens are resolved through the panel, so a gate stays readable while
        the underlying join is still on marker_id -- and a token shared by two
        modules is satisfied by either.
        """
        for n in names:
            if n in found or any(m in found for m in by_gene.get(n.lower(), ())):
                return True
        return False

    # ---- per-module completeness ----
    # A module with nothing to score reports pct = NA, never 0.0. Emitting 0.0
    # made "no marker defined" indistinguishable from "pathway absent" -- on a
    # digester report that reads as a failing plant.
    rows = []
    for m, d in sorted(modules.items()):
        mode = d['scoring']
        if mode == 'assumed':
            exp = fc = 0
            status = 'assumed-present'
        elif mode == 'dbcan':
            if a.dbcan:
                exp, fc = 1, (1 if cazyme_families else 0)
                status = 'scored'
            else:
                exp = fc = 0
                status = 'no-detector'        # run_dbcan output not supplied
        elif mode == 'pending':
            exp = fc = 0
            status = 'not-wired'
        else:
            exp = len(d['core'])
            fc = len(d['core'] & found)
            # Empty core tier: the module's genes are all accessory-tier, which
            # scoring does not read yet. Not a zero -- an absence of a metric.
            status = 'scored' if exp else 'no-core-tier'
        # How much of the core has a model at all: pct must be read against
        # this, not against exp. MEG-CORE can never exceed 10/11 while mvhA has
        # no detector, and 90.9% should not read as a missing gene.
        det = len(d['core'] & detectable)
        if status == 'scored' and exp and not det:
            # Module-level twin of the dead gate: every core marker lacks a
            # model, so 0.0% would be a measurement that never happened.
            exp = fc = 0
            status = 'no-detector-for-core'
        pct = round(100 * fc / exp, 1) if exp else None
        rows.append((m, d['branch'], exp, det, fc, pct, status))

    # ---- branch gates (own logic, replacing KEGG modules) ----
    # A gate whose requirements include a marker with NO model in the DB cannot
    # be evaluated at all, and printing it as an unticked box is a guaranteed
    # false negative rather than a cautious one. Butyrate oxidation was exactly
    # this: bcd and crt have no NCBIfam model, so two of its four AND-terms were
    # permanently false and the branch reported absent for every sample ever run.
    # Same rule as the module table -- never report absence when nothing was
    # measured. Requirements below mirror each boolean's structure: groups are
    # ANDed, members within a group are ORed.
    def _markers(n):
        return by_gene.get(n.lower()) or {n}

    def unmeasurable(groups):
        """Requirement groups with no detectable member -- the gate is blind."""
        return [g for g in groups
                if not any(_markers(n) & detectable for n in g)]

    C1 = ('fwdB','fmdB','ftr','mch','mtd','hmd','mer','mtrA')
    c1 = [g for g in C1 if has(g)]
    gate_defs = [
        ('methanogenesis: hydrogenotrophic',
         has('mcrA') and len(c1) >= 4, [('mcrA',), C1]),
        ('methanogenesis: acetoclastic',
         has('mcrA') and has('cdhA') and (has('acs') or (has('ackA') and has('pta'))),
         [('mcrA',), ('cdhA',), ('acs','ackA'), ('acs','pta')]),
        ('methanogenesis: methylotrophic',
         has('mcrA') and has('mtaB','mttB','mtbB','mtmB','mtsA'),
         [('mcrA',), ('mtaB','mttB','mtbB','mtmB','mtsA')]),
        ('acetogenesis: Wood-Ljungdahl',
         has('fhs') and (has('acsB') or has('cdhC')) and has('cooS','acsA'),
         [('fhs',), ('acsB','cdhC'), ('cooS','acsA')]),
        ('syntrophy: butyrate oxidation',
         has('bcd') and has('crt') and has('hbd') and has('thlA','atoB'),
         [('bcd',), ('crt',), ('hbd',), ('thlA','atoB')]),
        # Substrate capability, not a methanogenesis branch. Both ammonia-lyase
        # subunits required: EutB alone is not a functional enzyme, and the
        # downstream eut genes (eutD/eutE) only differ from housekeeping
        # acetate metabolism by operon context.
        ('substrate: ethanolamine (NH3 source)',
         has('eutB') and has('eutC'), [('eutB',), ('eutC',)]),
    ]
    # (label, True/False/None-if-unmeasurable, tokens with no detector)
    gates = [(label, None if blind else value, sorted({n for g in blind for n in g}))
             for label, value, reqs in gate_defs
             for blind in (unmeasurable(reqs),)]
    # Only offer a genus hint once the acetoclastic gate is actually met. acs,
    # ackA and pta are ubiquitous acetate-metabolism genes, so computing this
    # unconditionally put a methanogen genus on organisms with no mcrA at all --
    # the smoke test labelled E. coli "Methanothrix/Methanosaeta", and both
    # Clostridium ljungdahlii and Syntrophomonas wolfei "Methanosarcina".
    aceto = dict((lbl, v) for lbl, v, _b in gates).get('methanogenesis: acetoclastic')
    genus_hint = '-'
    if aceto:
        genus_hint = ('Methanothrix/Methanosaeta (acs)' if has('acs')
                      else 'Methanosarcina (ackA+pta)' if (has('ackA') and has('pta'))
                      else '-')

    # ---- write module table ----
    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['genome','module','branch','core_expected','core_detectable',
                    'core_found','pct_complete','status'])
        for m, b, exp, det, fc, pct, status in rows:
            w.writerow([a.name, m, b, exp, det, fc,
                        'NA' if pct is None else pct, status])

    # ---- stdout summary ----
    print(f"# {a.name} — AD/methanogenesis KEGG-less summary")
    print(f"mcrA (master methanogen marker): {'PRESENT' if has('mcrA') else 'absent'}")
    print(f"acetoclastic genus hint: {genus_hint}")
    print("branch gates:")
    for label, v, blind in gates:
        mark = '?' if v is None else ('x' if v else ' ')
        note = f"  — NOT ASSESSABLE: no detector for {', '.join(blind)}" if v is None else ''
        print(f"  [{mark}] {label}{note}")
    if aceto:
        # ACDS/CODH is reversible and autotrophic hydrogenotrophs run it in the
        # synthetic direction for carbon fixation, so cdhA presence does not
        # imply acetate is being consumed for methanogenesis. Gene content cannot
        # resolve direction; Methanothermobacter trips this gate in the smoke
        # test despite being an obligate hydrogenotroph.
        print("      NB: ACDS/CODH is bidirectional and also serves autotrophic carbon")
        print("          fixation — confirm against taxonomy (Methanosarcinaceae /")
        print("          Methanotrichaceae) before reporting acetoclastic methanogenesis.")
    print("diagnostic markers:")
    for mk in ('mcrA','fhs','fwdB','mtrA','cdhA','acs','mttB','mtaB','hydA','frhA','eutB'):
        print(f"  {mk:6} {'+' if has(mk) else '-'}")
    if cazyme_families:
        print(f"CAZyme families (hydrolysis): {', '.join(sorted(cazyme_families))}")

    unscored = [(m, s) for m, _b, _e, _d, _f, pct, s in rows if pct is None]
    if unscored:
        print(f"not scored ({len(unscored)}/{len(rows)} modules) — absent from the "
              f"percentages above, NOT zero:")
        for m, s in unscored:
            print(f"  {m:<13} {s}")
    print(f"\nmodule table -> {a.out}")

if __name__ == '__main__':
    main()
