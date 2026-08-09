#!/usr/bin/env python3
"""Aggregate per-MAG panel_scored.py output into one digester community profile.

panel_scored.py answers "what can this genome do". A digester report has to
answer "what can this community do, and what is at risk" -- which is a different
question, because the operationally interesting facts are distributional: which
methanogenesis route dominates, whether acetate has a consumer, whether anything
can oxidise propionate.

Inputs are the artifacts panel_scored.py already writes, one set per MAG:
  <name>.modules.tsv   per-module completeness
  <name>.summary.txt   gates, diagnostic markers, taxonomy

Emits a JSON profile (for rendering) and a human-readable text summary.

Design rule inherited from the rest of this branch: never let "not measured"
render as "absent". Every count here carries its denominator, and markers with
no detector are reported separately rather than folded into a zero.

Stdlib only.
"""
import argparse, json, os, re, sys

# Operational reading of route balance in an anaerobic digester. Acetate is the
# precursor of ~70% of methane in a healthy reactor, so losing its consumer is
# the classic precursor to VFA accumulation and souring.
ROUTE_GATES = {
    'acetoclastic':      'methanogenesis: acetoclastic',
    'hydrogenotrophic':  'methanogenesis: hydrogenotrophic',
    'methylotrophic':    'methanogenesis: methylotrophic',
}

def parse_summary(path):
    gates, markers, hint, lineage = {}, {}, '', ''
    section = None
    for line in open(path):
        s = line.rstrip('\n')
        if s.startswith('taxonomy:'):
            lineage = s.split(':', 1)[1].strip()
        elif s.startswith('acetoclastic genus hint:'):
            hint = s.split(':', 1)[1].strip()
        if s.startswith('branch gates:'):
            section = 'gates'; continue
        if s.startswith('diagnostic markers:'):
            section = 'markers'; continue
        if s.startswith(('CAZyme', 'not scored', 'module table')):
            section = None
        if section == 'gates':
            m = re.match(r'\s*\[([x? ])\]\s+(.*?)(?:\s+—\s+(.*))?$', s)
            if m:
                flag, label, note = m.group(1), m.group(2).strip(), (m.group(3) or '')
                gates[label] = {'value': None if flag == '?' else (flag == 'x'),
                                'note': note.strip()}
        elif section == 'markers':
            m = re.match(r'\s*(\S+)\s+([+-])\s*$', s)
            if m:
                markers[m.group(1)] = (m.group(2) == '+')
    return gates, markers, hint, lineage

def parse_modules(path):
    rows, hdr = [], None
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if hdr is None:
            hdr = f; continue
        rows.append(dict(zip(hdr, f)))
    return rows

def taxon(lineage, rank='g__'):
    for part in (lineage or '').split(';'):
        part = part.strip()
        if part.startswith(rank) and len(part) > len(rank):
            return part[len(rank):]
    return ''

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', required=True, help='directory of *.summary.txt / *.modules.tsv')
    ap.add_argument('--sample', default='sample', help='sample / digester name')
    ap.add_argument('--description', default='', help='what this sample actually is')
    ap.add_argument('--out-json', required=True)
    ap.add_argument('--out-txt')
    a = ap.parse_args()

    names = sorted(f[:-len('.summary.txt')] for f in os.listdir(a.dir)
                   if f.endswith('.summary.txt'))
    if not names:
        print(f"!! no *.summary.txt in {a.dir}", file=sys.stderr)
        return 2

    mags = []
    for n in names:
        gates, markers, hint, lineage = parse_summary(os.path.join(a.dir, n + '.summary.txt'))
        mod_path = os.path.join(a.dir, n + '.modules.tsv')
        modules = parse_modules(mod_path) if os.path.exists(mod_path) else []
        mags.append({
            'name': n, 'lineage': lineage,
            'genus': taxon(lineage, 'g__'), 'family': taxon(lineage, 'f__'),
            'phylum': taxon(lineage, 'p__'),
            'is_methanogen': bool(markers.get('mcrA')),
            'gates': gates, 'markers': markers, 'genus_hint': hint,
            'modules': modules,
        })

    n_tot = len(mags)
    methanogens = [m for m in mags if m['is_methanogen']]

    # ---- route inventory ----
    routes = {}
    for key, label in ROUTE_GATES.items():
        present, unassessable = [], []
        for m in mags:
            g = m['gates'].get(label)
            if not g:
                continue
            if g['value'] is None:
                unassessable.append(m['name'])
            elif g['value']:
                present.append({'mag': m['name'],
                                'genus': m['genus'],
                                'note': g['note']})
        routes[key] = {'present_in': present,
                       'n_present': len(present),
                       'unassessable_in': unassessable}

    # ---- other capabilities ----
    def gate_carriers(label):
        return [m['name'] for m in mags
                if (m['gates'].get(label) or {}).get('value') is True]
    def gate_blind(label):
        return [m['name'] for m in mags
                if (m['gates'].get(label) or {}).get('value') is None]

    capabilities = {}
    for key, label in [('wood_ljungdahl', 'acetogenesis: Wood-Ljungdahl'),
                       ('butyrate_oxidation', 'syntrophy: butyrate oxidation'),
                       ('ethanolamine', 'substrate: ethanolamine (NH3 source)')]:
        capabilities[key] = {'label': label,
                             'carriers': gate_carriers(label),
                             'unassessable_in': gate_blind(label)}

    # ---- module completeness across the community ----
    # Reported as "best MAG" and "how many MAGs carry it", never as a mean --
    # a community function only needs one organism to perform it.
    mod_stats = {}
    for m in mags:
        for r in m['modules']:
            mod = r.get('module')
            if not mod:
                continue
            st = mod_stats.setdefault(mod, {
                'branch': r.get('branch', ''), 'best_pct': None,
                'best_mag': '', 'n_scored': 0, 'n_unscoreable': 0,
                'statuses': set()})
            st['statuses'].add(r.get('status', ''))
            pct = r.get('pct_complete', 'NA')
            if pct == 'NA' or pct == '':
                st['n_unscoreable'] += 1
                continue
            st['n_scored'] += 1
            try: v = float(pct)
            except ValueError: continue
            if st['best_pct'] is None or v > st['best_pct']:
                st['best_pct'] = v; st['best_mag'] = m['name']
    for st in mod_stats.values():
        st['statuses'] = sorted(x for x in st['statuses'] if x)

    # ---- risk reading ----
    # Deliberately conservative: these are capability statements about DNA, not
    # measurements of activity. Wording matters in a client deliverable.
    aceto_n = routes['acetoclastic']['n_present']
    hydro_n = routes['hydrogenotrophic']['n_present']
    risks = []
    if methanogens and aceto_n == 0 and hydro_n > 0:
        risks.append({
            'level': 'attention',
            'title': 'No acetoclastic methanogen detected',
            'detail': (f'{hydro_n} hydrogenotrophic methanogen(s) but no acetate-consuming '
                       'methanogen. In a working digester acetate is the precursor of most '
                       'methane, so acetate turnover would have to run through syntrophic '
                       'acetate oxidation, which is slower and more ammonia- and '
                       'temperature-sensitive. Corroborate with measured VFA before acting.')})
    if aceto_n and hydro_n:
        risks.append({
            'level': 'ok',
            'title': 'Both acetate and H2/CO2 routes represented',
            'detail': 'Methane can be produced from acetate and from H2/CO2, which is the '
                      'expected configuration for a stable reactor.'})
    if capabilities['ethanolamine']['carriers']:
        risks.append({
            'level': 'note',
            'title': 'Ethanolamine utilisation present',
            'detail': (f"{len(capabilities['ethanolamine']['carriers'])} MAG(s) can catabolise "
                       'ethanolamine, which releases ammonia stoichiometrically. Relevant to '
                       'ammonia inhibition risk in protein- and lipid-rich feedstocks. '
                       'This is genomic capacity, not a flux measurement — read alongside '
                       'measured NH4+.')})
    if not methanogens:
        risks.append({
            'level': 'attention',
            'title': 'No methanogen detected in this MAG set',
            'detail': 'No MAG carries mcrA. Either the archaeal fraction was not recovered '
                      'by binning, or this is not a methanogenic community. Absence in a MAG '
                      'set is weak evidence of absence in the reactor.'})

    profile = {
        'sample': a.sample,
        'description': a.description,
        'n_mags': n_tot,
        'n_methanogens': len(methanogens),
        'mags': mags,
        'routes': routes,
        'capabilities': capabilities,
        'modules': mod_stats,
        'risks': risks,
    }
    with open(a.out_json, 'w') as fh:
        json.dump(profile, fh, indent=2, default=list)

    lines = []
    w = lines.append
    w(f"# {a.sample} — anaerobic digestion community profile")
    if a.description:
        w(f"  {a.description}")
    w(f"  {n_tot} MAG(s); {len(methanogens)} carrying mcrA")
    w("")
    w("methanogenesis routes:")
    for key in ROUTE_GATES:
        r = routes[key]
        who = ', '.join(f"{p['mag']}" + (f" ({p['genus']})" if p['genus'] else '')
                        for p in r['present_in']) or '—'
        w(f"  {key:<18} {r['n_present']}/{n_tot}  {who}")
        if r['unassessable_in']:
            w(f"  {'':<18} not assessable in {len(r['unassessable_in'])} MAG(s)")
    w("")
    w("other capabilities:")
    for key, c in capabilities.items():
        extra = (f"  (not assessable in {len(c['unassessable_in'])})"
                 if c['unassessable_in'] else '')
        w(f"  {key:<20} {len(c['carriers'])}/{n_tot}{extra}")
    w("")
    w("risk reading:")
    for r in risks:
        w(f"  [{r['level']}] {r['title']}")
    if a.out_txt:
        open(a.out_txt, 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    return 0

if __name__ == '__main__':
    sys.exit(main())
