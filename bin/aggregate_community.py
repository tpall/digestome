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

Emits a JSON profile (for rendering), a human-readable text summary and,
with --out-tables, four flat tables that are the client-facing "result tables":
genomes (one row per MAG), modules (one row per panel module), routes (one row
per route or capability) and risks. They are derived from the same in-memory
profile as the JSON, so a report built on them cannot drift from it.

With --abundance (a manifest from scripts/subset_catalogue_plant.py, or any TSV
with an accession/name column and rel_abundance_pct) each route and capability
is also reported as the share of mapped reads carried by the MAGs that have it.
Counts say how many organisms can do something; abundance says how much of the
community they are. A single abundant acetoclast and five rare ones are
different reactors. Abundance is reported alongside counts, never instead of
them, and the share of the community that was never scored is stated.

Design rule inherited from the rest of this branch: never let "not measured"
render as "absent". Every count here carries its denominator, and markers with
no detector are reported separately rather than folded into a zero.

Stdlib only.
"""
import argparse, csv, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from panel_scored import ACETATE_LINEAGES, load_acetate_lineages, lineage_in_role

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

def load_abundance(path):
    """-> ({name: pct}, {header key: value}) from a manifest / TSV.

    Header comments of the form '# key: value' are kept (subset_catalogue_plant.py
    writes the community totals there). The name column is 'accession' or
    'name'; the value column is 'rel_abundance_pct'."""
    ab, meta, hdr = {}, {}, None
    for line in open(path):
        s = line.rstrip('\n')
        if s.startswith('#'):
            m = re.match(r'#\s*([\w_]+):\s*(.*)$', s)
            if m:
                meta[m.group(1)] = m.group(2).strip()
            continue
        f = s.split('\t')
        if hdr is None:
            hdr = f
            try:
                ki = hdr.index('accession') if 'accession' in hdr else hdr.index('name')
                vi = hdr.index('rel_abundance_pct')
            except ValueError:
                sys.exit(f"!! {path}: need an 'accession' (or 'name') and a "
                         f"'rel_abundance_pct' column, got {hdr}")
            continue
        if len(f) <= max(ki, vi) or not f[ki]:
            continue
        try:
            ab[f[ki]] = float(f[vi])
        except ValueError:
            pass
    return ab, meta

def taxon(lineage, rank='g__'):
    for part in (lineage or '').split(';'):
        part = part.strip()
        if part.startswith(rank) and len(part) > len(rank):
            return part[len(rank):]
    return ''

# Modules whose secretion evidence is reported per genome. Hydrolysis is the
# one stage where "carries the family" and "exports the enzyme" differ in what
# they mean for the reactor; the same test runs on other modules but has no
# operational reading there.
EXPORT_MODULES_PREFIX = 'HYDROL-'

# panel_scored.py explains a gate call in a sentence; the genomes table needs
# a code, the same sentence 130 times over is noise in a spreadsheet. Matched
# on the sentence's opening words; an unrecognised note is passed through.
GATE_NOTE_CODES = [
    ('taxonomy-confirmed', 'taxonomy-confirmed'),
    ('markers present but', 'absent-on-taxonomy'),
    ('gene markers only',   'markers-only'),
    ('NOT ASSESSABLE',      'not-assessable'),
]

def gate_note_code(note):
    for prefix, code in GATE_NOTE_CODES:
        if note.startswith(prefix):
            return code
    return note

def write_tables(outdir, fmt, profile):
    """Four flat tables from the profile. One row per genome / module / route /
    risk, a 'sample' column first so files from repeated sampling concatenate.
    Unmeasured is written as NA, never 0, in every column."""
    os.makedirs(outdir, exist_ok=True)
    delim = '\t' if fmt == 'tsv' else ','
    def table(name, header, rows):
        path = os.path.join(outdir, f'{name}.{fmt}')
        with open(path, 'w', newline='') as fh:
            w = csv.writer(fh, delimiter=delim, lineterminator='\n',
                           quoting=csv.QUOTE_MINIMAL)
            w.writerow(header)
            w.writerows(rows)
        return path
    def na(v):
        return 'NA' if v is None else v
    def yn(v):
        return 'NA' if v is None else ('yes' if v else 'no')

    sample = profile['sample']
    mags = profile['mags']
    n_tot = profile['n_mags']
    gate_cols = [('acetoclastic', ROUTE_GATES['acetoclastic']),
                 ('hydrogenotrophic', ROUTE_GATES['hydrogenotrophic']),
                 ('methylotrophic', ROUTE_GATES['methylotrophic'])] + \
                [(k, c['label']) for k, c in profile['capabilities'].items()
                 if k != 'methane_oxidation']        # a lineage call, has its own column
    # panel order, taken from the first genome that has a module table
    mod_order = next(([r['module'] for r in m['modules']] for m in mags if m['modules']), [])
    export_mods = [x for x in mod_order if x.startswith(EXPORT_MODULES_PREFIX)]
    marker_order = next(([k for k in m['markers']] for m in mags if m['markers']), [])

    # ---- genomes ----
    hdr = (['sample', 'genome', 'abundance_pct', 'lineage', 'phylum', 'genus', 'mcrA', 'methane_oxidiser']
           + [k for k, _ in gate_cols] + ['gate_notes']
           + [f'{x}_pct_complete' for x in mod_order]
           + [c for x in export_mods for c in (f'{x}_families_found', f'{x}_families_exported')]
           + [f'marker_{k}' for k in marker_order])
    rows = []
    order = sorted(mags, key=lambda m: (-(m['abundance_pct'] if m['abundance_pct'] is not None else -1), m['name']))
    for m in order:
        mods = {r['module']: r for r in m['modules']}
        notes = '; '.join(f"{lab.split(': ', 1)[-1]}:{gate_note_code(g['note'])}"
                          for lab, g in m['gates'].items()
                          if g.get('note') and g.get('value') is not None)
        row = [sample, m['name'], na(m['abundance_pct']), m['lineage'], m['phylum'], m['genus'],
               yn(m['has_mcrA']), yn(m['is_methane_oxidiser'])]
        row += [yn((m['gates'].get(lab) or {}).get('value')) for _, lab in gate_cols]
        row.append(notes)
        row += [mods.get(x, {}).get('pct_complete', 'NA') or 'NA' for x in mod_order]
        for x in export_mods:
            r = mods.get(x, {})
            row += [r.get('core_found', 'NA') or 'NA', r.get('core_secreted', 'NA') or 'NA']
        row += [yn(m['markers'].get(k)) for k in marker_order]
        rows.append(row)
    table('genomes', hdr, rows)

    # ---- modules ----
    hdr = ['sample', 'module', 'branch', 'n_genomes', 'n_scored', 'n_unscoreable', 'status',
           'n_carriers_any', 'n_carriers_half', 'n_complete', 'carriers_half_abundance_pct',
           'best_pct_complete', 'best_genome', 'best_genome_genus',
           'n_exporters', 'exporters_abundance_pct']
    genus = {m['name']: m['genus'] for m in mags}
    rows = []
    for x in mod_order or sorted(profile['modules']):
        st = profile['modules'].get(x)
        if not st:
            continue
        scored = st['n_scored'] > 0
        exp = x.startswith(EXPORT_MODULES_PREFIX) and scored
        rows.append([sample, x, st['branch'], n_tot, st['n_scored'], st['n_unscoreable'],
                     '|'.join(st['statuses']),
                     na(st['n_any'] if scored else None),
                     na(st['n_half'] if scored else None),
                     na(st['n_full'] if scored else None),
                     na(st['abundance_half_pct'] if scored else None),
                     na(st['best_pct']), st['best_mag'], genus.get(st['best_mag'], ''),
                     na(st['n_exporters'] if exp else None),
                     na(st['exporters_abundance_pct'] if exp else None)])
    table('modules', hdr, rows)

    # ---- routes and capabilities ----
    hdr = ['sample', 'kind', 'key', 'label', 'n_present', 'n_genomes', 'n_unassessable',
           'abundance_pct', 'genomes']
    rows = []
    for k, r in profile['routes'].items():
        who = '; '.join(p['mag'] + (f" ({p['genus']})" if p['genus'] else '') for p in r['present_in'])
        rows.append([sample, 'methanogenesis route', k, ROUTE_GATES[k], r['n_present'], n_tot,
                     len(r['unassessable_in']), na(r['abundance_pct']), who])
    for k, c in profile['capabilities'].items():
        who = '; '.join(n + (f" ({genus[n]})" if genus.get(n) else '') for n in c['carriers'])
        rows.append([sample, 'capability', k, c['label'], len(c['carriers']), n_tot,
                     len(c['unassessable_in']), na(c['abundance_pct']), who])
    table('routes', hdr, rows)

    # ---- risks ----
    table('risks', ['sample', 'level', 'title', 'detail'],
          [[sample, r['level'], r['title'], r['detail']] for r in profile['risks']])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', required=True, help='directory of *.summary.txt / *.modules.tsv')
    ap.add_argument('--sample', default='sample', help='sample / digester name')
    ap.add_argument('--description', default='', help='what this sample actually is')
    ap.add_argument('--abundance', metavar='MANIFEST.TSV',
                    help='per-MAG relative abundance (%%); adds abundance-weighted shares')
    ap.add_argument('--acetate-lineages', default=ACETATE_LINEAGES, metavar='TSV',
                    help='lineage policy (default: assets/acetate_lineages.tsv next to bin/)')
    ap.add_argument('--out-json', required=True)
    ap.add_argument('--out-txt')
    ap.add_argument('--out-tables', metavar='DIR',
                    help='write genomes/modules/routes/risks tables into DIR')
    ap.add_argument('--table-format', choices=('tsv', 'csv'), default='tsv',
                    help='delimiter for --out-tables (default tsv; csv is RFC-4180 '
                         'comma-separated with quoting)')
    a = ap.parse_args()

    if not os.path.isfile(a.acetate_lineages):
        sys.exit(f'!! lineage policy not found: {a.acetate_lineages}')
    policy = load_acetate_lineages(a.acetate_lineages)
    abundance, ab_meta = ({}, {})
    if a.abundance:
        abundance, ab_meta = load_abundance(a.abundance)

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
            'has_mcrA': bool(markers.get('mcrA')),
            # mcr carriers that run the pathway in reverse (ANME, alkane oxidisers) are not
            # methanogens; the lineage policy names them
            'is_methane_oxidiser': bool(markers.get('mcrA')) and bool(lineage)
                                   and lineage_in_role(lineage, 'methane_oxidiser', policy),
            'gates': gates, 'markers': markers, 'genus_hint': hint,
            'modules': modules,
            'abundance_pct': abundance.get(n),
        })

    for m in mags:
        m['is_methanogen'] = m['has_mcrA'] and not m['is_methane_oxidiser']
    n_tot = len(mags)
    methanogens = [m for m in mags if m['is_methanogen']]
    methane_oxidisers = [m for m in mags if m['is_methane_oxidiser']]

    # ---- abundance bookkeeping ----
    # Shares are of reads mapped to the whole catalogue, as the manifest says, so
    # the scored MAGs never sum to 100. State what fraction they do cover.
    def share(names):
        return round(sum(abundance.get(n, 0.0) for n in names), 3) if abundance else None
    scored_share = share([m['name'] for m in mags])
    n_with_ab = sum(1 for m in mags if m['abundance_pct'] is not None)
    if abundance and n_with_ab < n_tot:
        print(f"!! {n_tot - n_with_ab} scored MAG(s) have no abundance in {a.abundance}; "
              f"they count as 0 in shares", file=sys.stderr)

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
                       'unassessable_in': unassessable,
                       'abundance_pct': share([p['mag'] for p in present])}

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
                             'unassessable_in': gate_blind(label),
                             'abundance_pct': share(gate_carriers(label))}
    # not a gate: a lineage call on mcr carriers (assets/acetate_lineages.tsv, role methane_oxidiser)
    capabilities['methane_oxidation'] = {
        'label': 'lineage: anaerobic methane / alkane oxidiser',
        'carriers': [m['name'] for m in methane_oxidisers],
        'unassessable_in': [],
        'abundance_pct': share([m['name'] for m in methane_oxidisers])}

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
                'statuses': set(),
                # carriers at three completeness levels; a community function
                # needs one organism, so these are counts and shares, not means
                'carriers_any': [], 'carriers_half': [], 'carriers_full': [],
                # genomes whose found core markers include an exported protein
                # (panel_scored.py --secretion); None until any genome was tested
                'exporters': [], 'n_secretion_tested': 0})
            st['statuses'].add(r.get('status', ''))
            sec = r.get('core_secreted', 'NA')
            if sec not in ('NA', ''):
                st['n_secretion_tested'] += 1
                try:
                    if float(sec) > 0:
                        st['exporters'].append(m['name'])
                except ValueError:
                    pass
            pct = r.get('pct_complete', 'NA')
            if pct == 'NA' or pct == '':
                st['n_unscoreable'] += 1
                continue
            st['n_scored'] += 1
            try: v = float(pct)
            except ValueError: continue
            if v > 0:    st['carriers_any'].append(m['name'])
            if v >= 50:  st['carriers_half'].append(m['name'])
            if v >= 100: st['carriers_full'].append(m['name'])
            if st['best_pct'] is None or v > st['best_pct']:
                st['best_pct'] = v; st['best_mag'] = m['name']
    for st in mod_stats.values():
        st['statuses'] = sorted(x for x in st['statuses'] if x)
        st['n_any'] = len(st['carriers_any'])
        st['n_half'] = len(st['carriers_half'])
        st['n_full'] = len(st['carriers_full'])
        st['abundance_half_pct'] = share(st['carriers_half'])
        tested = st.pop('n_secretion_tested') > 0
        st['n_exporters'] = len(st['exporters']) if tested else None
        st['exporters_abundance_pct'] = share(st['exporters']) if tested else None

    # ---- risk reading ----
    # Deliberately conservative: these are capability statements about DNA, not
    # measurements of activity, and the wording should not imply otherwise.
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
    if methane_oxidisers:
        risks.append({
            'level': 'note',
            'title': 'Anaerobic methane oxidisers present',
            'detail': (f'{len(methane_oxidisers)} MAG(s) carry mcr but belong to anaerobic methane or '
                       'alkane oxidiser lineages, which run the pathway in reverse. They are not counted '
                       'as methanogens. Their activity needs an electron acceptor: sulfate with partner '
                       'bacteria, or nitrate, iron or manganese (Methanoperedens).')})
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
        'abundance': ({'source': a.abundance,
                       'scored_mags_share_pct': scored_share,
                       'methanogens_share_pct': share([m['name'] for m in methanogens]),
                       'manifest': ab_meta} if abundance else None),
        'mags': mags,
        'routes': routes,
        'capabilities': capabilities,
        'modules': mod_stats,
        'risks': risks,
    }
    with open(a.out_json, 'w') as fh:
        json.dump(profile, fh, indent=2, default=list)
    if a.out_tables:
        write_tables(a.out_tables, a.table_format, profile)

    lines = []
    w = lines.append
    w(f"# {a.sample} — anaerobic digestion community profile")
    if a.description:
        w(f"  {a.description}")
    w(f"  {n_tot} MAG(s); {len(methanogens)} methanogen(s) (mcrA)"
      + (f"; {len(methane_oxidisers)} anaerobic methane/alkane oxidiser(s) (mcrA, reverse pathway)"
         if methane_oxidisers else ""))
    if abundance:
        w(f"  scored MAGs cover {scored_share:.1f} % of mapped reads; methanogens "
          f"{share([m['name'] for m in methanogens]):.1f} %"
          + (f"; {ab_meta['mags_not_deposited']} MAG(s) above the cutoff "
             f"({ab_meta.get('mags_not_deposited_pct', '?')} %) have no genome and were never scored"
             if 'mags_not_deposited' in ab_meta else ''))
    w("")
    w("methanogenesis routes:")
    def pct(v):
        return f"  {v:5.1f} % of reads" if v is not None else ''
    for key in ROUTE_GATES:
        r = routes[key]
        who = ', '.join(f"{p['mag']}" + (f" ({p['genus']})" if p['genus'] else '')
                        for p in r['present_in']) or '—'
        w(f"  {key:<18} {r['n_present']}/{n_tot}{pct(r['abundance_pct'])}  {who}")
        if r['unassessable_in']:
            w(f"  {'':<18} not assessable in {len(r['unassessable_in'])} MAG(s)")
    w("")
    w("other capabilities:")
    for key, c in capabilities.items():
        extra = (f"  (not assessable in {len(c['unassessable_in'])})"
                 if c['unassessable_in'] else '')
        w(f"  {key:<20} {len(c['carriers'])}/{n_tot}{pct(c['abundance_pct'])}{extra}")
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
