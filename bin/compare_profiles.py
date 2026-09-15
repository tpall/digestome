#!/usr/bin/env python3
"""Compare two or more community profiles of the same reactor over time.

aggregate_community.py answers "what can this community do" for one sample. A
monitoring report has to answer "what moved since last time", which is a
different table: rows are routes, capabilities, modules and genomes; columns
are the samples in time order; the last column is the direction.

Inputs are profile.json files in sampling order (oldest first). Optional
labels (--labels) replace the sample names, e.g. dates.

Outputs (--out DIR):
  routes_over_time.tsv     routes + capabilities: n genomes and share per sample, change
  modules_over_time.tsv    per module: carriers at >= 50 % completeness and their share
  genomes_over_time.tsv    every genome seen in any sample: share per sample, change,
                           status (new / lost / up / down / stable)
  changes.json             the same, machine-readable, plus the summary lines
  summary.txt              what moved, in the order a report would say it

Change calls are deliberately coarse. With two or three points there is no
trend to test, only a difference; a genome is "up" or "down" when its share
changed by at least --min-delta percentage points AND by at least --min-fold,
so 0.3 % -> 0.4 % is never a change and 0.1 % -> 3 % always is. Anything
below the abundance floor in every sample is dropped from the genome table.

Stdlib only.
"""
import argparse, csv, json, os, sys

ROUTE_ORDER = ['acetoclastic', 'hydrogenotrophic', 'methylotrophic']
CAP_ORDER = ['wood_ljungdahl', 'butyrate_oxidation', 'ethanolamine']

def load(path):
    with open(path) as fh:
        return json.load(fh)

def share(v):
    return 0.0 if v is None else float(v)

def direction(first, last, min_delta, min_fold):
    """Coarse change call between the first and the last sample."""
    if first is None and last is None:
        return 'NA'
    a, b = share(first), share(last)
    if a == 0 and b == 0:
        return 'absent'
    if a == 0:
        return 'new'
    if b == 0:
        return 'lost'
    d = b - a
    fold = b / a if a else float('inf')
    if abs(d) >= min_delta and (fold >= min_fold or fold <= 1 / min_fold):
        return 'up' if d > 0 else 'down'
    return 'stable'

def fmt(v, nd=1):
    return 'NA' if v is None else f'{v:.{nd}f}'

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('profiles', nargs='+', help='profile.json files, oldest first')
    ap.add_argument('--labels', nargs='*', help='column labels in the same order (e.g. dates)')
    ap.add_argument('--reactor', default='reactor', help='name for the report')
    ap.add_argument('--out', required=True, help='output directory')
    ap.add_argument('--min-delta', type=float, default=1.0,
                    help='percentage points a share must move to count as a change (default 1.0)')
    ap.add_argument('--min-fold', type=float, default=1.5,
                    help='fold change a share must show to count as a change (default 1.5)')
    ap.add_argument('--floor', type=float, default=0.5,
                    help='genomes below this share (%%) in every sample are left out of the genome table')
    a = ap.parse_args()
    if len(a.profiles) < 2:
        sys.exit('!! need at least two profiles')
    profs = [load(p) for p in a.profiles]
    labels = a.labels or [p['sample'] for p in profs]
    if len(labels) != len(profs):
        sys.exit('!! --labels must match the number of profiles')
    os.makedirs(a.out, exist_ok=True)
    def tsv(name, header, rows):
        with open(os.path.join(a.out, name), 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t', lineterminator='\n')
            w.writerow(header); w.writerows(rows)

    out = {'reactor': a.reactor, 'samples': labels,
           'n_mags': [p['n_mags'] for p in profs],
           'n_methanogens': [p['n_methanogens'] for p in profs],
           'methanogen_share': [(p.get('abundance') or {}).get('methanogens_share_pct') for p in profs],
           'scored_share': [(p.get('abundance') or {}).get('scored_mags_share_pct') for p in profs]}

    # ---- routes and capabilities ----
    rows, jrows = [], []
    for kind, keys, getter in (
            ('route', ROUTE_ORDER, lambda p, k: (p['routes'][k]['n_present'], p['routes'][k]['abundance_pct'],
                                                  len(p['routes'][k]['unassessable_in']))),
            ('capability', CAP_ORDER, lambda p, k: (len(p['capabilities'][k]['carriers']),
                                                    p['capabilities'][k]['abundance_pct'],
                                                    len(p['capabilities'][k]['unassessable_in'])))):
        for k in keys:
            vals = [getter(p, k) for p in profs]
            n_total = [p['n_mags'] for p in profs]
            if all(v[2] == nt for v, nt in zip(vals, n_total)):
                d = 'not assessable'
            else:
                d = direction(vals[0][1], vals[-1][1], a.min_delta, a.min_fold)
            row = [kind, k] + [f'{v[0]}' for v in vals] + [fmt(v[1]) for v in vals] + [d]
            rows.append(row)
            jrows.append({'kind': kind, 'key': k, 'n': [v[0] for v in vals],
                          'share_pct': [v[1] for v in vals], 'direction': d})
    tsv('routes_over_time.tsv',
        ['kind', 'key'] + [f'n_{l}' for l in labels] + [f'share_{l}' for l in labels] + ['direction'], rows)
    out['routes'] = jrows

    # ---- modules ----
    rows, jrows = [], []
    mods = list(profs[0]['modules'].keys())
    for m in mods:
        sts = [p['modules'].get(m) for p in profs]
        if any(st is None for st in sts):
            continue
        if all(st['n_scored'] == 0 for st in sts):
            d = 'not assessable'; nh = ['NA'] * len(sts); sh = ['NA'] * len(sts)
        else:
            nh = [st['n_half'] for st in sts]
            sh = [st['abundance_half_pct'] for st in sts]
            d = direction(sh[0], sh[-1], a.min_delta, a.min_fold)
        rows.append([m, sts[0]['branch']] + [str(x) for x in nh] + [fmt(x) if x != 'NA' else 'NA' for x in sh] + [d])
        jrows.append({'module': m, 'branch': sts[0]['branch'], 'n_half': nh, 'share_half_pct': sh, 'direction': d})
    tsv('modules_over_time.tsv',
        ['module', 'branch'] + [f'carriers_half_{l}' for l in labels] + [f'share_half_{l}' for l in labels] + ['direction'], rows)
    out['modules'] = jrows

    # ---- genomes ----
    seen = {}
    for i, p in enumerate(profs):
        for m in p['mags']:
            g = seen.setdefault(m['name'], {'lineage': m['lineage'], 'genus': m['genus'],
                                            'is_methanogen': m['is_methanogen'],
                                            'routes': set(), 'share': [None] * len(profs)})
            g['share'][i] = m['abundance_pct']
            g['is_methanogen'] = g['is_methanogen'] or m['is_methanogen']
            for k, lab in (('acetoclastic', 'methanogenesis: acetoclastic'),
                           ('hydrogenotrophic', 'methanogenesis: hydrogenotrophic'),
                           ('methylotrophic', 'methanogenesis: methylotrophic')):
                if (m['gates'].get(lab) or {}).get('value'):
                    g['routes'].add(k)
    rows, jrows = [], []
    for name, g in seen.items():
        sh = [share(x) for x in g['share']]
        if max(sh) < a.floor:
            continue
        d = direction(sh[0], sh[-1], a.min_delta, a.min_fold)
        rows.append([name, g['genus'] or g['lineage'], 'yes' if g['is_methanogen'] else 'no',
                     '|'.join(sorted(g['routes'])) or '-'] + [fmt(x, 2) for x in sh] + [fmt(sh[-1] - sh[0], 2), d])
        jrows.append({'genome': name, 'genus': g['genus'], 'lineage': g['lineage'],
                      'is_methanogen': g['is_methanogen'], 'routes': sorted(g['routes']),
                      'share_pct': sh, 'delta_pct': sh[-1] - sh[0], 'direction': d})
    order = {'new': 0, 'up': 1, 'down': 2, 'lost': 3, 'stable': 4, 'absent': 5}
    rows.sort(key=lambda r: (order.get(r[-1], 9), -abs(float(r[-2]))))
    jrows.sort(key=lambda r: (order.get(r['direction'], 9), -abs(r['delta_pct'])))
    tsv('genomes_over_time.tsv',
        ['genome', 'taxon', 'methanogen', 'routes'] + [f'share_{l}' for l in labels] + ['delta_pct', 'direction'], rows)
    out['genomes'] = jrows

    # ---- summary ----
    lines = [f"# {a.reactor} — {len(profs)} samples: {', '.join(labels)}",
             f"  genomes scored: {' -> '.join(str(n) for n in out['n_mags'])}; "
             f"methanogens: {' -> '.join(str(n) for n in out['n_methanogens'])} "
             f"({' -> '.join(fmt(x) for x in out['methanogen_share'])} % of reads)", '',
             'routes:']
    for r in out['routes']:
        if r['kind'] != 'route': continue
        lines.append(f"  {r['key']:<18} {' -> '.join(f'{n}/{fmt(s)} %' for n, s in zip(r['n'], r['share_pct']))}   {r['direction']}")
    lines.append(''); lines.append('capabilities:')
    for r in out['routes']:
        if r['kind'] != 'capability': continue
        lines.append(f"  {r['key']:<18} {' -> '.join(f'{n}/{fmt(s)} %' for n, s in zip(r['n'], r['share_pct']))}   {r['direction']}")
    moved = [g for g in out['genomes'] if g['direction'] in ('new', 'up', 'down', 'lost')]
    lines.append(''); lines.append(f"genomes that moved (>= {a.min_delta} pp and >= {a.min_fold}x; {len(moved)} of {len(out['genomes'])} above {a.floor} %):")
    for g in moved[:15]:
        tag = ' [methanogen: ' + ','.join(g['routes']) + ']' if g['is_methanogen'] else ''
        lines.append(f"  {g['direction']:<6} {g['genus'] or g['lineage'][:40]:<32} {' -> '.join(fmt(x, 2) for x in g['share_pct'])} %{tag}")
    stable_meth = [g for g in out['genomes'] if g['is_methanogen'] and g['direction'] == 'stable']
    if stable_meth:
        lines.append(f"  stable methanogens: {', '.join(g['genus'] or g['genome'] for g in stable_meth)}")
    text = '\n'.join(lines) + '\n'
    open(os.path.join(a.out, 'summary.txt'), 'w').write(text)
    out['summary'] = lines
    with open(os.path.join(a.out, 'changes.json'), 'w') as fh:
        json.dump(out, fh, indent=2)
    print(text)

if __name__ == '__main__':
    main()
