#!/usr/bin/env python3
"""Check panel_scored.py output for one genome against the expected biology.

Parses the summary that panel_scored.py prints and compares it to the rows for
this genome in expectations.tsv. Exits non-zero on any mismatch so the sbatch
job actually fails rather than printing a wall of text nobody reads.

  check_expectations.py --label Ecoli --summary out.txt --expect expectations.tsv

Only well-established biology is asserted, so a FAIL means either the panel is
wrong or the expectation is -- both worth knowing. Anything not asserted is
ignored here and stays visible in the summary itself.

Stdlib only.
"""
import argparse, csv, re, sys

def parse_summary(path):
    """-> ({gate label: True/False/None}, {marker: bool}, genus_hint)"""
    gates, markers, hint = {}, {}, ''
    section = None
    for line in open(path):
        s = line.rstrip('\n')
        if s.startswith('acetoclastic genus hint:'):
            hint = s.split(':', 1)[1].strip()
        if s.startswith('branch gates:'):
            section = 'gates';  continue
        if s.startswith('diagnostic markers:'):
            section = 'markers'; continue
        if s.startswith(('CAZyme', 'not scored', 'module table')):
            section = None
        if section == 'gates':
            m = re.match(r'\s*\[([x? ])\]\s+(.*?)(?:\s+—\s+NOT ASSESSABLE.*)?$', s)
            if m:
                flag, label = m.group(1), m.group(2).strip()
                gates[label] = None if flag == '?' else (flag == 'x')
        elif section == 'markers':
            m = re.match(r'\s*(\S+)\s+([+-])\s*$', s)
            if m:
                markers[m.group(1)] = (m.group(2) == '+')
    return gates, markers, hint

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--label', required=True)
    ap.add_argument('--summary', required=True)
    ap.add_argument('--expect', required=True)
    a = ap.parse_args()

    gates, markers, hint = parse_summary(a.summary)

    rows = []
    with open(a.expect) as fh:
        for r in csv.reader((l for l in fh if not l.startswith('#') and l.strip()),
                            delimiter='\t'):
            if len(r) >= 4 and r[0] == a.label:
                rows.append(r[:4])

    if not rows:
        print(f"  {a.label}: no expectations declared — reporting only")
        return 0

    fails = 0
    for _lbl, kind, name, want in rows:
        if kind == 'gate':
            got = gates.get(name, '<missing>')
            ok = (got is None) if want == 'na' else (got is (want == 'yes'))
            shown = 'NOT-ASSESSABLE' if got is None else got
        elif kind == 'marker':
            got = markers.get(name, '<missing>')
            ok = got is (want == 'yes')
            shown = got
        elif kind == 'hint':
            got = hint
            ok = (name.lower() in hint.lower()) == (want == 'yes')
            shown = repr(hint)
        else:
            print(f"  ?? unknown expectation kind '{kind}'", file=sys.stderr)
            fails += 1
            continue
        if not ok:
            fails += 1
        print(f"  [{'PASS' if ok else 'FAIL'}] {a.label:<16} {kind:<6} "
              f"{name:<38} want={want:<3} got={shown}")

    if fails:
        print(f"  {a.label}: {fails}/{len(rows)} FAILED", file=sys.stderr)
    return 1 if fails else 0

if __name__ == '__main__':
    sys.exit(main())
