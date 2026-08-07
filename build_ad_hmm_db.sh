#!/usr/bin/env bash
# build_ad_hmm_db.sh — assemble a license-clean HMM profile DB for the
# AD/methanogenesis panel from NCBIfam (public-domain) models, and fetch the
# EC->Rhea map (CC BY). Part of DRAM branch feature/kegg-less-ad.
#
# Outputs (under $OUT, default ./db):
#   ad_panel.hmm(.h3*)   pressed HMMER DB for hmmsearch
#   ad_panel_map.tsv     model_accession <tab> model_name <tab> gene <tab> modules
#   rhea2ec.tsv          EC -> Rhea reactions, for reaction-level reporting
#
# Requires: curl, tar, HMMER3 (hmmfetch, hmmpress), python3.
# CAZymes (HYDROL-CARB rows) come from run_dbcan separately — see README.md.
# SCAFFOLD: verify NCBIfam URLs/column names and review EC-based matches.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PANEL="${PANEL:-$HERE/AD_methanogenesis_panel.tsv}"
OUT="${OUT:-$HERE/db}"
WORK="${WORK:-$OUT/_work}"
NCBIFAM_BASE="${NCBIFAM_BASE:-https://ftp.ncbi.nlm.nih.gov/hmm/current}"
RHEA_EC_URL="${RHEA_EC_URL:-https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv}"

for t in curl tar hmmfetch hmmpress python3; do
  command -v "$t" >/dev/null || { echo "!! missing dependency: $t" >&2; exit 1; }
done
mkdir -p "$OUT" "$WORK"

echo "==> 1. Download NCBIfam HMMs + metadata"
[ -s "$WORK/hmm_PGAP.HMM.tgz" ] || curl -L --fail -C - -o "$WORK/hmm_PGAP.HMM.tgz" "$NCBIFAM_BASE/hmm_PGAP.HMM.tgz"
[ -s "$WORK/hmm_PGAP.tsv" ]     || curl -L --fail -C - -o "$WORK/hmm_PGAP.tsv"     "$NCBIFAM_BASE/hmm_PGAP.tsv"

echo "==> 2. Concatenate models into one fetchable, indexed HMM file"
if [ ! -s "$WORK/ncbifam_all.hmm" ]; then
  # Verified 2026-08: the archive is ~19k individual hmm_PGAP/NF*.HMM files.
  # Layout has changed before, so still tolerate a flat dir or lowercase .hmm.
  find_models() { find "$WORK" -name "$1" ! -name 'ncbifam_all.hmm' "${@:2}"; }

  # Skip the (slow) extraction if models are already unpacked.
  [ -n "$(find_models '*.HMM' -print -quit)" ] || [ -n "$(find_models '*.hmm' -print -quit)" ] || \
    tar -xzf "$WORK/hmm_PGAP.HMM.tgz" -C "$WORK"

  # NB: do NOT probe with `find ... | grep -q .`. grep exits on the first match,
  # find is killed by SIGPIPE partway through 19k paths, and under `pipefail`
  # the pipeline reports failure even though the models are present. Small test
  # dirs hide this -- their output fits the 64K pipe buffer. Use -print -quit.
  if   [ -n "$(find_models '*.HMM' -print -quit)" ]; then pat='*.HMM'
  elif [ -n "$(find_models '*.hmm' -print -quit)" ]; then pat='*.hmm'
  else
    echo "!! could not locate extracted HMM models under $WORK" >&2; exit 1
  fi
  echo "    concatenating $(find_models "$pat" | wc -l) models matching $pat"
  find_models "$pat" -print0 | xargs -0 cat > "$WORK/ncbifam_all.hmm"
fi
hmmfetch --index "$WORK/ncbifam_all.hmm"

echo "==> 3. Match panel genes/EC to NCBIfam models -> keylist + map"
PANEL="$PANEL" NCBIFAM_TSV="$WORK/hmm_PGAP.tsv" KEYS="$WORK/keys.txt" MAP="$OUT/ad_panel_map.tsv" \
python3 - <<'PY'
import csv, os, re
panel_path = os.environ['PANEL']; ncbi_path = os.environ['NCBIFAM_TSV']
keys_path  = os.environ['KEYS'];  map_path  = os.environ['MAP']

SKIP = {'multi','core','use','dbcan','see','pep','gh','lip','est'}
def gene_tokens(field):
    field = re.sub(r'\(.*?\)', '', field)            # drop parenthetical notes
    out = []
    for p in re.split(r'[\/\+\s]+', field):          # split on / + whitespace
        for q in p.split('-'):                       # loose split of compounds
            q = q.strip()
            if len(q) >= 3 and re.search(r'[a-z]', q) and q.lower() not in SKIP:
                out.append(q)
    return out

# --- load panel ---
panel_genes = {}   # lower -> (orig, set(modules))
panel_ec    = {}   # ec -> set(orig genes)
for r in csv.DictReader(open(panel_path), delimiter='\t'):
    mod = r['module']
    toks = gene_tokens(r['gene'])
    for g in toks:
        e = panel_genes.setdefault(g.lower(), (g, set())); e[1].add(mod)
    for ec in re.split(r'[;, ]+', r.get('ec','') or ''):
        ec = ec.strip()
        if ec and ec[0].isdigit():
            panel_ec.setdefault(ec, set()).update(toks)

# --- locate columns in NCBIfam metadata by header name ---
reader = csv.DictReader(open(ncbi_path), delimiter='\t')
fn = reader.fieldnames or []
def find(*cands):
    for c in cands:
        for f in fn:
            if f.lower().lstrip('#') == c: return f
    for c in cands:
        for f in fn:
            if c in f.lower(): return f
    return None
acc_col  = find('ncbi_accession','accession')
gene_col = find('gene_symbol','gene')
ec_col   = find('ec_numbers','ec')
name_col = find('label','source_identifier','hmm_name','name')

keys, maprows = set(), []
for r in reader:
    acc  = (r.get(acc_col)  or '').strip() if acc_col else ''
    gsym = (r.get(gene_col) or '').strip() if gene_col else ''
    name = (r.get(name_col) or '').strip() if name_col else ''
    ecs  = [e.strip() for e in re.split(r'[;, ]+', (r.get(ec_col) or '') if ec_col else '') if e.strip()]
    gene, mods, basis = None, set(), ''
    if gsym and gsym.lower() in panel_genes:                    # reliable: gene-symbol match
        gene, mods = panel_genes[gsym.lower()]
        basis = 'gene_symbol'
    else:                                                       # weaker: EC match (review!)
        for e in ecs:
            if e in panel_ec:
                gs = sorted(panel_ec[e]); gene = gs[0] if gs else gsym
                for g in panel_ec[e]:
                    if g.lower() in panel_genes: mods |= panel_genes[g.lower()][1]
                basis = 'ec:' + e
                break
    if gene and acc:
        keys.add(acc)
        maprows.append((acc, name, gene, ';'.join(sorted(mods)), basis))

with open(keys_path, 'w') as fh:
    fh.write('\n'.join(sorted(keys)) + '\n')
with open(map_path, 'w') as fh:
    fh.write('#model_accession\tmodel_name\tgene\tmodules\tmatch_basis\n')
    for row in sorted(set(maprows)):
        fh.write('\t'.join(row) + '\n')
n_sym = sum(1 for r in maprows if r[4] == 'gene_symbol')
print(f"   resolved columns: acc={acc_col} gene={gene_col} ec={ec_col} name={name_col}")
print(f"   matched {len(keys)} models -> {len(maprows)} map rows "
      f"({n_sym} by gene symbol, {len(maprows) - n_sym} by EC -- review the EC ones)")
PY

echo "==> 4. Fetch matched models -> pressed DB"
# hmm_PGAP.tsv indexes more models than the .tgz ships -- some accessions in the
# metadata (seen 2026-08: 21 of 337, all .7 versions) have no model in this
# release, and hmmfetch -f aborts on the first one it cannot find. Intersect the
# key list with the accessions actually present, and name the panel genes that
# lost a model instead of dropping them silently.
grep '^ACC' "$WORK/ncbifam_all.hmm" | awk '{print $2}' | LC_ALL=C sort -u > "$WORK/available_acc.txt"
LC_ALL=C sort -u "$WORK/keys.txt" > "$WORK/keys_sorted.txt"
LC_ALL=C comm -12 "$WORK/keys_sorted.txt" "$WORK/available_acc.txt" > "$WORK/keys_present.txt"
LC_ALL=C comm -23 "$WORK/keys_sorted.txt" "$WORK/available_acc.txt" > "$WORK/keys_missing.txt"

n_want=$(wc -l < "$WORK/keys_sorted.txt"); n_have=$(wc -l < "$WORK/keys_present.txt")
n_missing=$(wc -l < "$WORK/keys_missing.txt")
echo "    $n_have/$n_want matched accessions present in this NCBIfam release"
if [ "$n_missing" -gt 0 ]; then
  echo "    !! $n_missing absent from the archive (accession / gene / match_basis):"
  awk -F'\t' 'NR==FNR{miss[$1];next} FNR>1 && ($1 in miss){print "       " $1 "\t" $3 "\t" $5}' \
      "$WORK/keys_missing.txt" "$OUT/ad_panel_map.tsv" | LC_ALL=C sort -k2,2
fi

hmmfetch -f "$WORK/ncbifam_all.hmm" "$WORK/keys_present.txt" > "$OUT/ad_panel.hmm"
hmmpress -f "$OUT/ad_panel.hmm"

# Keep the map honest: it must describe only what is in the pressed DB.
{ head -1 "$OUT/ad_panel_map.tsv"
  awk -F'\t' 'NR==FNR{ok[$1];next} FNR>1 && ($1 in ok)' "$WORK/keys_present.txt" "$OUT/ad_panel_map.tsv"
} > "$OUT/ad_panel_map.tsv.tmp"
mv "$OUT/ad_panel_map.tsv.tmp" "$OUT/ad_panel_map.tsv"

echo "==> 5. EC -> Rhea map (CC BY)"
[ -s "$OUT/rhea2ec.tsv" ] || curl -L --fail -o "$OUT/rhea2ec.tsv" "$RHEA_EC_URL"

echo
echo "==> Done."
echo "    DB:  $OUT/ad_panel.hmm"
echo "    map: $OUT/ad_panel_map.tsv"
echo "    Run: hmmsearch --cut_nc --tblout hits.tblout $OUT/ad_panel.hmm proteins.faa"
