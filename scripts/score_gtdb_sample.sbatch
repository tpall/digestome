#!/usr/bin/env bash
#SBATCH --job-name=gtdb-panel
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=logs/gtdb-panel-%j.out
#SBATCH --error=logs/gtdb-panel-%j.err
#
# Score the phylum-stratified GTDB representative sample.
#
#   export DB=/path/to/databases
#   mkdir -p logs && sbatch -A <account> -p <partition> scripts/score_gtdb_sample.sbatch
#
# Runs directly under bash too; the SLURM variables all have fallbacks.
#
# Why this set: the digester MAG benchmark establishes specificity on genomes
# assembled and binned by one group with one pipeline. These are isolate-derived
# reference genomes, not binned at all, so agreement across both cannot be an
# artefact of anyone's binning.
#
# It doubles as the specificity test for PF05369 (MtmB). That family is NOT in
# the panel: monomethylamine methyltransferase carries an in-frame amber codon,
# and whether it is detectable depends on the gene caller upstream. It is scored
# here separately so the decision to include it rests on a large negative set
# rather than on the seven reference genomes.
set -uo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
[ -f "$ROOT/bin/panel_scored.py" ] || {
  echo "!! cannot find bin/panel_scored.py under $ROOT -- submit from the repo root" >&2; exit 1; }

: "${DB:?set DB to a database root}"
WORK="${WORK:-$DB/ad_panel_test/gtdb}"
PROT="${PROTEOMES:-$WORK/proteomes_prodigal}"
OUT="${OUT:-$WORK/scored}"          # override to score the same input two ways
CPUS="${SLURM_CPUS_PER_TASK:-4}"
HMM="$DB/ad_panel/ad_panel.hmm"
MAP="$DB/ad_panel/ad_panel_map.tsv"
PANEL="$ROOT/assets/AD_methanogenesis_panel.tsv"
TAX="$WORK/gtdbtk_sample.tsv"

command -v hmmsearch >/dev/null 2>&1 || module load "${HMMER_MODULE:-hmmer}" 2>/dev/null || true
command -v hmmsearch >/dev/null 2>&1 || { echo "!! hmmsearch not available" >&2; exit 1; }
for f in "$HMM" "$MAP" "$PANEL" "$TAX"; do
  [ -s "$f" ] || { echo "!! missing: $f" >&2; exit 1; }
done
n=$(ls "$PROT"/*.faa 2>/dev/null | wc -l)
[ "$n" -gt 0 ] || { echo "!! no proteomes in $PROT -- run fetch_gtdb_sample.sh first" >&2; exit 1; }
mkdir -p "$OUT"
echo "==> $n proteomes, $CPUS cpus, started $(date)"

# ---- score each genome -----------------------------------------------------
# CONTIG_LEVEL=1 requires every term of a gate on one contig. Meant for
# whole-assembly input; on binned genomes it only costs recall.
export CONTIG_FLAG="${CONTIG_LEVEL:+--contig-level}"

score_one() {
  local faa="$1" out="$2" hmm="$3" map="$4" panel="$5" tax="$6" sec="$7"
  local acc; acc=$(basename "$faa" .faa)
  [ -s "$out/$acc.modules.tsv" ] && return 0
  local tbl="$out/$acc.tblout"
  hmmsearch --cpu 1 --cut_nc --tblout "$tbl" "$hmm" "$faa" > /dev/null 2>&1 || return 0
  python3 "$out/../panel_scored.py" --tblout "$tbl" --map "$map" --panel "$panel" \
      --secretion "$sec" --gtdbtk "$tax" --name "$acc" ${CONTIG_FLAG:-} \
      --acetate-lineages "$ROOT/assets/acetate_lineages.tsv" \
      --out "$out/$acc.modules.tsv" > "$out/$acc.summary.txt" 2>> "$out/scoring.err"
  rm -f "$tbl"
}
export -f score_one
export ROOT
cp "$ROOT/bin/panel_scored.py" "$WORK/panel_scored.py"

ls "$PROT"/*.faa | xargs -P "$CPUS" -I{} bash -c \
  'score_one "$@"' _ {} "$OUT" "$HMM" "$MAP" "$PANEL" "$TAX" "$ROOT/assets/secretion_modules.tsv"
echo "==> scored $(ls "$OUT"/*.modules.tsv 2>/dev/null | wc -l) genomes"
n_empty=$(find "$OUT" -name '*.summary.txt' -empty | wc -l)
if [ "$n_empty" -gt 0 ]; then
  echo "!! $n_empty empty summaries -- the scorer failed; last errors:" >&2
  tail -20 "$OUT/scoring.err" >&2 2>/dev/null
  exit 1
fi

# ---- the specificity question ----------------------------------------------
python3 - "$WORK" <<'PY'
import sys, pathlib, csv, collections
work = pathlib.Path(sys.argv[1])
stratum, order = {}, {}
with open(work/'gtdb_sample.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        stratum[r['accession']] = r['stratum']; order[r['accession']] = r['order']

routes = collections.defaultdict(set)
for f in (work/'scored').glob('*.summary.txt'):
    acc = f.name[:-len('.summary.txt')]
    for line in f.read_text().splitlines():
        s = line.strip()
        if s.startswith('methanogenesis:') and ' YES' in s.upper():
            routes[acc].add(s.split(':')[1].split()[0])

scored = {f.name[:-len('.modules.tsv')] for f in (work/'scored').glob('*.modules.tsv')}
tally = collections.Counter()
for acc in scored:
    st = stratum.get(acc, 'unknown')
    tally[(st, bool(routes.get(acc)))] += 1

print('\n================ GTDB SPECIFICITY ================')
print(f'  genomes scored: {len(scored)}')
for st in ('background', 'methanogen'):
    withr = tally[(st, True)]; without = tally[(st, False)]
    print(f'  {st:11s}  route assigned: {withr:5d}   no route: {without:5d}')
fp = sorted(a for a in scored if stratum.get(a) == 'background' and routes.get(a))
print(f'\n  FALSE POSITIVES (background genomes given a route): {len(fp)}')
for a in fp[:25]:
    print(f'    {a}  {order.get(a,"?")}  -> {sorted(routes[a])}')
fn = sorted(a for a in scored if stratum.get(a) == 'methanogen' and not routes.get(a))
print(f'\n  METHANOGENS WITH NO ROUTE: {len(fn)}')
byorder = collections.Counter(order.get(a,'?') for a in fn)
for o, c in byorder.most_common(12):
    print(f'    {o}: {c}')
PY

# ---- PF05369 (MtmB) specificity, separate from the panel -------------------
if [ -s "$WORK/PF05369.hmm" ]; then
  echo
  echo "==> PF05369 (MtmB) across the same set"
  : > "$WORK/pf05369_hits.txt"
  for faa in "$PROT"/*.faa; do
    acc=$(basename "$faa" .faa)
    h=$(hmmsearch --cpu 1 --cut_ga "$WORK/PF05369.hmm" "$faa" 2>/dev/null | grep -c '^>>')
    [ "${h:-0}" -gt 0 ] && echo -e "$acc\t$h" >> "$WORK/pf05369_hits.txt"
  done
  python3 - "$WORK" <<'PY'
import sys, pathlib, csv, collections
work = pathlib.Path(sys.argv[1])
stratum, order = {}, {}
with open(work/'gtdb_sample.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        stratum[r['accession']] = r['stratum']; order[r['accession']] = r['order']
hits = {}
p = work/'pf05369_hits.txt'
if p.exists():
    for line in p.read_text().splitlines():
        if line.strip(): a, n = line.split('\t'); hits[a] = int(n)
c = collections.Counter(stratum.get(a, 'unknown') for a in hits)
print(f'  genomes with a PF05369 hit: {len(hits)}')
print(f'    in methanogens: {c["methanogen"]}')
print(f'    in background:  {c["background"]}   <- false positives if > 0')
byorder = collections.Counter(order.get(a,'?') for a in hits)
for o, n in byorder.most_common(12):
    print(f'      {o}: {n}')
PY
fi
echo "==> done $(date)"
