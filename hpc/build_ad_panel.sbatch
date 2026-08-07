#!/usr/bin/env bash
#SBATCH --job-name=build-ad-panel
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/build-ad-panel-%j.out
#SBATCH --error=logs/build-ad-panel-%j.err
# First LIVE build of the KEGG-less AD panel HMM DB against current NCBIfam.
#
# The point of this run is verification, not just artifacts. build_ad_hmm_db.sh
# was written as a scaffold that GUESSES at two things it could not check
# offline: the hmm_PGAP.tsv column names, and the .tgz archive layout. Sections
# A-D below report what was actually found so those guesses can be confirmed or
# corrected, and so the EC-derived matches can be reviewed before anything built
# here is trusted in a report.
#
# Prereqs (login node):
#   bash kegg-less/hpc/prefetch_ncbifam.sh
#   mkdir -p logs
# Submit from the worktree root:
#   sbatch kegg-less/hpc/build_ad_panel.sbatch
set -euo pipefail
module purge
module load hmmer/3.3.2

# SLURM copies the batch script into its spool dir, so BASH_SOURCE does NOT
# point into the repo under sbatch -- use the submit dir there, and fall back to
# BASH_SOURCE when this is run directly with bash.
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  KEGGLESS="$SLURM_SUBMIT_DIR/kegg-less"
else
  KEGGLESS="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
[ -f "$KEGGLESS/build_ad_hmm_db.sh" ] || {
  echo "!! cannot find build_ad_hmm_db.sh under $KEGGLESS" >&2
  echo "   submit from the worktree root: sbatch kegg-less/hpc/build_ad_panel.sbatch" >&2
  exit 1
}
DB="${DB:-/gpfs/space/projects/preterm/databases}"
export PANEL="$KEGGLESS/AD_methanogenesis_panel.tsv"
export WORK="${WORK:-$DB/ncbifam/_work}"
export OUT="${OUT:-$DB/ad_panel}"

echo "==> host=$(hostname) job=${SLURM_JOB_ID:-none} $(date)"
echo "==> hmmer: $(hmmpress -h 2>&1 | grep -m1 -i hmmer || echo unknown)"
echo "==> panel=$PANEL"
echo "==> work=$WORK  out=$OUT"
echo

bash "$KEGGLESS/build_ad_hmm_db.sh"

echo
echo "================ A. NCBIfam metadata columns ================"
echo "(build script guessed: ncbi_accession / gene_symbol / ec_numbers / label)"
head -1 "$WORK/hmm_PGAP.tsv" | tr '\t' '\n' | nl
echo
echo "================ B. archive layout ================"
# Materialise the listing before slicing it: `tar -tzf ... | head` kills tar
# with SIGPIPE, which under `set -o pipefail` fails the whole job -- and it
# saves decompressing 424M twice.
TARLIST="$WORK/archive_listing.txt"
[ -s "$TARLIST" ] || tar -tzf "$WORK/hmm_PGAP.HMM.tgz" > "$TARLIST"
head -5 "$TARLIST"
echo "  ... total entries: $(wc -l < "$TARLIST")"
echo "  models available:  $(wc -l < "$WORK/available_acc.txt")"
echo
echo "================ C. panel coverage ================"
MAP="$OUT/ad_panel_map.tsv" python3 - <<'PY'
import csv, os, re
# gene_tokens mirrors build_ad_hmm_db.sh -- keep in sync if that splitter changes.
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

covered = {(r.get('gene') or '').lower()
           for r in csv.DictReader(open(os.environ['MAP']), delimiter='\t')}

by_module, missing = {}, []
for r in csv.DictReader(open(os.environ['PANEL']), delimiter='\t'):
    mod = r['module']; toks = gene_tokens(r['gene'])
    ok = any(t.lower() in covered for t in toks)
    n_tot, n_ok = by_module.get(mod, (0, 0))
    by_module[mod] = (n_tot + 1, n_ok + (1 if ok else 0))
    if not ok:
        missing.append((mod, r['gene'], (r.get('tier') or '').strip()))

print(f"{'module':<28} {'genes':>6} {'with model':>11}")
for mod, (tot, ok) in sorted(by_module.items()):
    flag = '' if ok == tot else '   <-- gaps'
    print(f"{mod:<28} {tot:>6} {ok:>11}{flag}")
tot = sum(t for t, _ in by_module.values()); okc = sum(o for _, o in by_module.values())
print(f"\nTOTAL {okc}/{tot} panel genes resolved to >=1 NCBIfam model")
if missing:
    print(f"\nUNRESOLVED ({len(missing)}) -- need Pfam fallback, dbCAN, or panel fix:")
    for mod, gene, tier in missing:
        print(f"  [{mod}] {gene}" + (f"  (tier {tier})" if tier else ""))
PY
echo
echo "================ D. EC-derived matches (REVIEW THESE) ================"
ECLIST="$WORK/ec_derived.txt"
awk -F'\t' 'NR>1 && $5 ~ /^ec:/ {print "  " $1 "\t" $3 "\t" $5 "\t" $2}' "$OUT/ad_panel_map.tsv" \
  | LC_ALL=C sort -k2,2 > "$ECLIST"
head -60 "$ECLIST"
echo "  ... EC-derived total: $(wc -l < "$ECLIST")"
echo
echo "==> artifacts:"
ls -lh "$OUT"
echo "==> done $(date)"
