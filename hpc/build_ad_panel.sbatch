#!/usr/bin/env bash
#SBATCH --job-name=build-ad-panel
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/build-ad-panel-%j.out
#SBATCH --error=logs/build-ad-panel-%j.err
# No --partition/--account here: those are site-specific. Pass them at submit
# time (sbatch -p <partition> -A <account> ...) or export SBATCH_PARTITION and
# SBATCH_ACCOUNT, which SLURM picks up automatically.
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
# hmmsearch may come from a module, a conda env or a container. Only touch the
# module system when hmmsearch is not already on PATH, so this runs unchanged on
# sites with different module names, or none at all.
if ! command -v hmmsearch >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
  module load "${HMMER_MODULE:-hmmer}" 2>/dev/null || true
fi
command -v hmmsearch >/dev/null 2>&1 || {
  echo "!! hmmsearch not found on PATH." >&2
  echo "   Load HMMER first (module load <your hmmer module>, or activate a conda" >&2
  echo "   env providing it), or set HMMER_MODULE to your site's module name." >&2
  exit 1; }

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
: "${DB:?set DB to a database root with room for the panel, e.g. export DB=/path/to/databases}"
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
MAP="$OUT/ad_panel_map.tsv" python3 - <<'PYCOV'
import csv, os
# Joins on marker_id -- no gene-name re-tokenisation here, so this cannot
# disagree with what the build actually recorded.
covered = {(r.get('marker_id') or '').strip()
           for r in csv.DictReader(open(os.environ['MAP']), delimiter='\t')}

by_module, missing = {}, []
for r in csv.DictReader(open(os.environ['PANEL']), delimiter='\t'):
    mod, mid = r['module'], (r.get('marker_id') or '').strip()
    ok = mid in covered
    n_tot, n_ok = by_module.get(mod, (0, 0))
    by_module[mod] = (n_tot + 1, n_ok + (1 if ok else 0))
    if not ok:
        missing.append((mod, r['gene'], (r.get('tier') or '').strip(), mid))

print(f"{'module':<28} {'markers':>8} {'with model':>11}")
for mod, (tot, ok) in sorted(by_module.items()):
    flag = '' if ok == tot else '   <-- gaps'
    print(f"{mod:<28} {tot:>8} {ok:>11}{flag}")
tot = sum(t for t, _ in by_module.values()); okc = sum(o for _, o in by_module.values())
print(f"\nTOTAL {okc}/{tot} panel markers resolved to >=1 NCBIfam model")
if missing:
    print(f"\nUNRESOLVED ({len(missing)}) -- need Pfam fallback, dbCAN, or panel fix:")
    for mod, gene, tier, mid in missing:
        print(f"  {mid:<24} {gene}" + (f"  (tier {tier})" if tier else ""))
PYCOV
echo
echo "================ D. EC-derived matches (REVIEW THESE) ================"
ECLIST="$WORK/ec_derived.txt"
awk -F'\t' 'NR>1 && $6 ~ /^ec:/ {print "  " $1 "\t" $3 "\t" $6 "\t" $2}' "$OUT/ad_panel_map.tsv" \
  | LC_ALL=C sort -k2,2 > "$ECLIST"
head -60 "$ECLIST"
echo "  ... EC-derived total: $(wc -l < "$ECLIST")"
echo
echo "================ E. CURATED accessions vs what NCBIfam says they are ================"
# Pinned rows have no symbol/EC fallback, so a wrong accession silently swaps in
# an unrelated protein and still reports as "resolved". Four of the original
# eight were wrong this way (mcrG -> a phosphonate transporter, mch -> catechol
# 2,3-dioxygenase, mtrA -> an EPS transport protein, fhs -> a glutaredoxin).
# Eyeball this table whenever the curated set changes.
awk -F'\t' '
  NR==FNR { if (FNR>1) prod[$1]=$11; next }              # hmm_PGAP.tsv: acc -> product_name
  FNR>1 && $6=="curated" { printf "  %-24s %-14s %-18s %s\n", $3, $1, $2, (($1 in prod) ? prod[$1] : "?? not in metadata") }
' "$WORK/hmm_PGAP.tsv" "$OUT/ad_panel_map.tsv" | LC_ALL=C sort
echo
echo "==> artifacts:"
ls -lh "$OUT"
echo "==> done $(date)"
