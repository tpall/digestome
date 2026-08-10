#!/usr/bin/env bash
#SBATCH --job-name=ad-community
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=logs/ad-community-%j.out
#SBATCH --error=logs/ad-community-%j.err
# No --partition/--account here: those are site-specific. Pass them at submit
# time (sbatch -p <partition> -A <account> ...) or export SBATCH_PARTITION and
# SBATCH_ACCOUNT, which SLURM picks up automatically.
# Score a set of MAG proteomes and aggregate them into one community profile --
# the shape a client-facing "digester microbiome health check" would take.
#
#   PROT=<dir of *.faa> SAMPLE=<name> [GTDBTK=<summary.tsv>] sbatch community_report.sbatch
#
# Runs two sets by default:
#   digester_mags/  real MAGs from PRJEB31310 (small, genuinely from a digester)
#   proteomes/      the 7 reference genomes, as a labelled mock community
# The two are kept in separate profiles on purpose: mixing real bins with
# reference genomes into one "sample" would misrepresent both.
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

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  KEGGLESS="$SLURM_SUBMIT_DIR/kegg-less"
else
  KEGGLESS="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
[ -f "$KEGGLESS/aggregate_community.py" ] || {
  echo "!! run from the worktree root" >&2; exit 1; }

: "${DB:?set DB to a database root with room for the panel, e.g. export DB=/path/to/databases}"
PANEL="$KEGGLESS/AD_methanogenesis_panel.tsv"
HMM="$DB/ad_panel/ad_panel.hmm"
MAP="$DB/ad_panel/ad_panel_map.tsv"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

score_set() {                       # $1=proteome dir  $2=sample name  $3=out dir  $4=desc  $5=gtdbtk(optional)
  local prot="$1" sample="$2" out="$3" desc="$4" gtdb="${5:-}"
  mkdir -p "$out"
  local n=0
  for faa in "$prot"/*.faa; do
    [ -s "$faa" ] || continue
    local label; label=$(basename "$faa" .faa)
    hmmsearch --cpu "$THREADS" --cut_nc \
              --tblout "$out/${label}.tblout" "$HMM" "$faa" > /dev/null
    local gt=(); [ -n "$gtdb" ] && gt=(--gtdbtk "$gtdb")
    python3 "$KEGGLESS/panel_scored.py" \
        --tblout "$out/${label}.tblout" --map "$MAP" --panel "$PANEL" \
        --name "$label" --out "$out/${label}.modules.tsv" "${gt[@]}" \
        --secretion "$KEGGLESS/secretion_modules.tsv" \
        > "$out/${label}.summary.txt"
    n=$((n+1))
  done
  echo "  scored $n MAG(s) from $prot"
  python3 "$KEGGLESS/aggregate_community.py" \
      --dir "$out" --sample "$sample" --description "$desc" \
      --out-json "$out/profile.json" --out-txt "$out/profile.txt"
}

echo "==> host=$(hostname) job=${SLURM_JOB_ID:-none} $(date)"
echo "==> db=$HMM ($(grep -c '^NAME ' "$HMM") models)"
echo

echo "================ real digester MAGs (PRJEB31310) ================"
score_set "$DB/ad_panel_test/digester_mags" \
          "PRJEB31310 digester MAGs" \
          "$DB/ad_panel_test/community_digester" \
          "Real anaerobic-digester MAGs from PRJEB31310. No GTDB-Tk taxonomy available, so acetoclastic calls fall back to gene markers only."

echo
echo "================ reference mock community ================"
score_set "$DB/ad_panel_test/proteomes" \
          "Mock digester community" \
          "$DB/ad_panel_test/community_mock" \
          "MOCK community of 7 reference genomes, not a real sample. Used to show the report shape when methanogens of each route are present." \
          "$KEGGLESS/tests/gtdbtk_fixture.tsv"

echo
echo "==> done $(date)"
