#!/usr/bin/env bash
#SBATCH --job-name=ad-panel-smoke
#SBATCH --partition=main
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=logs/ad-panel-smoke-%j.out
#SBATCH --error=logs/ad-panel-smoke-%j.err
# End-to-end smoke test: reference proteomes -> hmmsearch -> panel_scored.py,
# checked against known biology.
#
# Reference genomes rather than biogas MAGs on purpose. A MAG has no ground
# truth -- asserting anything about it validates the panel against its own
# output. These organisms have published metabolism, so the test can actually
# fail. E. coli is the specificity control: it carries pta, ackA, folD, atoB and
# the eut operon, so it exercises nearly every symbol collision fixed on this
# branch, and any methanogenesis gate firing on it is a live bug.
#
# Prereqs (login node):
#   bash kegg-less/hpc/fetch_test_proteomes.sh
#   mkdir -p logs
# Submit from the worktree root:
#   sbatch kegg-less/hpc/smoke_test.sbatch
set -euo pipefail
module purge
module load hmmer/3.3.2

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  KEGGLESS="$SLURM_SUBMIT_DIR/kegg-less"
else
  KEGGLESS="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
[ -f "$KEGGLESS/panel_scored.py" ] || {
  echo "!! cannot find panel_scored.py under $KEGGLESS -- submit from the worktree root" >&2
  exit 1; }

DB="${DB:-/gpfs/space/projects/preterm/databases}"
PANEL="$KEGGLESS/AD_methanogenesis_panel.tsv"
HMM="$DB/ad_panel/ad_panel.hmm"
MAP="$DB/ad_panel/ad_panel_map.tsv"
PROT="${PROTEOMES:-$DB/ad_panel_test/proteomes}"
WORK="${SMOKE_WORK:-$DB/ad_panel_test/run}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

for f in "$HMM" "$MAP" "$PANEL"; do
  [ -s "$f" ] || { echo "!! missing $f -- run build_ad_panel.sbatch first" >&2; exit 1; }
done
mkdir -p "$WORK"

echo "==> host=$(hostname) job=${SLURM_JOB_ID:-none} $(date)"
echo "==> db=$HMM ($(grep -c '^NAME ' "$HMM") models)"
echo

fails=0; n=0
while IFS=$'\t' read -r label acc organism role; do
  case "$label" in ''|\#*) continue;; esac
  faa="$PROT/${label}.faa"
  [ -s "$faa" ] || { echo "!! missing proteome $faa -- run fetch_test_proteomes.sh" >&2; exit 1; }
  n=$((n+1))
  echo "================ $label — $organism ================"
  echo "  ($role)"

  # --cut_nc applies each model's own NCBIfam-curated cutoff. Never substitute a
  # global -T here: curated cutoffs across this DB span ~65 to ~1400 bits.
  hmmsearch --cpu "$THREADS" --cut_nc \
            --tblout "$WORK/${label}.tblout" "$HMM" "$faa" > /dev/null

  python3 "$KEGGLESS/panel_scored.py" \
      --tblout "$WORK/${label}.tblout" --map "$MAP" --panel "$PANEL" \
      --name "$label" --out "$WORK/${label}.modules.tsv" \
      --gtdbtk "$KEGGLESS/tests/gtdbtk_fixture.tsv" \
      | tee "$WORK/${label}.summary.txt"

  echo "  --- assertions ---"
  if python3 "$KEGGLESS/tests/check_expectations.py" \
        --label "$label" --summary "$WORK/${label}.summary.txt" \
        --expect "$KEGGLESS/tests/expectations.tsv"; then :; else
    fails=$((fails+1))
  fi
  echo
done < "$KEGGLESS/tests/genomes.tsv"

echo "================ RESULT ================"
if [ "$fails" -eq 0 ]; then
  echo "  all $n genomes passed their assertions"
else
  echo "  $fails of $n genomes FAILED"
  echo "  $fails of $n genomes FAILED" >&2
fi
echo "==> artifacts in $WORK"
echo "==> done $(date)"
exit $([ "$fails" -eq 0 ] && echo 0 || echo 1)
