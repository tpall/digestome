#!/usr/bin/env bash
#SBATCH --job-name=ad-catalogue
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=logs/ad-catalogue-%j.out
#SBATCH --error=logs/ad-catalogue-%j.err
# Score the whole biogas MAG catalogue and aggregate it into one community
# profile. Pass -p/-A at submit time; nothing site-specific is hardcoded.
#
# ~1400 MAGs x (126 HMMs vs ~2500 proteins) is small per genome but adds up, so
# genomes are scored in parallel with xargs rather than serially.
set -euo pipefail

if ! command -v hmmsearch >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
  module load "${HMMER_MODULE:-hmmer}" 2>/dev/null || true
fi
command -v hmmsearch >/dev/null 2>&1 || {
  echo "!! hmmsearch not on PATH; load HMMER or set HMMER_MODULE" >&2; exit 1; }

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
  KEGGLESS="$SLURM_SUBMIT_DIR/kegg-less"
else
  KEGGLESS="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
: "${DB:?set DB to the database root}"

CAT="${CATALOGUE:-$DB/ad_panel_test/biogas_catalogue}"
PROT="$CAT/proteomes"
TAX="$CAT/catalogue_taxonomy.tsv"
OUTD="${CATOUT:-$CAT/scored}"
PANEL="$KEGGLESS/AD_methanogenesis_panel.tsv"
HMM="$DB/ad_panel/ad_panel.hmm"
MAP="$DB/ad_panel/ad_panel_map.tsv"

for f in "$HMM" "$MAP" "$PANEL" "$TAX"; do
  [ -s "$f" ] || { echo "!! missing $f" >&2; exit 1; }
done
mkdir -p "$OUTD"

NCPU="${SLURM_CPUS_PER_TASK:-8}"
PAR="${PAR:-$(( NCPU / 2 ))}"; [ "$PAR" -lt 1 ] && PAR=1
echo "==> host=$(hostname) job=${SLURM_JOB_ID:-none} $(date)"
echo "==> $(ls "$PROT"/*.faa 2>/dev/null | wc -l) proteomes; $PAR concurrent, 2 threads each"
echo "==> panel: $(grep -c '^NAME ' "$HMM") models"

export KEGGLESS PANEL HMM MAP TAX OUTD
score_one() {
  faa="$1"
  label=$(basename "$faa" .faa)
  # Skip work already done, so the job is resumable after a timeout.
  [ -s "$OUTD/${label}.summary.txt" ] && return 0
  hmmsearch --cpu 2 --cut_nc --tblout "$OUTD/${label}.tblout" "$HMM" "$faa" > /dev/null 2>&1 || return 0
  python3 "$KEGGLESS/panel_scored.py" \
      --tblout "$OUTD/${label}.tblout" --map "$MAP" --panel "$PANEL" \
      --name "$label" --gtdbtk "$TAX" \
      --out "$OUTD/${label}.modules.tsv" > "$OUTD/${label}.summary.txt" 2>/dev/null || true
  rm -f "$OUTD/${label}.tblout"     # ~1400 tblouts is a lot of small files
}
export -f score_one

ls "$PROT"/*.faa | xargs -P "$PAR" -I{} bash -c 'score_one "$@"' _ {}

echo "==> scored $(ls "$OUTD"/*.summary.txt 2>/dev/null | wc -l) genomes"
echo
python3 "$KEGGLESS/aggregate_community.py" \
    --dir "$OUTD" \
    --sample "Biogas microbiome catalogue (PRJNA602310)" \
    --description "Campanaro/Treu MAG catalogue: MHQ+HQ bins from 134 anaerobic digester metagenomes across multiple plants, feedstocks and temperatures. Taxonomy derived from NCBI organism names, not GTDB-Tk." \
    --out-json "$CAT/profile.json" --out-txt "$CAT/profile.txt"

echo
echo "==> done $(date)"
