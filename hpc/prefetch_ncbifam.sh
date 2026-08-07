#!/usr/bin/env bash
# Fetch NCBIfam (public domain) + Rhea (CC BY) on the LOGIN node.
#
# Compute nodes are not guaranteed outbound internet, so the downloads happen
# here and build_ad_hmm_db.sh finds them already in place -- its `[ -s file ] ||
# curl` guards make the sbatch step fully offline. Idempotent and resumable
# (curl -C -), so re-running after an interrupted transfer is safe.
set -euo pipefail

DB="${DB:-/gpfs/space/projects/preterm/databases}"
WORK="${WORK:-$DB/ncbifam/_work}"
OUT="${OUT:-$DB/ad_panel}"
NCBIFAM_BASE="${NCBIFAM_BASE:-https://ftp.ncbi.nlm.nih.gov/hmm/current}"
RHEA_EC_URL="${RHEA_EC_URL:-https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv}"

mkdir -p "$WORK" "$OUT"

echo "==> NCBIfam models"
[ -s "$WORK/hmm_PGAP.HMM.tgz" ] || curl -L --fail -C - -o "$WORK/hmm_PGAP.HMM.tgz" "$NCBIFAM_BASE/hmm_PGAP.HMM.tgz"
echo "==> NCBIfam metadata"
[ -s "$WORK/hmm_PGAP.tsv" ]     || curl -L --fail -C - -o "$WORK/hmm_PGAP.tsv"     "$NCBIFAM_BASE/hmm_PGAP.tsv"
echo "==> EC -> Rhea map"
[ -s "$OUT/rhea2ec.tsv" ]       || curl -L --fail -o "$OUT/rhea2ec.tsv" "$RHEA_EC_URL"

echo
echo "==> Fetched:"
ls -lh "$WORK/hmm_PGAP.HMM.tgz" "$WORK/hmm_PGAP.tsv" "$OUT/rhea2ec.tsv"
echo
echo "Next: mkdir -p logs && sbatch kegg-less/hpc/build_ad_panel.sbatch"
