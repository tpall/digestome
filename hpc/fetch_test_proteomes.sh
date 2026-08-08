#!/usr/bin/env bash
# Fetch reference proteomes for the smoke test, on the LOGIN node.
#
# Accessions are pinned in tests/genomes.tsv so the test is reproducible; the
# exact assembly directory is resolved from the FTP listing because it carries
# the assembly name, which the accession alone does not give us.
#
# Idempotent: already-downloaded proteomes are left alone.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GENOMES="${GENOMES:-$HERE/tests/genomes.tsv}"
OUT="${PROTEOMES:-/gpfs/space/projects/preterm/databases/ad_panel_test/proteomes}"
BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all

mkdir -p "$OUT"

acc_dir() {                       # GCF_000970025.1 -> GCF/000/970/025
  local digits=${1#GCF_}; digits=${digits%%.*}
  printf 'GCF/%s/%s/%s' "${digits:0:3}" "${digits:3:3}" "${digits:6:3}"
}

n=0
while IFS=$'\t' read -r label acc organism role; do
  case "$label" in ''|\#*) continue;; esac
  dest="$OUT/${label}.faa"
  if [ -s "$dest" ]; then
    echo "==> $label already present ($(grep -c '^>' "$dest") proteins)"
    n=$((n+1)); continue
  fi
  echo "==> $label  $acc  ($organism)"
  listing=$(curl -s --fail --max-time 60 "$BASE/$(acc_dir "$acc")/") || {
    echo "!! could not list directory for $acc" >&2; exit 1; }
  dir=$(printf '%s' "$listing" | grep -o "${acc}_[^\"/<>]*" | head -1)
  [ -n "$dir" ] || { echo "!! no assembly directory matching $acc" >&2; exit 1; }
  url="$BASE/$(acc_dir "$acc")/${dir}/${dir}_protein.faa.gz"
  curl -s --fail --max-time 300 "$url" | gunzip -c > "$dest" || {
    echo "!! download failed: $url" >&2; rm -f "$dest"; exit 1; }
  echo "    $(grep -c '^>' "$dest") proteins -> $dest"
  n=$((n+1))
done < "$GENOMES"

echo
echo "==> $n proteomes in $OUT"
echo "Next: mkdir -p logs && sbatch kegg-less/hpc/smoke_test.sbatch"
