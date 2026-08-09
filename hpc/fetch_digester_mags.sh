#!/usr/bin/env bash
# Fetch real anaerobic-digester MAGs (proteomes) on the LOGIN node.
#
# These are genuine digester bins from PRJEB31310, "Investigating the response of
# digester metagenomes to various short-chain fatty acids using genome-centric
# metatranscriptomics" -- as opposed to the reference proteomes in
# fetch_test_proteomes.sh, which exist to assert known biology.
#
# NCBI carries only a handful of ANNOTATED genome-sized assemblies under taxid
# 1263854 (anaerobic digester metagenome); most entries there are whole
# metagenome assemblies of 300+ Mbp, not individual MAGs. So this is a small
# real sample, useful for exercising the community report on true digester data
# rather than for claiming community-wide coverage.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LIST="${MAGLIST:-$HERE/tests/digester_mags.tsv}"
OUT="${MAGPROT:-/gpfs/space/projects/preterm/databases/ad_panel_test/digester_mags}"
BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all

mkdir -p "$OUT"

acc_dir() {                        # GCA_900750605.1 -> GCA/900/750/605
  local pre=${1%%_*}               # GCA or GCF
  local digits=${1#*_}; digits=${digits%%.*}
  printf '%s/%s/%s/%s' "$pre" "${digits:0:3}" "${digits:3:3}" "${digits:6:3}"
}

n=0
while IFS=$'\t' read -r label acc note; do
  case "$label" in ''|\#*) continue;; esac
  dest="$OUT/${label}.faa"
  if [ -s "$dest" ]; then
    echo "==> $label already present ($(grep -c '^>' "$dest") proteins)"
    n=$((n+1)); continue
  fi
  echo "==> $label  $acc"
  listing=$(curl -s --fail --max-time 60 "$BASE/$(acc_dir "$acc")/") || {
    echo "!! could not list directory for $acc" >&2; exit 1; }
  dir=$(printf '%s' "$listing" | grep -o "${acc}_[^\"/<>]*" | head -1)
  [ -n "$dir" ] || { echo "!! no assembly directory matching $acc" >&2; exit 1; }
  url="$BASE/$(acc_dir "$acc")/${dir}/${dir}_protein.faa.gz"
  if curl -s --fail --max-time 300 "$url" | gunzip -c > "$dest" 2>/dev/null && [ -s "$dest" ]; then
    echo "    $(grep -c '^>' "$dest") proteins -> $dest"
  else
    rm -f "$dest"
    echo "!! no protein FASTA for $acc (unannotated assembly); skipping" >&2
    continue
  fi
  n=$((n+1))
done < "$LIST"

echo
echo "==> $n digester MAG proteomes in $OUT"
