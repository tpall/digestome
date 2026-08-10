#!/usr/bin/env bash
# Download genome sequences for the GTDB sample and call genes with Prodigal.
#
#   export DB=/path/to/databases
#   bash scripts/call_gtdb_genes.sh            # login node: needs internet
#
# Why not just download NCBI's protein sets, which is one step instead of two:
#
#   1. Coverage. Roughly a third of GTDB representatives are GenBank assemblies
#      with no NCBI annotation and therefore no protein file. The gap is not
#      random: it falls on uncultured, MAG-derived lineages, which are exactly
#      the ones a phylum-stratified sample exists to include. Fetching proteomes
#      directly lost 28 phyla outright and 44% of Methanomassiliicoccales.
#
#   2. Comparability. The digester MAG catalogue was gene-called with Prodigal.
#      Scoring one benchmark from Prodigal calls and the other from NCBI PGAP
#      compares annotation pipelines as much as genomes.
#
#   3. Validity of the MtmB test. PGAP does not emit pyrrolysine proteins at all:
#      in M. barkeri it yields no monomethylamine methyltransferase, so PF05369
#      scores zero regardless of the genome. A specificity test for that family
#      run on PGAP proteomes would return a meaningless zero everywhere.
#
# Output: $WORK/proteomes_prodigal/<accession>.faa
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
: "${DB:?set DB to a database root}"
WORK="${WORK:-$DB/ad_panel_test/gtdb}"
FNA="$WORK/genomes"
OUT="$WORK/proteomes_prodigal"
JOBS="${JOBS:-8}"
TABLE="${TABLE:-11}"           # 11 = standard bacterial/archaeal. See the amber note below.

cd "$WORK"
[ -s tmp/sample_ftp.tsv ] || { echo "!! run fetch_gtdb_sample.sh first (needs tmp/sample_ftp.tsv)" >&2; exit 1; }
mkdir -p "$FNA" "$OUT"

command -v prodigal >/dev/null 2>&1 || module load prodigal/2.6.3 2>/dev/null || true
command -v prodigal >/dev/null 2>&1 || { echo "!! prodigal not available" >&2; exit 1; }

# Genome sequence, then genes. The .fna is removed once called: 3,000 genomes of
# nucleotide sequence is tens of GB and nothing downstream reads it again.
#
# Translation table 11 treats TAG as a stop, so pyrrolysine proteins come out
# split at the amber codon. Both fragments still clear the Pfam threshold for
# MtmB, so this is sufficient for detection; TABLE=6 reads TAG through and
# reconstructs full-length proteins, at the cost of over-extending ordinary
# genes. Detection is the question here, so the default stays at 11.
call_one() {
  local acc="$1" ftp="$2" fna="$3" out="$4" table="$5"
  [ -s "$out/$acc.faa" ] && return 0
  local base; base=$(basename "$ftp")
  if [ ! -s "$fna/$acc.fna" ]; then
    curl -sfL --max-time 300 --retry 2 "$ftp/${base}_genomic.fna.gz" 2>/dev/null \
      | gunzip -c > "$fna/$acc.fna.part" 2>/dev/null
    [ -s "$fna/$acc.fna.part" ] && mv "$fna/$acc.fna.part" "$fna/$acc.fna" || { rm -f "$fna/$acc.fna.part"; return 0; }
  fi
  # -p meta for short/fragmented assemblies, which single-genome mode handles badly
  local size; size=$(stat -c %s "$fna/$acc.fna" 2>/dev/null || echo 0)
  local mode=single; [ "$size" -lt 200000 ] && mode=meta
  prodigal -i "$fna/$acc.fna" -a "$out/$acc.faa.part" -p "$mode" -g "$table" -q >/dev/null 2>&1
  if [ -s "$out/$acc.faa.part" ]; then mv "$out/$acc.faa.part" "$out/$acc.faa"; else rm -f "$out/$acc.faa.part"; fi
  rm -f "$fna/$acc.fna"
}
export -f call_one

todo=$(awk -v o="$OUT" -F'\t' '{ if (system("[ -s "o"/"$1".faa ]") != 0) print }' tmp/sample_ftp.tsv | wc -l)
echo "==> $(wc -l < tmp/sample_ftp.tsv) accessions, $todo still to call, concurrency $JOBS"

awk -F'\t' -v f="$FNA" -v o="$OUT" -v t="$TABLE" '{print $1"\t"$2"\t"f"\t"o"\t"t}' tmp/sample_ftp.tsv \
  | xargs -P "$JOBS" -n 5 bash -c 'call_one "$0" "$1" "$2" "$3" "$4"'

have=$(ls "$OUT"/*.faa 2>/dev/null | wc -l)
echo "==> $have proteomes in $OUT"
echo "==> $(( $(wc -l < tmp/sample_ftp.tsv) - have )) accessions produced no proteome"
rmdir "$FNA" 2>/dev/null || echo "==> note: $FNA still holds $(ls "$FNA" | wc -l) uncalled .fna files"
