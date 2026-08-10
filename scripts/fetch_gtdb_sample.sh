#!/usr/bin/env bash
# Fetch proteomes for a GTDB representative sample, on a LOGIN node (needs internet).
#
#   export DB=/path/to/databases
#   bash scripts/fetch_gtdb_sample.sh
#
# Input:  $WORK/gtdb_sample.tsv   from sample_gtdb_reps.py
# Output: $WORK/proteomes/<accession>.faa
#         $WORK/gtdbtk_sample.tsv  GTDB-Tk-shaped taxonomy for --gtdbtk
#
# Resolution goes through NCBI's assembly summary rather than the `datasets` CLI.
# The CLI has moved download paths more than once and is not present on every
# cluster; the summary files give an ftp_path per accession directly and need
# nothing installed. Note the summaries live under ASSEMBLY_REPORTS/, not at the
# top of genomes/, which is a 404.
#
# Idempotent: an accession whose .faa is already present is skipped, so an
# interrupted run resumes by being re-invoked.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
: "${DB:?set DB to a database root}"
WORK="${WORK:-$DB/ad_panel_test/gtdb}"
SAMPLE="$WORK/gtdb_sample.tsv"
PROT="$WORK/proteomes"
JOBS="${JOBS:-8}"                 # concurrent downloads; keep modest, NCBI is shared
BASE=https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS

[ -s "$SAMPLE" ] || { echo "!! no sample at $SAMPLE -- run sample_gtdb_reps.py first" >&2; exit 1; }
mkdir -p "$PROT" "$WORK/tmp"
cd "$WORK"

# ---- 1. GTDB-Tk-shaped taxonomy -------------------------------------------
# The acetoclastic gate is confirmed against lineage, so this is not optional.
awk -F'\t' 'NR==1 {print "user_genome\tclassification"; next} {print $1"\t"$4}' \
    "$SAMPLE" > "$WORK/gtdbtk_sample.tsv"
echo "==> taxonomy: $(( $(wc -l < "$WORK/gtdbtk_sample.tsv") - 1 )) genomes"

# ---- 2. accession -> ftp_path ---------------------------------------------
for f in assembly_summary_refseq.txt assembly_summary_genbank.txt; do
  [ -s "$WORK/${f/assembly_summary_/as_}" ] || \
    curl -sfL --retry 3 --retry-delay 5 "$BASE/$f" -o "$WORK/${f/assembly_summary_/as_}" || {
      echo "!! could not fetch $f" >&2; exit 1; }
done
if [ ! -s "$WORK/acc2ftp.tsv" ]; then
  cat as_refseq.txt as_genbank.txt \
    | awk -F'\t' '!/^#/ && $20!="na" && $20!="" {print $1"\t"$20}' | sort -u > acc2ftp.tsv
fi
awk -F'\t' 'NR>1 {print $1}' "$SAMPLE" | sort -u > tmp/wanted.txt
join -t"$(printf '\t')" tmp/wanted.txt acc2ftp.tsv > tmp/sample_ftp.tsv
echo "==> resolved $(wc -l < tmp/sample_ftp.tsv) of $(wc -l < tmp/wanted.txt) accessions to a download path"
comm -23 tmp/wanted.txt <(cut -f1 acc2ftp.tsv) > "$WORK/unresolved_accessions.txt"

# ---- 3. download -----------------------------------------------------------
# Not every assembly carries a protein set: GenBank entries without annotation
# have no _protein.faa.gz. A miss here is expected and is counted, not retried.
fetch_one() {
  local acc="$1" ftp="$2" out="$3"
  [ -s "$out/$acc.faa" ] && return 0
  local base; base=$(basename "$ftp")
  curl -sfL --max-time 180 --retry 2 "$ftp/${base}_protein.faa.gz" 2>/dev/null \
    | gunzip -c > "$out/$acc.faa.part" 2>/dev/null
  if [ -s "$out/$acc.faa.part" ]; then mv "$out/$acc.faa.part" "$out/$acc.faa"; else rm -f "$out/$acc.faa.part"; fi
}
export -f fetch_one

need=$(awk -v p="$PROT" -F'\t' '{ f=p"/"$1".faa"; if ((getline line < f) <= 0) print }' tmp/sample_ftp.tsv | wc -l)
echo "==> $need to download at concurrency $JOBS (this takes a while)"
awk -F'\t' -v p="$PROT" '{print $1"\t"$2"\t"p}' tmp/sample_ftp.tsv \
  | xargs -P "$JOBS" -n 3 bash -c 'fetch_one "$0" "$1" "$2"'

have=$(ls "$PROT"/*.faa 2>/dev/null | wc -l)
echo "==> $have proteomes in $PROT"
echo "==> $(wc -l < "$WORK/unresolved_accessions.txt") accessions had no assembly record"
echo "==> $(( $(wc -l < tmp/sample_ftp.tsv) - have )) resolved accessions carried no protein set"
