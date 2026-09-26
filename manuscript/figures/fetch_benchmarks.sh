#!/usr/bin/env bash
# Copy the scored benchmark outputs the figures are drawn from off the cluster, into results/benchmarks/
# (gitignored). Only the reduced tables in manuscript/figures/data/ are committed.
#
#   bash manuscript/figures/fetch_benchmarks.sh
#
# The benchmarks were scored on the academic account (public data only). Each directory is named for
# the run the manuscript reports; see make_figure_data.py for which table each one feeds.
set -euo pipefail
cd "$(dirname "$0")/../.."
HOST=${HOST:-taavi74@login1.hpc.ut.ee}
SRC=${SRC:-/gpfs/space/projects/preterm/databases/ad_panel_test}
DST=results/benchmarks
mkdir -p "$DST/gtdb" "$DST/catalogue"

rsync -a "$HOST:$SRC/gtdb/"{scored,scored_2026-09-26,scored_contig_2026-09-26,gtdb_sample.tsv,gtdbtk_sample.tsv,pf05369_hits.txt} "$DST/gtdb/"
rsync -a "$HOST:$SRC/biogas_catalogue/"{scored_2026-09-26,catalogue_taxonomy.tsv} "$DST/catalogue/"

# Contigs per genome, counted from the Prodigal proteomes on the cluster: the proteomes themselves are
# large and only this count is needed. Prodigal names proteins <contig>_<n>.
ssh "$HOST" "cd $SRC/gtdb/proteomes_prodigal && for f in *.faa; do
  printf '%s\t%s\n' \"\${f%.faa}\" \"\$(grep '>' \"\$f\" | cut -d' ' -f1 | sed -E 's/^>//; s/_[0-9]+\$//' | sort -u | wc -l)\"
done" > "$DST/gtdb/contig_counts.tsv"

echo "fetched into $DST:"
du -sh "$DST"/*
