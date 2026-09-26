#!/usr/bin/env bash
# Copy the scored benchmark outputs the figures are drawn from off the cluster, into results/benchmarks/
# (gitignored). Only the reduced tables in manuscript/figures/data/ are committed.
#
#   bash manuscript/figures/fetch_benchmarks.sh
#
# The benchmarks were scored on Magrittr OÜ's commercial allocation at the University of Tartu HPC
# Centre, from inputs rebuilt there from public sources (GTDB r226 files, NCBI). Remote directory
# names are mapped to the local ones make_figure_data.py reads.
set -euo pipefail
cd "$(dirname "$0")/../.."
HOST=${HOST:-magrittr-hpc}
SRC=${SRC:-databases/ad_panel_test}
DST=results/benchmarks
mkdir -p "$DST/gtdb" "$DST/catalogue"

rsync -a "$HOST:$SRC/gtdb/"{gtdb_sample.tsv,gtdbtk_sample.tsv} "$DST/gtdb/"
rsync -a --delete "$HOST:$SRC/gtdb/scored_final/" "$DST/gtdb/final/"
rsync -a --delete "$HOST:$SRC/gtdb/scored_contig_final/" "$DST/gtdb/contig_final/"
rsync -a --delete "$HOST:$SRC/gtdb_aug11/scored/" "$DST/gtdb/family_level_aug11/"   # scorer at d2daee1
rsync -a "$HOST:$SRC/biogas_catalogue/catalogue_taxonomy.tsv" "$DST/catalogue/"
rsync -a --delete "$HOST:$SRC/biogas_catalogue/scored_final_run1/" "$DST/catalogue/final/"

# Contigs per genome, counted from the Prodigal proteomes on the cluster: the proteomes themselves are
# large and only this count is needed. Prodigal names proteins <contig>_<n>.
ssh "$HOST" "cd $SRC/gtdb/proteomes_prodigal && for f in *.faa; do
  printf '%s\t%s\n' \"\${f%.faa}\" \"\$(grep '>' \"\$f\" | cut -d' ' -f1 | sed -E 's/^>//; s/_[0-9]+\$//' | sort -u | wc -l)\"
done" > "$DST/gtdb/contig_counts.tsv"

echo "fetched into $DST:"
du -sh "$DST"/*
