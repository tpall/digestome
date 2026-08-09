#!/usr/bin/env bash
# Fetch the Campanaro/Treu biogas microbiome MAG catalogue, on the LOGIN node.
#
#   Campanaro, Treu et al. (2020) "New insights from the biogas microbiome by
#   comprehensive genome-resolved metagenomics of nearly 1600 species originating
#   from multiple anaerobic digesters". Biotechnol Biofuels 13:25.
#   MHQ/HQ MAGs deposited under NCBI BioProject PRJNA602310.
#
# This is THE reference catalogue for the anaerobic-digestion microbiome: ~1441
# annotated MAGs assembled from 134 biogas reactor metagenomes. Unlike the
# handful of digester genomes findable by taxid, these are real reactor bins
# across many plants, feedstocks and temperatures -- the right basis both for a
# community report and for benchmarking the panel at scale.
#
# Proteomes come straight from the deposited annotation, so no gene calling.
#
# Writes:
#   $OUT/proteomes/<accession>.faa
#   $OUT/catalogue_taxonomy.tsv   GTDB-Tk-summary-SHAPED, built from NCBI organism
#                                 names -- see the warning written into its header.
set -euo pipefail

: "${DB:?set DB to a database root, e.g. export DB=/path/to/databases}"
OUT="${CATALOGUE:-$DB/ad_panel_test/biogas_catalogue}"
BIOPROJECT="${BIOPROJECT:-PRJNA602310}"
API=https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome

command -v datasets >/dev/null 2>&1 || {
  echo "!! the NCBI 'datasets' CLI is not on PATH." >&2
  echo "   conda activate an env providing it, or install from" >&2
  echo "   https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/" >&2
  exit 1; }

mkdir -p "$OUT/proteomes"
cd "$OUT"

# ---- 1. page through the report to collect accession + organism ----
if [ ! -s "$OUT/catalogue_index.tsv" ]; then
  echo "==> listing $BIOPROJECT"
  python3 - "$API" "$BIOPROJECT" "$OUT/catalogue_index.tsv" <<'PY'
import json, subprocess, sys
api, bioproject, dest = sys.argv[1], sys.argv[2], sys.argv[3]
rows, token, page = [], None, 0
while True:
    url = f"{api}/bioproject/{bioproject}/dataset_report?page_size=1000"
    if token:
        url += f"&page_token={token}"
    raw = subprocess.run(["curl", "-s", "--max-time", "180", url],
                         capture_output=True, text=True).stdout
    d = json.loads(raw)
    for r in d.get('reports', []):
        org = (r.get('organism') or {}).get('organism_name', '')
        st = r.get('assembly_stats') or {}
        rows.append((r['accession'], org, str(st.get('total_sequence_length', ''))))
    page += 1
    token = d.get('next_page_token')
    print(f"    page {page}: {len(rows)} cumulative", file=sys.stderr)
    if not token:
        break
# The bioproject carries both GenBank and RefSeq deposits of some genomes
# (GCA_012517695.1 and GCF_012517695.1 are one MAG, not two). Counting both
# double-weights those genomes in every community statistic, so keep one --
# RefSeq where it exists, since that is the curated copy.
by_id = {}
for a, o, b in rows:
    key = a.split('_', 1)[1].split('.')[0]
    prev = by_id.get(key)
    if prev is None or (a.startswith('GCF_') and not prev[0].startswith('GCF_')):
        by_id[key] = (a, o, b)
deduped = sorted(by_id.values())
n_dup = len(rows) - len(deduped)
with open(dest, 'w') as fh:
    fh.write("#accession\torganism\ttotal_bp\n")
    for a, o, b in deduped:
        fh.write(f"{a}\t{o}\t{b}\n")
print(f"    {len(rows)} assemblies indexed, {n_dup} GCA/GCF duplicates dropped, "
      f"{len(deduped)} kept", file=sys.stderr)
PY
fi
n_idx=$(( $(wc -l < "$OUT/catalogue_index.tsv") - 1 ))
echo "==> $n_idx assemblies in $BIOPROJECT"

# ---- 2. download the proteomes, in batches ----
# One request for all 1441 fails: the NCBI datasets CLI aborts mid-stream
# ("stream error ... INTERNAL_ERROR") after a few hundred MB. Batching keeps each
# request small enough to succeed, makes the whole fetch resumable, and lets a
# single flaky batch be retried without redoing the rest.
BATCH="${BATCH:-200}"
RETRIES="${RETRIES:-3}"
if [ ! -s "$OUT/proteomes/.complete" ]; then
  awk -F'\t' 'NR>1{print $1}' "$OUT/catalogue_index.tsv" > "$OUT/accessions.txt"
  rm -rf "$OUT/_batches"; mkdir -p "$OUT/_batches"
  split -l "$BATCH" -d -a 3 "$OUT/accessions.txt" "$OUT/_batches/b"
  n_batch=$(ls "$OUT/_batches" | wc -l)
  echo "==> downloading $n_idx proteomes in $n_batch batches of $BATCH"
  bi=0
  for bf in "$OUT/_batches"/b*; do
    bi=$((bi+1))
    # Already have every accession in this batch? then skip it.
    missing=0
    while read -r acc; do
      [ -s "$OUT/proteomes/${acc}.faa" ] || missing=1
    done < "$bf"
    [ "$missing" -eq 0 ] && { echo "    batch $bi/$n_batch already present"; continue; }

    ok=0
    for try in $(seq 1 "$RETRIES"); do
      rm -f "$OUT/_b.zip"
      if datasets download genome accession --inputfile "$bf" \
            --include protein --filename "$OUT/_b.zip" --no-progressbar 2>/dev/null \
         && unzip -tq "$OUT/_b.zip" >/dev/null 2>&1; then
        ok=1; break
      fi
      echo "    batch $bi attempt $try failed; retrying" >&2
      sleep $(( try * 5 ))
    done
    [ "$ok" -eq 1 ] || { echo "!! batch $bi failed after $RETRIES attempts" >&2; continue; }

    rm -rf "$OUT/_unzip"; mkdir -p "$OUT/_unzip"
    unzip -q -o "$OUT/_b.zip" -d "$OUT/_unzip"
    for d in "$OUT/_unzip"/ncbi_dataset/data/*/; do
      acc=$(basename "$d")
      [ -s "$d/protein.faa" ] && mv "$d/protein.faa" "$OUT/proteomes/${acc}.faa"
    done
    rm -rf "$OUT/_unzip" "$OUT/_b.zip"
    echo "    batch $bi/$n_batch -> $(ls "$OUT/proteomes"/*.faa 2>/dev/null | wc -l) proteomes so far"
  done
  rm -rf "$OUT/_batches"
  ls "$OUT/proteomes"/*.faa 2>/dev/null | wc -l > "$OUT/proteomes/.complete"
fi
echo "==> proteomes: $(ls "$OUT/proteomes"/*.faa 2>/dev/null | wc -l)"

# ---- 3. taxonomy table, shaped like a GTDB-Tk summary ----
# panel_scored.py --gtdbtk expects user_genome + classification. These MAGs were
# not classified with GTDB-Tk, so the lineage here is derived from the NCBI
# organism name and carries only genus and species. That is enough for the
# acetoclastic clade check, which matches at genus and family, but it is NOT a
# GTDB-Tk result and must not be presented as one.
python3 - "$OUT/catalogue_index.tsv" "$OUT/catalogue_taxonomy.tsv" <<'PY'
import re, sys
src, dest = sys.argv[1], sys.argv[2]
ARCHAEAL = ('Methano', 'Thermoplasma', 'Halo', 'Sulfolobus', 'Archaeo',
            'Nitroso', 'Bathyarch', 'Thermoproteus', 'Candidatus Methano')
out = []
for line in open(src):
    if line.startswith('#'):
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 2:
        continue
    acc, org = f[0], f[1]
    name = re.sub(r'^Candidatus\s+', '', org).strip()
    parts = name.split()
    genus = parts[0] if parts else ''
    if not genus:
        continue
    domain = 'Archaea' if any(genus.startswith(a.replace('Candidatus ', ''))
                              for a in ARCHAEAL) else 'Bacteria'
    lineage = f"d__{domain};g__{genus}"
    if len(parts) > 1:
        lineage += f";s__{name}"
    out.append((acc, lineage))
with open(dest, 'w') as fh:
    fh.write("#\tDERIVED FROM NCBI ORGANISM NAMES -- NOT GTDB-Tk OUTPUT.\n")
    fh.write("#\tGenus and species only; no family/order/class ranks. Sufficient for the\n")
    fh.write("#\tacetoclastic clade check (which matches genus and family), but do not\n")
    fh.write("#\tpresent these as GTDB-Tk classifications.\n")
    fh.write("user_genome\tclassification\n")
    for acc, lin in out:
        fh.write(f"{acc}\t{lin}\n")
print(f"    {len(out)} taxonomy rows -> {dest}", file=sys.stderr)
PY

echo
echo "==> catalogue ready in $OUT"
echo "Next: mkdir -p logs && sbatch -p <partition> kegg-less/hpc/score_catalogue.sbatch"
