#!/usr/bin/env bash
# Fetch the Pfam HMMs named in the panel, on the LOGIN node.
#
# Pfam is CC0, so it is as licence-clean as NCBIfam and is the right fallback
# when NCBIfam has no usable model. Several markers are only reachable this way:
# NCBIfam wraps some Pfam families but does not ship the wrapper (mtmB's
# NF017207.7 is indexed in hmm_PGAP.tsv yet absent from the archive), so going to
# Pfam directly is not a workaround -- it is the primary source for those rows.
#
# Only accessions the panel explicitly pins are fetched. Pfam families are often
# superfamily-wide and would wreck specificity if pulled in wholesale; each one
# has to be opted into per row and its specificity measured, exactly as for
# NCBIfam accessions.
#
# CUTOFF NORMALISATION. The panel is searched with `hmmsearch --cut_nc` because
# that is NCBIfam's curated threshold. Pfam's curated threshold is GA, and its NC
# sits below GA (PF05369: GA 25, NC 22.7), so applying --cut_nc to an unmodified
# Pfam model would be looser than Pfam intends. On import we therefore set
# NC := GA, so one --cut_nc run applies each source's own intended cutoff. This
# is a deliberate, documented rewrite -- the alternative, a second search pass
# with --cut_ga, buys nothing and doubles the moving parts.
set -euo pipefail

: "${DB:?set DB to a database root, e.g. export DB=/path/to/databases}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PANEL="${PANEL:-$HERE/AD_methanogenesis_panel.tsv}"
OUT="${PFAMDIR:-$DB/pfam}"
API="${INTERPRO_API:-https://www.ebi.ac.uk/interpro/api}"

mkdir -p "$OUT/models"

accs=$(grep -o 'PF[0-9]\{5\}' "$PANEL" | sort -u)
[ -n "$accs" ] || { echo "==> no PF accessions pinned in the panel; nothing to do"; exit 0; }
echo "==> panel pins $(echo "$accs" | wc -w) Pfam accession(s)"

for acc in $accs; do
  dest="$OUT/models/${acc}.hmm"
  if [ -s "$dest" ]; then
    echo "    $acc already cached"
    continue
  fi
  tmp="$OUT/models/.${acc}.gz"
  if ! curl -s --fail --max-time 90 "$API/entry/pfam/${acc}/?annotation=hmm" -o "$tmp"; then
    echo "!! could not fetch $acc from InterPro" >&2; rm -f "$tmp"; continue
  fi
  # InterPro serves the HMM gzipped; tolerate it being served plain.
  if ! gunzip -c "$tmp" > "$dest.raw" 2>/dev/null; then cp "$tmp" "$dest.raw"; fi
  rm -f "$tmp"
  grep -q '^HMMER3' "$dest.raw" || { echo "!! $acc did not return an HMM" >&2; rm -f "$dest.raw"; continue; }

  # NC := GA, per the note above. Keep GA and TC untouched so the original
  # thresholds remain visible in the model.
  awk '
    /^GA / { ga=$2; sub(/;$/,"",ga) }
    /^NC / { if (ga != "") { printf "NC    %s %s;\n", ga, ga; next } }
    { print }
  ' "$dest.raw" > "$dest"
  rm -f "$dest.raw"
  printf "    %-10s %s  (GA %s -> NC)\n" "$acc" \
    "$(awk '/^DESC /{ $1=""; print substr($0,2,52) }' "$dest")" \
    "$(awk '/^GA /{ g=$2; sub(/;$/,"",g); print g; exit }' "$dest")"
done

cat "$OUT/models"/*.hmm > "$OUT/panel_pfam.hmm" 2>/dev/null || true
n=$(grep -c '^NAME ' "$OUT/panel_pfam.hmm" 2>/dev/null || echo 0)
echo
echo "==> $n Pfam model(s) in $OUT/panel_pfam.hmm"
echo "Next: rebuild the panel (build_ad_panel.sbatch) to merge them in."
