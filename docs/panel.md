# The panel: detection and curation

Which database answers which row, how to resolve a new accession, and what a single marker tells you. Maintainer-facing.

## Detection strategy (which tool for which row)
| Row type | Detect with | Why |
|---|---|---|
| Carbohydrate, protein and lipid hydrolysis | **Pfam catalytic families** (CC0) | curated one family per marker. CAZy/dbCAN and MEROPS are avoided: both are non-commercial in their public form |
| Peptidases / lipases | Pfam (avoid PROSITE/SMART, non-commercial) | family level is the right resolution for capability |
| Everything else | **NCBIfam (incl. TIGRFAM) gene HMMs** first, Pfam as fallback | NCBIfam HMMs are gene-specific + public domain |
| Reactions/EC cross-ref | **Rhea** (CC BY 4.0) via EC | reaction-level identifiers |

Results are reported in **gene symbol, EC, Pfam/NCBIfam and Rhea** identifiers. KEGG KO
codes and map images are not emitted, since they are KEGG content.

---

## Single-gene diagnostic markers (the quick fingerprint)
| Marker | Tells you |
|---|---|
| **mcrA** | total methanogen abundance; phylo-bin it to know *which* methanogens dominate (the #1 health signal) |
| **fhs** | acetogen / homoacetogenesis potential |
| **fwdB + mtrA** | hydrogenotrophic capacity |
| **cdhA + acs/ackA** | acetoclastic capacity (and Methanothrix vs Methanosarcina) |
| **mttB / mtaB** | methylamine / methanol (methylotrophic) capacity |
| **hydA / mvhA / frhA** | electron flow & H2 balance (instability context) |

**Interpretation for operators:** a healthy stable digester usually shows a balanced hydrogenotrophic + acetoclastic signal. A shift toward hydrogenotrophic dominance with rising VFA markers and falling acetoclastic (Methanothrix) signal is a classic **acidification / overload early-warning** pattern — that shift is the earliest actionable signal the panel produces.

---

## Resolver — auto-populate & validate accessions (free sources)
```bash
# 1. NCBIfam (gene-specific, public domain): HMMs + their EC/gene mapping
#    https://ftp.ncbi.nlm.nih.gov/hmm/current/   (hmm_PGAP.HMM + hmm_PGAP.tsv)

# 2. EC -> Rhea (CC BY 4.0) for reaction-level reporting:
#    https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv

# 3. EC/gene -> Pfam validation via InterPro API, e.g.:
curl -s 'https://www.ebi.ac.uk/interpro/api/entry/pfam/?search=methyl-coenzyme%20M%20reductase' | jq '.results[].metadata | {accession,name}'

# 4. Hydrolysis: curated Pfam catalytic families, pinned per marker. dbCAN/CAZy is
#    deliberately not used; its public terms are non-commercial.
```
