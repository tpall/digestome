# Anaerobic-Digestion / Methanogenesis marker-gene panel

> **Part of DRAM's experimental KEGG-less workflow** (branch `feature/kegg-less-ad`).
> Goal: distill microbial metabolism using only
> commercially-licensable databases (Pfam CC0, NCBIfam public-domain, Rhea
> CC BY, dbCAN) — no KEGG/KOfam, no MetaCyc. This panel is the curated scaffold
> the KEGG-less distill step scores against.

A **license-clean** functional scaffold for biogas/AD microbiome analysis: which genes mark each step of anaerobic digestion, how to detect them with free HMMs, and how to score pathway completeness for a "digester microbiome health check" report. No KEGG, no MetaCyc.

Companion data file: **`AD_methanogenesis_panel.tsv`**.

---

## ⚠️ Read this before production use
- **Biology is authoritative** — gene symbols, EC numbers, pathway branches, and the module-completeness rules below are curated and correct.
- **Accessions need a validation pass.** Pfam/NCBIfam IDs drift between releases and a few in the TSV are best-known starters or shared domains (not gene-specific). Before billing a client, resolve/validate every `pfam`/`ncbifam_tigrfam` against your *actual* HMM database (see "Resolver" below). The well-established methanogenesis accessions (mcrABG, ftr, mch, mtd, mer, mtrA) are reliable; treat the rest as candidates.
- Some Pfam domains are **shared** across enzymes (e.g. PF00871 = acetate *and* butyrate kinase; PF00374 = several [NiFe]-hydrogenases; PF06253 = whole MttB superfamily). For those, gene-level calls need either NCBIfam gene-specific HMMs, sequence-similarity to references, or phylogenetic placement — not the Pfam hit alone.

---

## Detection strategy (which tool for which row)
| Row type | Detect with | Why |
|---|---|---|
| Carbohydrate hydrolysis (`HYDROL-CARB`) | **dbCAN / run_dbcan** | CAZy families, not single HMMs |
| Peptidases / lipases | Pfam (avoid PROSITE/SMART — NC license) | clan-level is fine for "potential" |
| Everything else | **NCBIfam (incl. TIGRFAM) gene HMMs** first, Pfam as fallback | NCBIfam HMMs are gene-specific + public domain |
| Reactions/EC cross-ref | **Rhea** (CC BY 4.0) via EC | for reaction-level reporting in client outputs |

Report in **gene symbol + EC + Pfam/NCBIfam + Rhea** identifiers. Never emit `K#####` KO codes or `mapXXXXX` images in a paid deliverable (KEGG IP).

---

## Module-completeness scoring (the "% pathway present" figure)
Score each module as `core genes found / core genes expected`. A branch is "present" when its gate condition is met. These replace KEGG module scoring with rules you own.

**Hydrogenotrophic methanogenesis (CO2+H2 -> CH4)** — gate: `mcrA` present AND ≥4 of the C1 carriers.
Core: `fwdB/fmdB, ftr, mch, mtd (or hmd), mer, mtrA, mcrA`  (+ `frhA`, `mvhA`, `hdrA/B` for electron supply)

**Acetoclastic methanogenesis (acetate -> CH4)** — gate: `mcrA` AND `cdhA` AND acetate activation.
Core: `(ackA + pta) OR acs`, `cdhA`, `cdhC`, `mtrA`, `mcrA`
Genus hint: `acs` (high-affinity) -> *Methanothrix/Methanosaeta*; `ackA+pta` (low-affinity) -> *Methanosarcina*.

**Methylotrophic methanogenesis (methanol/methylamines/methylsulfides -> CH4)** — gate: `mcrA` AND ≥1 substrate methyltransferase.
Core: `mtaB (methanol) OR mttB/mtbB/mtmB (methylamines) OR mtsA (DMS)`, corrinoid partner, `mcrA`
Note: H2-dependent methylotrophs (Methanomassiliicoccales) lack the Mtr/H4MPT C1 oxidation set — `mtrA` may be absent; don't penalise.

**Wood-Ljungdahl / homoacetogenesis (CO2 -> acetate)** — gate: `fhs` AND CODH/ACS.
Core: `fhs, folD, metF, acsB/cdhC, cooS/acsA`

**Syntrophic VFA oxidation** — gate: organism has the oxidation set AND a methanogenic/electron partner.
Butyrate: `bcd-etfAB, crt, hbd, thlA` (run oxidatively) — *Syntrophomonas*
Propionate: methylmalonyl-CoA set (`mcmA, mce`) — *Syntrophobacter*

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

**Interpretation for operators:** a healthy stable digester usually shows a balanced hydrogenotrophic + acetoclastic signal. A shift toward hydrogenotrophic dominance with rising VFA markers and falling acetoclastic (Methanothrix) signal is a classic **acidification / overload early-warning** pattern — that narrative is the sellable output.

---

## Resolver — auto-populate & validate accessions (free sources)
```bash
# 1. NCBIfam (gene-specific, public domain): get HMMs + their EC/gene mapping
#    https://ftp.ncbi.nlm.nih.gov/hmm/current/   (hmm_PGAP.HMM + hmm_PGAP.tsv)
#    Filter hmm_PGAP.tsv by the gene symbols / EC numbers in the panel,
#    then `hmmfetch` those models into a custom AD profile DB.

# 2. EC -> Rhea (CC BY 4.0) for reaction-level reporting:
#    https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv
#    join on the panel `ec` column to attach Rhea reaction IDs.

# 3. EC/gene -> Pfam validation via InterPro API, e.g.:
curl -s 'https://www.ebi.ac.uk/interpro/api/entry/pfam/?search=methyl-coenzyme%20M%20reductase' | jq '.results[].metadata | {accession,name}'

# 4. CAZymes: run_dbcan (dbCAN) on predicted proteins for HYDROL-CARB rows.
```
Build a custom `hmmsearch` DB from the validated NCBIfam models, run against MAG/contig proteins (Prodigal), then apply the module rules above. This is your DRAM-free, KEGG-free functional engine.

### Scripts
- **`build_ad_hmm_db.sh`** — downloads NCBIfam HMMs + `rhea2ec.tsv`, matches the panel (by gene symbol / EC), and emits `db/ad_panel.hmm` (pressed HMMER DB), `db/ad_panel_map.tsv` (model → gene/module), and `db/rhea2ec.tsv`.
- **`panel_scored.py`** — takes `hmmsearch --tblout` (+ optional `run_dbcan` overview) and the map, writes per-module completeness and prints a branch-gate + diagnostic-marker summary.

```bash
bash build_ad_hmm_db.sh
hmmsearch --cut_nc --tblout hits.tblout db/ad_panel.hmm proteins.faa
python3 panel_scored.py --tblout hits.tblout --map db/ad_panel_map.tsv \
  --panel AD_methanogenesis_panel.tsv --name MAG001 --out MAG001.modules.tsv
```
Both are **scaffolds** — validate matched accessions (gene-symbol matches are reliable; EC matches need review) and tune cutoffs before production.

---

## How this plugs into the pipeline
1. `aftekas` -> MAGs + GTDB-Tk taxonomy (who's there).
2. Prodigal proteins -> `hmmsearch` vs the validated AD panel HMMs + run_dbcan (what they can do).
3. Apply module-completeness rules -> per-MAG and community pathway profile.
4. Report: methanogen composition, methanogenesis branch balance, VFA/syntrophy markers, hydrolysis (CAZyme) capacity -> digester health narrative.

The curated panel is **your IP** — it is the defensible, license-clean core of the biogas offering.
