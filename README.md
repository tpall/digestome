# Anaerobic-Digestion / Methanogenesis marker-gene panel

> **Part of DRAM's experimental KEGG-less workflow** (branch `feature/kegg-less-ad`).
> Goal: distill microbial metabolism using only
> commercially-licensable databases (Pfam CC0, NCBIfam public-domain, Rhea
> CC BY, dbCAN) — no KEGG/KOfam, no MetaCyc. This panel is the curated scaffold
> the KEGG-less distill step scores against.

A **license-clean** functional scaffold for biogas/AD microbiome analysis: which genes mark each step of anaerobic digestion, how to detect them with free HMMs, and how to score pathway completeness for a "digester microbiome health check" report. No KEGG, no MetaCyc.

Companion data file: **`AD_methanogenesis_panel.tsv`**.

---

## Quick start

Four steps: fetch the source databases, build the HMM panel, score genomes, aggregate to a community profile.

```bash
# 0. Where the databases live. Required, and site-specific — nothing is
#    defaulted, so this must be set before anything else.
export DB=/path/to/databases

# 1. Download NCBIfam + Rhea (needs internet; run on a login node)
bash kegg-less/hpc/prefetch_ncbifam.sh

# 2. Build the pressed HMM panel -> $DB/ad_panel/
mkdir -p logs
sbatch -p <partition> kegg-less/hpc/build_ad_panel.sbatch

# 3. Score one genome
hmmsearch --cut_nc --tblout hits.tblout $DB/ad_panel/ad_panel.hmm proteins.faa
python3 kegg-less/panel_scored.py \
    --tblout hits.tblout \
    --map    $DB/ad_panel/ad_panel_map.tsv \
    --panel  kegg-less/AD_methanogenesis_panel.tsv \
    --name   MAG001 \
    --gtdbtk gtdbtk.bac120.summary.tsv \
    --out    MAG001.modules.tsv          # summary goes to stdout

# 4. Aggregate a directory of scored MAGs into one digester profile
python3 kegg-less/aggregate_community.py \
    --dir scored/ --sample digester_A \
    --out-json profile.json --out-txt profile.txt
```

Use `--cut_nc`: every model in the panel carries a curated NCBIfam noise cutoff,
and the build asserts this. Do not substitute a global `-T`; the curated cutoffs
across this panel span roughly 65 to 1400 bits.

**`--gtdbtk` is not optional in practice.** The acetoclastic gate is confirmed
against GTDB-Tk lineage, because ACDS/CODH is reversible and gene content alone
cannot distinguish acetyl-CoA cleavage (acetoclastic methanogenesis) from
acetyl-CoA synthesis (autotrophic carbon fixation). Without it, obligate
hydrogenotrophs carrying ACDS are called acetoclastic. See the
`Mthermautotroph` case in `tests/expectations.tsv`.

---

## Scripts

| script | what it does |
|---|---|
| `hpc/prefetch_ncbifam.sh` | Downloads NCBIfam HMMs + `rhea2ec.tsv`. Login node (compute nodes may lack outbound internet). Idempotent and resumable. |
| `build_ad_hmm_db.sh` | Matches the panel against NCBIfam, emits `ad_panel.hmm` (pressed), `ad_panel_map.tsv` (model → gene/module) and `rhea2ec.tsv`. Audits curated accessions and asserts every pressed model has an NC cutoff. |
| `panel_scored.py` | Per genome: `hmmsearch --tblout` (+ optional dbCAN, GTDB-Tk) → per-module completeness table on `--out`, gates + diagnostic markers on stdout. |
| `aggregate_community.py` | Many scored genomes → one community profile. Route balance, acetate-consumer presence, syntrophic capacity. JSON + text. |
| `audit_symbol_matches.py` | Build-time audit: flags panel rows whose gene symbol collides with an unrelated enzyme. Run when editing the panel. |

### `panel_scored.py`

| flag | |
|---|---|
| `--tblout` | **required.** `hmmsearch --tblout` output. |
| `--map` | **required.** `ad_panel_map.tsv` from the build. |
| `--panel` | **required.** `AD_methanogenesis_panel.tsv`. |
| `--out` | **required.** Per-module completeness TSV. |
| `--name` | Genome label for the report. Default `genome`. |
| `--gtdbtk` | GTDB-Tk `*.summary.tsv`. Repeatable (bac120 + ar53). Confirms the acetoclastic call and enables the genus hint. |
| `--dbcan` | `run_dbcan` `overview.txt`, satisfies the `HYDROL-CARB` module. Without it that module reports `NA` / `no-detector`. |
| `--evalue` | Optional extra E-value cut on top of `--cut_nc`. Off by default; the curated cutoffs are the intended filter. |

### `aggregate_community.py`

| flag | |
|---|---|
| `--dir` | **required.** Directory of `*.summary.txt` + `*.modules.tsv` from `panel_scored.py`. |
| `--sample` | Sample / digester name. |
| `--description` | What the sample actually is; carried into the report. |
| `--out-json` | **required.** Machine-readable profile. |
| `--out-txt` | Human-readable summary. |

Stdlib only, so it runs anywhere Python 3 does.

---

## Running on a cluster

The `hpc/` scripts are the intended entry points. They are site-neutral:

- **`DB` is required** and everything derives from it (`ncbifam/`, `ad_panel/`, `ad_panel_test/`). Nothing is defaulted to a particular filesystem.
- **No `--partition` or `--account`** in the sbatch headers. Pass them at submit time (`sbatch -p <partition> -A <account>`) or export `SBATCH_PARTITION` / `SBATCH_ACCOUNT`.
- **HMMER** is used from `PATH` if present (conda, container, spack); otherwise the module system is tried, with the module name overridable via `HMMER_MODULE`.

| script | run where | purpose |
|---|---|---|
| `hpc/prefetch_ncbifam.sh` | login | download NCBIfam + Rhea |
| `hpc/build_ad_panel.sbatch` | batch | build + press the panel |
| `hpc/fetch_test_proteomes.sh` | login | 7 reference proteomes for the smoke test |
| `hpc/smoke_test.sbatch` | batch | end-to-end run + assertions |
| `hpc/fetch_digester_mags.sh` | login | real digester MAGs (PRJEB31310) |
| `hpc/community_report.sbatch` | batch | score a MAG set + aggregate to a profile |

---

## Tests

`hpc/smoke_test.sbatch` runs the full path on seven reference proteomes and
asserts against published biology, so a failure means a real disagreement
rather than a drifted snapshot.

```bash
export DB=/path/to/databases
bash kegg-less/hpc/fetch_test_proteomes.sh
mkdir -p logs && sbatch -p <partition> kegg-less/hpc/smoke_test.sbatch
```

- `tests/genomes.tsv` — the 7 genomes and their accessions
- `tests/expectations.tsv` — the assertions; only well-established biology is asserted, anything uncertain is reported but not enforced
- `tests/gtdbtk_fixture.tsv` — lineages for the reference set, so the taxonomy-confirmed gate is exercised without running GTDB-Tk
- `tests/check_expectations.py` — the checker
- `tests/digester_mags.tsv` — real digester MAGs for `community_report.sbatch`

The set is chosen to be falsifiable: *M. barkeri* (generalist), *M. soehngenii*
(obligate acetoclastic), *M. bourgensis* and *M. thermautotrophicus* (obligate
hydrogenotrophs), *C. ljungdahlii* (acetogen), *S. wolfei* (syntroph), and
**E. coli as a specificity control** — it carries `pta`, `ackA`, `folD`, `atoB`
and the full `eut` operon, so it exercises nearly every gene-symbol collision
fixed on this branch. Any methanogenesis gate firing on E. coli is a live bug.

---

## Maturity

The biology is curated and the panel is tested, but read the boundaries:

- **Gene symbols, EC numbers, pathway branches and the module rules are authoritative.** They are curated, not inferred.
- **Curated accessions are audited at build time.** Four wrong accessions were corrected, nine colliding gene symbols pinned, and the build asserts every pressed model carries a curated NC cutoff. This is no longer an open validation task.
- **Non-curated (symbol- and EC-matched) rows are audited but weaker.** `audit_symbol_matches.py` flags collisions; gene-symbol matches are reliable, EC matches still warrant review before a client deliverable.
- **Some Pfam domains are shared across enzymes** (PF00871 = acetate *and* butyrate kinase; PF00374 = several [NiFe]-hydrogenases; PF06253 = the whole MttB superfamily). Gene-level calls on those need the NCBIfam gene-specific model, not the Pfam hit alone.
- **Coverage is honest about its gaps.** Modules and gates with no license-clean detector report `NA` / `NOT ASSESSABLE`, never `0.0`.

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

### The `scoring` column — how a module is (or is not) scored
Detection policy lives in the panel's `scoring` column, not in module names hardcoded in `panel_scored.py`. Per module:

| mode | meaning | reported as |
|---|---|---|
| `scored` | default — HMM scoring off the core tier | `pct_complete` |
| `assumed` | universal central metabolism; present by assumption, no discriminating marker exists | `NA`, status `assumed-present` |
| `dbcan` | satisfied by any CAZyme call from `run_dbcan` | `pct_complete`, or `NA`/`no-detector` if `--dbcan` was not passed |
| `pending` | no license-clean detector wired yet | `NA`, status `not-wired` |

**A module that cannot be scored reports `NA`, never `0.0`.** Emitting `0.0` made "no marker is defined for this module" indistinguishable from "this pathway is absent from the genome" — on a digester health report that reads as a failing plant. `panel_scored.py` prints the unscored modules and their reason under the branch gates so they cannot be missed.

Current non-default modules: `ACID-GLYC` = `assumed` (Embden-Meyerhof-Parnas is in essentially every organism), `HYDROL-CARB` = `dbcan`, `HYDROL-PROT` + `HYDROL-LIP` = `pending` (peptidases/lipases would need MEROPS, whose commercial terms are unverified).

A module whose rows are **all accessory-tier** (`ACID-ETOH`, `ACID-LACT`, `DIET`, `ENERGY`) also reports `NA`, status `no-core-tier`: scoring reads only the core tier today. Do not "fix" this by falling back to the accessory tier until `ahaA` is curated — it currently matches ~120 ATP-synthase models and would push `ENERGY` to a false ~100% in any genome.

**Hydrogenotrophic methanogenesis (CO2+H2 -> CH4)** — gate: `mcrA` present AND ≥4 of the C1 carriers.
Core: `fwdB/fmdB, ftr, mch, mtd (or hmd), mer, mtrA, mcrA`  (+ `frhA`, `mvhA`, `hdrA/B` for electron supply)

**Acetoclastic methanogenesis (acetate -> CH4)** — gate: `mcrA` AND `cdhA` AND acetate activation, **confirmed against GTDB-Tk lineage**.
Core: `(ackA + pta) OR acs`, `cdhA`, `cdhC`, `mtrA`, `mcrA`
Genus hint (gated behind the confirmed call): `acs` (high-affinity) -> *Methanothrix/Methanosaeta*; `ackA+pta` (low-affinity) -> *Methanosarcina*.

**Methylotrophic methanogenesis (methanol/methylamines/methylsulfides -> CH4)** — gate: `mcrA` AND ≥1 substrate methyltransferase.
Core: `mtaB (methanol) OR mttB/mtbB/mtmB (methylamines) OR mtsA (DMS)`, corrinoid partner, `mcrA`
Note: H2-dependent methylotrophs (Methanomassiliicoccales) lack the Mtr/H4MPT C1 oxidation set — `mtrA` may be absent; don't penalise.

**Wood-Ljungdahl / homoacetogenesis (CO2 -> acetate)** — gate: `fhs` AND CODH/ACS.
Core: `fhs, folD, metF, acsB/cdhC, cooS/acsA`

**Syntrophic VFA oxidation** — gate: organism has the oxidation set AND a methanogenic/electron partner.
Butyrate: `bcd-etfAB, crt, hbd, thlA` (run oxidatively) — *Syntrophomonas*. `bcd` and `crt` have no license-clean model, so this gate reports **NOT ASSESSABLE** rather than absent.
Propionate: methylmalonyl-CoA set (`mcmA, mce`) — *Syntrophobacter*

---

### Substrate signals (not methanogenesis branches)
**Ethanolamine utilisation (`ACID-ETA`)** — gate: `eutB` AND `eutC`, both subunits of the ammonia-lyase.
Core: `eutB, eutC` · accessory: `eutA` (reactivase), `eutD` (eut-specific phosphotransacetylase), `eutM` (microcompartment shell)

Read this as **two signals at once**: ethanolamine catabolism releases NH₃ stoichiometrically with acetaldehyde, so it is an ammonia-load indicator — and ammonia inhibition is a primary digester failure mode. It also points at feedstock, since ethanolamine comes from phosphatidylethanolamine, abundant in manure and food waste.

Every row here is **pinned**, and that is a requirement rather than a preference: most NCBIfam models whose product name mentions "ethanolamine" are phosphoethanolamine–lipid A transferases (MCR colistin resistance), entirely unrelated. Name or symbol matching would pull antibiotic-resistance genes into the panel.

`eutD` deserves its own note. It catalyses the same reaction as housekeeping `pta` and was originally matching the `pta` markers through the shared EC 2.3.1.8. It is excluded there and given a marker here instead, so a hit reads as ethanolamine catabolism rather than inflating acetate activation — which matters because the acetoclastic gate requires `ackA AND pta`.

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
# 1. NCBIfam (gene-specific, public domain): HMMs + their EC/gene mapping
#    https://ftp.ncbi.nlm.nih.gov/hmm/current/   (hmm_PGAP.HMM + hmm_PGAP.tsv)

# 2. EC -> Rhea (CC BY 4.0) for reaction-level reporting:
#    https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv

# 3. EC/gene -> Pfam validation via InterPro API, e.g.:
curl -s 'https://www.ebi.ac.uk/interpro/api/entry/pfam/?search=methyl-coenzyme%20M%20reductase' | jq '.results[].metadata | {accession,name}'

# 4. CAZymes: run_dbcan (dbCAN) on predicted proteins for HYDROL-CARB rows.
```

## How this plugs into the pipeline
1. `aftekas` -> MAGs + GTDB-Tk taxonomy (who's there).
2. Prodigal proteins -> `hmmsearch` vs the panel + run_dbcan (what they can do).
3. `panel_scored.py` per MAG -> module completeness + gates.
4. `aggregate_community.py` -> per-sample digester profile.
5. Report: methanogen composition, branch balance, VFA/syntrophy markers, hydrolysis capacity -> digester health narrative.

The curated panel is **your IP** — it is the defensible, license-clean core of the biogas offering.
