# digestome

**Licence-clean functional profiling of anaerobic digestion microbiomes.**

> **Two entry points.** The command-line path needs only HMMER and Python 3
> (standard library), with no licensed database at any stage. A Nextflow workflow
> wraps the same tools for scale and executor portability. Inspired by DRAM's idea
> of distilling annotations into pathway-level statements, but it shares no code
> with DRAM and does not depend on it.

A licence-clean functional scaffold for biogas and anaerobic digestion microbiome
analysis: which genes mark each step, how to detect them with freely usable HMMs,
and how to score pathway completeness for a digester health-check report. No KEGG,
no MetaCyc, no CAZy/dbCAN.

Companion data file: **`assets/AD_methanogenesis_panel.tsv`**.

## Licence

Code is MIT (`LICENSE`). The curated data files, including the panel itself, are
CC BY 4.0 (`LICENSE-DATA`): the panel is a dataset rather than software, and
attribution is a licence term rather than a hope. No upstream database is
redistributed here; the build scripts fetch NCBIfam, Pfam and Rhea from source.

---

## Layout

```
bin/         panel_scored.py, aggregate_community.py — the analysis tools.
             Nextflow adds this to PATH in every task, so the workflow calls
             them by name and needs no path into the project directory.
scripts/     one-time setup: database staging, panel construction, curation
             audits, and sbatch wrappers. None requires a scheduler.
manuscript/  MANUSCRIPT.qmd and its bibliography.
tests/       expectations and fixtures.
assets/      the curated panel and the secretion modules — the data this
             project exists to provide. Both are CC BY 4.0; see LICENSE-DATA.
```

## Quick start

Two ways to run the same analysis. Either needs the panel built once first.

### 0. Build the panel (once)

Staging stays as scripts: it is a one-time operation that needs outbound network
access, and it gains nothing from a task graph.

```bash
export DB=/path/to/databases          # site-specific, required, never defaulted

bash scripts/prefetch_ncbifam.sh          # NCBIfam + Rhea   (login node: needs internet)
bash scripts/prefetch_pfam.sh             # pinned Pfam families
mkdir -p logs && sbatch -p <partition> scripts/build_ad_panel.sbatch
```

That writes `$DB/ad_panel/` with the pressed HMM database and the model map.

### Option A: Nextflow, for many genomes

Per-genome parallelism, resume, and executor independence: the same command runs
on a laptop, a SLURM cluster or a cloud batch service.

```bash
nextflow run . \
    --proteomes 'mags/*.faa' \
    --db        $DB \
    --gtdbtk    gtdbtk.bac120.summary.tsv \
    --sample_name 'digester A' \
    --outdir    results \
    -profile    slurm            # or: standard | conda | singularity | docker
```

Writes `results/profile.txt` and `results/profile.json`, plus a Nextflow trace,
timeline and report under `results/pipeline_info/`.

A samplesheet works instead of a glob when names matter:

```bash
nextflow run . --input samples.csv --db $DB     # columns: sample,faa
```

To check the wiring without tools, data or a database:

```bash
nextflow run . -profile test -stub-run
```

### Option B: command line, for a handful of genomes

No workflow engine, no Java. HMMER and Python 3 standard library is the whole
dependency list, and this is exactly what the workflow calls underneath.

```bash
# 1. search one proteome against the panel
hmmsearch --cut_nc --tblout MAG001.tblout $DB/ad_panel/ad_panel.hmm proteins.faa

# 2. score it
python3 bin/panel_scored.py \
    --tblout    MAG001.tblout \
    --map       $DB/ad_panel/ad_panel_map.tsv \
    --panel     assets/AD_methanogenesis_panel.tsv \
    --secretion assets/secretion_modules.tsv \
    --gtdbtk    gtdbtk.bac120.summary.tsv \
    --name      MAG001 \
    --out       MAG001.modules.tsv          # summary goes to stdout

# 3. aggregate a directory of scored genomes into one community profile
python3 bin/aggregate_community.py \
    --dir scored/ --sample 'digester A' \
    --out-json profile.json --out-txt profile.txt
```

### Which to use

Option B for a few genomes, for debugging, or where installing a workflow engine
is not worth it. Option A once there are enough genomes that parallelism and
resume matter, or when the same analysis has to move between a cluster and a
single machine. They produce identical results because they run the same code.

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

## The tools

| script | what it does |
|---|---|
| `scripts/prefetch_ncbifam.sh` | Downloads NCBIfam HMMs + `rhea2ec.tsv`. Login node (compute nodes may lack outbound internet). Idempotent and resumable. |
| `scripts/build_ad_hmm_db.sh` | Matches the panel against NCBIfam, emits `ad_panel.hmm` (pressed), `ad_panel_map.tsv` (model → gene/module) and `rhea2ec.tsv`. Audits curated accessions and asserts every pressed model has an NC cutoff. |
| `bin/panel_scored.py` | Per genome: `hmmsearch --tblout` (+ optional dbCAN, GTDB-Tk) → per-module completeness table on `--out`, gates + diagnostic markers on stdout. |
| `bin/aggregate_community.py` | Many scored genomes → one community profile. Route balance, acetate-consumer presence, syntrophic capacity. JSON + text. |
| `main.nf`, `nextflow.config` | Nextflow workflow over the analysis path: HMMSEARCH, SCORE, AGGREGATE. Calls the same command-line tools rather than reimplementing them. |
| `audit_symbol_matches.py` | Build-time audit: flags panel rows whose gene symbol collides with an unrelated enzyme. Run when editing the panel. |

### `panel_scored.py`

| flag | |
|---|---|
| `--tblout` | **required.** `hmmsearch --tblout` output. |
| `--map` | **required.** `ad_panel_map.tsv` from the build. |
| `--panel` | **required.** `assets/AD_methanogenesis_panel.tsv`. |
| `--out` | **required.** Per-module completeness TSV. |
| `--name` | Genome label for the report. Default `genome`. |
| `--gtdbtk` | GTDB-Tk `*.summary.tsv`. Repeatable (bac120 + ar53). Confirms the acetoclastic call and enables the genus hint. |
| `--dbcan` | Optional `run_dbcan` `overview.txt`. Only affects rows whose `scoring` column is `dbcan`; no panel row is currently set that way, since hydrolysis is detected with Pfam families. |
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

## Setup and maintenance

These live in `scripts/` alongside the panel builder. Despite the `.sbatch` extension on some of
them, none requires a scheduler: each falls back to sensible defaults when the SLURM
variables are absent, so they run directly with `bash` and can equally be submitted
with `sbatch`. They are site-neutral:

- **`DB` is required** and everything derives from it (`ncbifam/`, `ad_panel/`, `ad_panel_test/`). Nothing is defaulted to a particular filesystem.
- **No `--partition` or `--account`** in the sbatch headers. Pass them at submit time (`sbatch -p <partition> -A <account>`) or export `SBATCH_PARTITION` / `SBATCH_ACCOUNT`.
- **HMMER** is used from `PATH` if present (conda, container, spack); otherwise the module system is tried, with the module name overridable via `HMMER_MODULE`.

| script | run where | purpose |
|---|---|---|
| `scripts/prefetch_ncbifam.sh` | login | download NCBIfam + Rhea |
| `scripts/build_ad_panel.sbatch` | batch | build + press the panel |
| `scripts/fetch_test_proteomes.sh` | login | 7 reference proteomes for the smoke test |
| `scripts/smoke_test.sbatch` | batch | end-to-end run + assertions |
| `scripts/fetch_digester_mags.sh` | login | real digester MAGs (PRJEB31310) |
| `scripts/community_report.sbatch` | batch | score a MAG set + aggregate to a profile |

---

## Tests

`scripts/smoke_test.sbatch` runs the full path on seven reference proteomes and
asserts against published biology, so a failure means a real disagreement
rather than a drifted snapshot.

```bash
export DB=/path/to/databases
bash scripts/fetch_test_proteomes.sh
mkdir -p logs && sbatch -p <partition> scripts/smoke_test.sbatch
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
- **Non-curated (symbol- and EC-matched) rows are audited but weaker.** `audit_symbol_matches.py` flags collisions; gene-symbol matches are reliable, EC matches still warrant review before results are relied upon.
- **Some Pfam domains are shared across enzymes** (PF00871 = acetate *and* butyrate kinase; PF00374 = several [NiFe]-hydrogenases; PF06253 = the whole MttB superfamily). Gene-level calls on those need the NCBIfam gene-specific model, not the Pfam hit alone.
- **Coverage is honest about its gaps.** Modules and gates with no license-clean detector report `NA` / `NOT ASSESSABLE`, never `0.0`.

---

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

Only one panel row is currently non-default: `ACID-GLYC` is `assumed`, because
Embden-Meyerhof-Parnas is present in essentially every organism and so is not
diagnostic. The `dbcan` mode remains available for users who have run `run_dbcan`
themselves, but no panel row uses it: all three hydrolysis modules are detected with
Pfam catalytic families instead.

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

# 4. Hydrolysis: curated Pfam catalytic families, pinned per marker. dbCAN/CAZy is
#    deliberately not used; its public terms are non-commercial.
```

## Where this sits in an analysis

digestome consumes predicted proteins from metagenome-assembled genomes, together
with taxonomic assignments. It does not perform assembly, binning or gene calling;
any workflow that produces MAGs and proteins can feed it.

1. Assembly and binning produce MAGs.
2. Taxonomic classification (GTDB-Tk) produces lineages. These are not optional in
   practice: the acetoclastic call is confirmed against lineage, and without it that
   call falls back to gene evidence and over-calls.
3. Gene calling (for example Prodigal) produces one protein FASTA per MAG.
4. `hmmsearch` against the panel, then `panel_scored.py` per MAG: module
   completeness, branch gates, and secretion evidence for hydrolysis.
5. `aggregate_community.py` combines the scored genomes into one community profile.
