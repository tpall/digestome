# Command reference

Flags and scripts. For what the numbers mean, see [scoring.md](scoring.md); for how markers were chosen, [panel.md](panel.md).

## The tools

| script | what it does |
|---|---|
| `scripts/prefetch_ncbifam.sh` | Downloads NCBIfam HMMs + `rhea2ec.tsv`. Login node (compute nodes may lack outbound internet). Idempotent and resumable. |
| `scripts/build_ad_hmm_db.sh` | Matches the panel against NCBIfam, emits `ad_panel.hmm` (pressed), `ad_panel_map.tsv` (model → gene/module) and `rhea2ec.tsv`. Audits curated accessions and asserts every pressed model has an NC cutoff. |
| `bin/panel_scored.py` | Per genome: `hmmsearch --tblout` (+ optional dbCAN, GTDB-Tk) → per-module completeness table on `--out`, gates + diagnostic markers on stdout. |
| `bin/aggregate_community.py` | Many scored genomes → one community profile. Route balance, acetate-consumer presence, syntrophic capacity. JSON + text. |
| `bin/compare_profiles.py` | Two or more profiles of the same reactor in time order → what moved: routes, modules and genomes over time with a coarse change call (≥ 1 pp and ≥ 1.5-fold by default). TSV + JSON + text. |
| `main.nf`, `nextflow.config` | Nextflow workflow over the analysis path: HMMSEARCH, SCORE, AGGREGATE. Calls the same command-line tools rather than reimplementing them. |
| `scripts/audit_symbol_matches.py` | Build-time audit: flags panel rows whose gene symbol collides with an unrelated enzyme. Run when editing the panel. |

### `panel_scored.py`

| flag | |
|---|---|
| `--tblout` | **required.** `hmmsearch --tblout` output. |
| `--map` | **required.** `ad_panel_map.tsv` from the build. |
| `--panel` | **required.** `assets/AD_methanogenesis_panel.tsv`. |
| `--out` | **required.** Per-module completeness TSV. |
| `--name` | Genome label for the report. Default `genome`. |
| `--gtdbtk` | GTDB-Tk `*.summary.tsv`. Repeatable (bac120 + ar53). Confirms the acetoclastic call and enables the genus hint. |
| `--secretion` | Extracellular-targeting modules (`assets/secretion_modules.tsv`). Adds the `core_secreted` column: a hydrolysis marker counts as exported only when the same protein also carries a dockerin, cohesin, CBM or SLH domain. Omit it and the column reports `NA`, not zero. |
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
| `--out-tables DIR` | Result tables, one row per genome / module / route / risk (below). |
| `--table-format` | `tsv` (default) or `csv` (RFC 4180, quoted). TSV is the repository convention and survives the commas and semicolons in lineages and notes; use `csv` when a client's tooling insists. |

Stdlib only, so it runs anywhere Python 3 does.

### `compare_profiles.py`

    python3 bin/compare_profiles.py t1/profile.json t2/profile.json t3/profile.json \
        --labels 2026-03 2026-06 2026-09 --reactor "Reactor 1" --out monitoring/

| flag | |
|---|---|
| profiles | **required.** `profile.json` files in sampling order, oldest first. |
| `--labels` | Column labels in the same order, e.g. dates. Default: the sample names. |
| `--reactor` | Name used in the summary. |
| `--out DIR` | **required.** Writes `routes_over_time.tsv`, `modules_over_time.tsv`, `genomes_over_time.tsv`, `changes.json`, `summary.txt`. |
| `--min-delta` | Percentage points a share must move to count as a change (default 1.0). |
| `--min-fold` | Fold change it must also show (default 1.5). |
| `--floor` | Genomes below this share in every sample are left out of the genome table (default 0.5 %). |

With two or three samples there is no trend to test, only a difference. The change call is deliberately coarse so that 0.3 % → 0.4 % is never reported as movement and 0.1 % → 3 % always is; modules with no scored genome stay `not assessable`. Directions are `new`, `lost`, `up`, `down`, `stable`, `absent`.

#### Result tables (`--out-tables`)

Flat files derived from the same profile as the JSON, so a report built on them cannot disagree with it. Every table starts with a `sample` column so files from repeated sampling of one reactor concatenate. Unmeasured is written as `NA`, never `0`.

| file | one row per | columns |
|---|---|---|
| `genomes` | MAG, sorted by abundance | `abundance_pct`, lineage/phylum/genus, `mcrA`, one yes/no/NA column per gated route and capability, `gate_notes` (why a call went the way it did, as `gate:code` — `taxonomy-confirmed`, `absent-on-taxonomy` for a genome that carries the markers but is not a known acetoclastic lineage, `markers-only` when no taxonomy was supplied; the full sentence stays in `profile.json`), `<MODULE>_pct_complete` for every panel module, `HYDROL-*_families_found` / `_families_exported` (the secretion test), `marker_<gene>` for the diagnostic markers |
| `modules` | panel module | genomes scored vs unscoreable and the `status` that explains the NA, carriers at any / ≥ 50 % / 100 % completeness, share of reads carried by the ≥ 50 % carriers, best genome and its completeness, and for hydrolysis modules the number and share of genomes that **export** the enzyme |
| `routes` | methanogenesis route or capability | present / total / not assessable, share of reads, the genomes with genus |
| `risks` | risk statement | level (`ok` / `note` / `attention`), title, detail |

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
| `scripts/prefetch_pfam.sh` | login | pinned Pfam families, incl. the secretion modules |
| `scripts/fetch_biogas_catalogue.sh` | login | the Campanaro biogas MAG catalogue, for benchmarking |
| `scripts/score_catalogue.sbatch` | batch | score a whole catalogue; the specificity benchmark |
| `scripts/subset_catalogue_plant.py` | login (first run) | cut one plant's community out of the catalogue by per-sample abundance (Additional file 8); output scores with `score_catalogue.sbatch` |
