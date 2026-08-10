# Command reference

Flags and scripts. For what the numbers mean, see [scoring.md](scoring.md); for how markers were chosen, [panel.md](panel.md).

## The tools

| script | what it does |
|---|---|
| `scripts/prefetch_ncbifam.sh` | Downloads NCBIfam HMMs + `rhea2ec.tsv`. Login node (compute nodes may lack outbound internet). Idempotent and resumable. |
| `scripts/build_ad_hmm_db.sh` | Matches the panel against NCBIfam, emits `ad_panel.hmm` (pressed), `ad_panel_map.tsv` (model → gene/module) and `rhea2ec.tsv`. Audits curated accessions and asserts every pressed model has an NC cutoff. |
| `bin/panel_scored.py` | Per genome: `hmmsearch --tblout` (+ optional dbCAN, GTDB-Tk) → per-module completeness table on `--out`, gates + diagnostic markers on stdout. |
| `bin/aggregate_community.py` | Many scored genomes → one community profile. Route balance, acetate-consumer presence, syntrophic capacity. JSON + text. |
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
| `scripts/prefetch_pfam.sh` | login | pinned Pfam families, incl. the secretion modules |
| `scripts/fetch_biogas_catalogue.sh` | login | the Campanaro biogas MAG catalogue, for benchmarking |
| `scripts/score_catalogue.sbatch` | batch | score a whole catalogue; the specificity benchmark |
