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
docs/        the reference material: flags, scoring rules, panel curation.
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

Writes `results/profile.txt`, `results/profile.json` and the result tables in
`results/tables/` (genomes, modules, routes, risks; see `docs/cli.md`), plus a Nextflow trace,
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
    --out-json profile.json --out-txt profile.txt --out-tables tables/
```

### Which to use

Option B for a few genomes, for debugging, or where installing a workflow engine
is not worth it. Option A once there are enough genomes that parallelism and
resume matter, or when the same analysis has to move between a cluster and a
single machine. They produce identical results because they run the same code.

Use `--cut_nc`: every model in the panel carries a curated NCBIfam noise cutoff,
and the build asserts this. Do not substitute a global `-T`; the curated cutoffs
across this panel span roughly 20 to 1400 bits.

**`--gtdbtk` is not optional in practice.** The acetoclastic gate is confirmed
against GTDB-Tk lineage, because ACDS/CODH is reversible and gene content alone
cannot distinguish acetyl-CoA cleavage (acetoclastic methanogenesis) from
acetyl-CoA synthesis (autotrophic carbon fixation). Without it, obligate
hydrogenotrophs carrying ACDS are called acetoclastic. See the
`Mthermautotroph` case in `tests/expectations.tsv`.

---

## Documentation

The detail lives in `docs/`, so this page stays the getting-started page:

| | |
|---|---|
| [docs/cli.md](docs/cli.md) | Every flag of both tools, and what each script in `scripts/` is for. |
| [docs/scoring.md](docs/scoring.md) | How module completeness is computed, what `NA` means, and the condition behind each branch gate. |
| [docs/panel.md](docs/panel.md) | Which database detects which row, the single-marker fingerprint, and how to resolve a new accession. |
| [docs/methanogenesis-markers.md](docs/methanogenesis-markers.md) | Branch gate definitions and methanogenesis marker accessions, for reuse or comparison against another marker set. |

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
