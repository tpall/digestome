# Methanogenesis markers and gate definitions

Reference for anyone reusing this panel's methanogenesis logic outside `digestome`,
including comparison against a KEGG-based marker set. Generated from
`assets/AD_methanogenesis_panel.tsv`; the gates are transcribed from
`bin/panel_scored.py`, which remains authoritative.

There is deliberately no KO mapping here. The panel exists to avoid KEGG, and
publishing a KO crosswalk would reintroduce the dependency it removes.

## Branch gates

Each gate is a conjunction of requirement groups; within a group the alternatives
are ORed. A gate evaluates to true, false, or **not assessable** when a required
group has no detectable member in the built database.

**Hydrogenotrophic** — `mcrA` AND at least 4 present among the C1 carriers:

    fwdB, fmdB, ftr, mch, mtd, hmd, mer, mtrA

Note the threshold counts **gene symbols, not pathway steps**: `fwdB`/`fmdB` are
alternatives for one enzyme, as are `mtd`/`hmd`, so the eight symbols cover six
steps. Four of eight is the coded condition.

**Acetoclastic** — `mcrA` AND `cdhA` AND (`acs` OR (`ackA` AND `pta`)), then
confirmed against GTDB-Tk lineage. Gene content alone cannot resolve this: the
ACDS/CODH complex is reversible and autotrophic hydrogenotrophs carry it to fix
carbon. Without taxonomy the gate falls back to gene evidence and over-calls.

**Methylotrophic** — `mcrA` AND at least one of:

    mtaB, mttB, mtbB, mtmB, mtsA, MEG-METHYL:comMT

`comMT` is the terminal MT2 step shared by every methylotrophic route. Its model is
subfamily-level and hits ~30 bacterial genomes alone, so it is sound **only** behind
the mandatory `mcrA` term. Do not lift it out of the gate.

## Why a lone C1 carrier cannot signal methanogenesis

Every gate carries a mandatory `mcrA` term, and `mcrA` is pinned exclusively to a
gene-specific model. A rise in a single shared carrier such as `mch` therefore
cannot produce a methanogenesis call: `mch` is one of eight symbols in a group
requiring four, behind a gate requiring `mcrA`. This matters when comparing against
marker sets that score genes independently, where a shared H4MPT step can move
without the pathway moving.

## A caveat on binned input

These markers are scored per genome, so prevalence measured over MAGs understates
markers carried by taxa that do not bin. A parallel analysis of bile-acid
7alpha-dehydroxylation found prevalence rising from 2.1% to 17.9% of samples when
whole-assembly proteomes were searched instead of binned ones, with 31% of carrying
contigs never binned. The trade-off is real in both directions: whole-assembly
proteomes give better prevalence, but cannot support the gates or the
taxonomy-confirmed acetoclastic call, since neither is meaningful for a mixed
community. Use binned input for capability and pathway completeness; use assembly
proteomes if the question is whether a marker is present in the sample at all.

## Pyrrolysine dependence

`mtmB`, `mtbB` and `mttB` all carry in-frame amber codons. They are recoverable from
Prodigal calls (standard table 11 truncates each into fragments that independently
clear the curated cutoffs) and absent from NCBI PGAP proteomes, which emit no protein
for them at all. Detection of these three depends on the upstream gene caller.

## Markers

### MEG-CORE

| marker_id | gene | enzyme | detector | tier |
|---|---|---|---|---|
| `MEG-CORE:mcrA` | mcrA | methyl-CoM reductase alpha | TIGR03256 | core |
| `MEG-CORE:mcrB` | mcrB | methyl-CoM reductase beta | TIGR03257 | core |
| `MEG-CORE:mcrG` | mcrG | methyl-CoM reductase gamma | TIGR03259 | core |
| `MEG-CORE:mcrC` | mcrC | McrC (assembly/activation) | TIGR03264 | accessory |
| `MEG-CORE:mcrD` | mcrD | McrD (assembly) | Pfam PF03296 | accessory |
| `MEG-CORE:hdrA` | hdrA | heterodisulfide reductase A | NF041778; NF041891 | core |
| `MEG-CORE:hdrB` | hdrB | heterodisulfide reductase B | Pfam PF02754 | core |
| `MEG-CORE:hdrC` | hdrC | heterodisulfide reductase C | Pfam PF13183;PF13534 | accessory |
| `MEG-CORE:mvhA` | mvhA | F420-non-reducing hydrogenase (large) | Pfam PF00374 | core |
| `MEG-CORE:mtrA` | mtrA | N5-methyl-H4MPT:CoM methyltransferase A | TIGR01111 | core |
| `MEG-CORE:mtrH` | mtrH | methyl-H4MPT:CoM methyltransferase H | Pfam PF01571 | core |

### MEG-HYDRO

| marker_id | gene | enzyme | detector | tier |
|---|---|---|---|---|
| `MEG-HYDRO:fwdB` | fwdB/fmdB | formylmethanofuran dehydrogenase (Mo/W subunit) | TIGR03129 | core |
| `MEG-HYDRO:fwdA` | fwdA/fmdA | formylmethanofuran dehydrogenase A | TIGR03121 | accessory |
| `MEG-HYDRO:ftr` | ftr | formyl-MFR:H4MPT formyltransferase | TIGR03119 | core |
| `MEG-HYDRO:mch` | mch | methenyl-H4MPT cyclohydrolase | TIGR03120 | core |
| `MEG-HYDRO:mtd` | mtd | F420-dep methylene-H4MPT dehydrogenase | NF002162 | core |
| `MEG-HYDRO:hmd` | hmd | H2-forming methylene-H4MPT dehydrogenase | Pfam PF03201 | accessory |
| `MEG-HYDRO:mer` | mer | F420-dep methylene-H4MPT reductase | TIGR03555 | core |
| `MEG-HYDRO:frhA` | frhA | F420-reducing [NiFe]-hydrogenase (large) | Pfam PF00374 | core |
| `MEG-HYDRO:frhB` | frhB | F420-reducing hydrogenase (F420-binding) | Pfam PF04432 | accessory |

### MEG-METHYL

| marker_id | gene | enzyme | detector | tier |
|---|---|---|---|---|
| `MEG-METHYL:mtaB` | mtaB | methanol:corrinoid methyltransferase | NF040651 | core |
| `MEG-METHYL:mtaC` | mtaC | methanol corrinoid protein | Pfam PF02310 | accessory |
| `MEG-METHYL:mtaA` | mtaA | methylcobalamin:CoM methyltransferase | Pfam PF03291 | accessory |
| `MEG-METHYL:mttB` | mttB | trimethylamine methyltransferase | TIGR02369 | core |
| `MEG-METHYL:mtbB` | mtbB | dimethylamine methyltransferase | TIGR02368 | accessory |
| `MEG-METHYL:mtmB` | mtmB | monomethylamine methyltransferase | PF05369 | accessory |
| `MEG-METHYL:mtsA` | mtsA | methylthiol:CoM methyltransferase | Pfam PF03291 | accessory |
| `MEG-METHYL:comMT` | MT2 | methylcobalamin:CoM methyltransferase (MT2 family) | NF004889 | core |

### MEG-ACETO

| marker_id | gene | enzyme | detector | tier |
|---|---|---|---|---|
| `MEG-ACETO:acs` | acs | acetyl-CoA synthetase (AMP-forming) | Pfam PF00501;PF13193 | core |
| `MEG-ACETO:ackA` | ackA | acetate kinase | Pfam PF00871 | core |
| `MEG-ACETO:pta` | pta | phosphotransacetylase | TIGR00651; NF004167; NF008852 | core |
| `MEG-ACETO:cdhA` | cdhA | CO dehydrogenase/acetyl-CoA decarbonylase, alpha | Pfam PF03063 | core |
| `MEG-ACETO:cdhC` | cdhC | ACDS beta (acetyl-CoA synthase) | Pfam PF03598 | core |
| `MEG-ACETO:cdhE` | cdhE/cdhD | ACDS corrinoid (delta/gamma) | Pfam PF00198;PF02607 | accessory |
