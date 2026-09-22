# Module plan: electron acceptors (sulfate and nitrate reduction)

_Status: planned, not implemented (2026-09-22). Owner: TP._

## Why

The panel sees methanogens and the fermentation chain, but not the organisms that compete with them
for the same substrates. For a digester operator the story, in order of how often it matters:

1. **Competition.** At high sulfate (distillery and other industrial wastewaters, seaweed, some
   co-substrates) sulfate reducers outcompete methanogens for H₂ and acetate and lower the gas yield.
2. **H₂S in the biogas.** Sulfate reducers make it; operators pay for desulfurisation.
3. **Methane loss.** Anaerobic methane oxidisers (now reported apart from methanogens, role
   `methane_oxidiser`) need an acceptor: sulfate via partner bacteria (ANME-1, ANME-2a/b/c), or
   nitrate, iron or manganese (ANME-2d, *Methanoperedens*). Without the acceptor side the report can
   name them but not say whether they can be active. Rare in digesters, but costly where present.

## Scope

| process | core markers (one condition per step) | notes |
|---|---|---|
| dissimilatory sulfate reduction | `sat`, `aprA`\|`aprB`, `dsrA`, `dsrB`, `dsrD` | `dsrD` separates reductive from oxidative DsrAB |
| sulfur oxidation via reverse Dsr (negative control) | `dsrA`, `dsrB`, `dsrE`/`dsrF`/`dsrH` | same DsrAB, other direction |
| dissimilatory nitrate reduction | `narG` \| `napA` | respiratory, not assimilatory (`nasA` excluded) |
| nitrite to ammonium (DNRA) | `nrfA` | adds to NH₄⁺ load; relevant to the ammonia rules |
| Fe(III) / Mn(IV) reduction | none | multiheme cytochromes have no reliable single marker: reported as **not assessed** |

Candidate models (NCBIfam / TIGRFAM, public domain), **to verify and pin by accession** as the other
markers were: TIGR00339 (`sat`), TIGR02061 / TIGR02060 (`aprA` / `aprB`), TIGR02064 / TIGR02066
(`dsrA` / `dsrB`), Pfam PF08679 (`dsrD`), TIGR01580 (`narG`), TIGR01706 (`napA`), TIGR03152 (`nrfA`).
Check each against the NCBIfam release in `assets/` before use; pin to the gene-specific model, not
a superfamily.

## The direction problem, again

DsrAB and AprAB run both ways: sulfate reducers reduce, many sulfur oxidisers (e.g. *Chlorobiaceae*,
some *Proteobacteria*) run them in reverse. Same pattern as ACDS/CODH and the C1 pathway, same
answer: a gene gate plus a discriminating marker (`dsrD` for reduction, `dsrEFH` for oxidation) and,
where that is not enough, a lineage rule in `assets/acetate_lineages.tsv` (new roles, e.g.
`sulfate_reducer` allow list at family level for *Desulfobacterota*, `sulfur_oxidiser` deny list).
Fail closed: an unresolved direction is reported as "sulfur-cycle genes, direction not resolved",
not as a sulfate reducer.

## Report output (renderer)

- Result row "Sulfate reducers" and "Nitrate reducers", count and share; status note when present.
- Findings: sulfate reducers + sulfate in feed → H₂S and competition; methane oxidisers + a matching
  acceptor organism → "methane oxidation possible" (ANME-1/2abc with sulfate reducers; ANME-2d with
  nitrate); methane oxidisers without one → "acceptor not detected in the genomes; check the feed".
- Suggestion rules (in `report/suggestion_rules.tsv`): measure H₂S in the gas and sulfate in the feed
  when sulfate reducers are present; `process_sulfate_mg_per_l`, `process_nitrate_mg_per_l` (already
  declared) as context; no numeric setpoints.
- Literature per statement, as for the genus notes.

## Validation before release

1. Smoke test: add *Desulfovibrio vulgaris* (sulfate reduction yes), a green sulfur bacterium
   (reverse Dsr: sulfate reduction **no**), *E. coli* (nitrate reduction yes, sulfate reduction no).
2. GTDB r226 sample: sulfate-reduction calls by lineage; expect *Desulfobacterota* and a few
   Firmicutes (*Desulfotomaculum*); every call outside known sulfate-reducer lineages is reviewed.
3. Campanaro catalogue: counts per plant; no call without `dsrD`.
4. Manuscript: one paragraph and a row in the specificity table, or a follow-up note.

## Effort

About three working days: models and pinning 1, gates and lineage roles 0.5, benchmarks and smoke
test 1 (HPC), renderer rules and tests 0.5.
