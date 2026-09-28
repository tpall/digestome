# bioRxiv submission — everything to paste, 2026-09-28

_Prepared for the account holder to submit at <https://submit.biorxiv.org>. **Submitted 2026-09-28**
as BIORXIV/2026/755039 (New Results, Bioinformatics, CC BY 4.0). Before approving: the author
block must show taavi@magrittr.ee (the account email was being switched from the UT address). The PDF is `manuscript/MANUSCRIPT.pdf` (16 pages with three figures, rendered from `MANUSCRIPT.qmd` at
the commit below in the quarto-preprint style: `quarto render MANUSCRIPT.qmd --to preprint-typst`,
which needs Quarto ≥ 1.9.36; `~/opt/quarto-1.10.18/bin/quarto` on the laptop)._

## Before you start

| | |
|---|---|
| File to upload | `~/Projects/digestome/manuscript/MANUSCRIPT.pdf` |
| Rendered from | `MANUSCRIPT.qmd`, digestome commit `359044c` |
| Software cited in it | digestome v0.1.2 (tag `v0.1.2`), Zenodo DOI 10.5281/zenodo.22913615 (all versions) |

## Field by field

**Type of submission:** New Results.

**Title**

> digestome: a licence-clean marker-gene panel for functional profiling of anaerobic digestion microbiomes

**Author**

| field | value |
|---|---|
| Name | Taavi Päll |
| ORCID | 0000-0001-7606-675X |
| E-mail | taavi@magrittr.ee |
| Corresponding | yes |
| Affiliation 1 | Magrittr OÜ, Tartu, Estonia |

**Abstract** (236 words, plain text — bioRxiv strips formatting, so paste from
`/tmp/abstract.txt` or the manuscript and check that the subscripts read as "H2/CO2" rather than
broken characters)

**Subject area:** *Bioinformatics*.
Alternative: *Microbiology*. Bioinformatics fits better — the contribution is a method and a
benchmark, and that is where a reader looking for AD profiling tools will browse. The category can
be changed on a later version, so it is not worth agonising over.

**Licence:** CC BY 4.0 recommended.
It matches the code (MIT) and the panel (CC BY 4.0), and lets anyone reuse the text with
attribution. CC BY-NC-ND would restrict commercial reuse of the *text* — it would not protect the
method, which is public anyway, and it looks defensive on a paper whose whole argument is
licence-freedom. CC0 gives away attribution, which we want to keep.

**Competing interest statement** (bioRxiv asks for this explicitly; it is also in the PDF)

> The author is the owner of Magrittr OÜ, which commercially provides anaerobic digestion microbiome
> analyses using the panel described herein. The panel, scoring methodology and validation scripts
> are openly available to facilitate independent replication and verification of the results.

**Funding statement**

> No external funding. Computation was bought from the University of Tartu High Performance
> Computing Centre under a commercial service agreement.

**Has this been submitted to a journal?** No.

**Data and code availability:** in the manuscript (GitHub + the Zenodo DOI above). No separate
supplementary files.

## What happens next

bioRxiv screens submissions before posting, normally within a day or two; the DOI appears when it
goes live. After that:

1. add the preprint DOI to the digestome README and to `CITATION.cff` as a `preferred-citation`;
2. add it to the Zenodo record as a related identifier (`isDocumentedBy`);
3. the vendor outreach can cite the preprint rather than only the software DOI — which is the reason
   this was moved forward (TP, 2026-09-24).

## One judgement call worth making consciously

Posting is public and permanent, and it puts the method in front of competitors as well as customers.
That was already true of the public repository and the Zenodo release; the preprint mainly makes it
*findable* by people searching the literature rather than GitHub. The argument for posting is that
the commercial offer rests on being the one who can show the method is correct — which requires
showing it.
