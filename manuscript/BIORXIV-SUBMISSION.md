# bioRxiv submission — everything to paste, 2026-09-24

_Prepared for the account holder to submit at <https://submit.biorxiv.org>. Nothing here has been
submitted. The PDF is `manuscript/MANUSCRIPT.pdf` (13 pages, rendered from `MANUSCRIPT.qmd` at the
commit below)._

## Before you start

| | |
|---|---|
| File to upload | `~/Projects/digestome/manuscript/MANUSCRIPT.pdf` |
| Rendered from | `MANUSCRIPT.qmd`, digestome commit `f904c48` |
| Software cited in it | Zenodo v0.1.0 DOI 10.5281/zenodo.22913616 (already live) |

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
| Affiliation 2 | University of Tartu, Tartu, Estonia |

**Abstract** (219 words, plain text — bioRxiv strips formatting, so paste from
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

> The author is the owner of Magrittr OÜ, which sells analyses of anaerobic digestion microbiomes to
> biogas plant operators using the panel described here. The panel, the scorer and the validation
> scripts are released openly so that the results can be reproduced and checked independently of
> that commercial interest.

**Funding statement**

> No external funding. The work was carried out by the author, and the compute for the benchmarks ran
> on the University of Tartu HPC centre, in part on a paid commercial allocation held by Magrittr OÜ.

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
