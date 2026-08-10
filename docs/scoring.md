# How a module is scored

What `pct_complete`, `NA` and each branch gate actually mean. Command flags are in [cli.md](cli.md).

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
