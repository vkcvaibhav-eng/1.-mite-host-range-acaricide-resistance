# Protocol

Status: approved by V. K. Chaudhari on 13 November 2026 at the stage 1 checkpoint. Fixed from this date; later changes go in the change log.

## Hypotheses

| Code | Hypothesis | Response variable | Test |
| --- | --- | --- | --- |
| H1 | Resistance is clustered on the mite phylogeny | Resistant or not; number of IRAC groups resisted | Fritz–Purvis D; Pagel's λ; Blomberg's K |
| H2 | Broader host range predicts more IRAC groups resisted | Number of IRAC groups with at least one case | Phylogenetic mixed model for counts, quadratic host-breadth term |
| H3 | Broader host range shortens time to first resistance | Years from first registration to first case, per species and active ingredient | Mixed-effects Cox model with phylogenetic covariance |
| H4 | Mite and insect host-range effects run in the same direction | Standardised effect sizes | Same models on published insect data; joint model with group × host-range interaction |

## Species set

All phytophagous mites (Tetranychidae, Tenuipalpidae, Eriophyidae) recorded on cultivated crops with usable sequence data, including species with no resistance record (entered as zeros; research effort modelled).

## Resistance case definition

Main analysis: a case as recorded in the Arthropod Pesticide Resistance Database (one species, one active ingredient, one location and report), excluding populations resistant only after laboratory selection and duplicate database entries of the same paper. Cases added from the literature for species absent from the database follow the same rule and, where reported, must show a resistance ratio of 10 or more against a susceptible reference or a field control failure confirmed in a bioassay.

Sensitivity runs: (a) field-evolved cases only; (b) cases supported by a journal paper, excluding FAO reports and personal communications.

## Variables

| Variable | Definition | Role |
| --- | --- | --- |
| Resistant | At least one accepted case | Response (H1) |
| IRAC groups resisted | Number of IRAC groups with at least one case | Response (H1, H2) |
| Time to resistance | Years from first registration to first case | Response (H3) |
| Host breadth | Number of host plant families; host phylogenetic diversity; rarefied breadth | Main predictor |
| Voltinism | Generations per year | Covariate |
| Reproductive mode | Arrhenotoky or thelytoky | Covariate |
| Geographic range | Number of countries with records | Covariate |
| Pest status | Number of crops listed as host in pest literature | Covariate |
| Research effort | Number of publications naming the species | Bias control |

## Exclusion rules

- Host records from laboratory rearing only.
- COI sequences with stop codons under the invertebrate mitochondrial code.
- Sequences nesting inside another genus (flagged as probable misidentification).
- Unidentified or genus-only records.

## Models

| Hypothesis | Model | Software | Details |
| --- | --- | --- | --- |
| H1 | Fritz–Purvis D (binary); Pagel's λ and Blomberg's K (count) | `caper`, `phytools` | 1,000 permutations each |
| H2 | Bayesian phylogenetic mixed model, Poisson or negative binomial | `MCMCglmm` or `brms` | Predictors scaled; host breadth with a quadratic term |
| H3 | Mixed-effects Cox model | `coxme` | Species and active ingredient as random effects; phylogenetic covariance; pairs with no case by the cut-off year right-censored |
| H4 | H2 and H3 structures refitted on published insect data; joint model with group × host-breadth interaction | as above | Standardised coefficients compared with intervals |

Every model includes research effort and geographic range, and is run over 100 posterior trees with estimates pooled across trees.

## Sensitivity runs

- Without *Tetranychus urticae*
- Tetranychidae only
- Taxonomy-only tree
- Field-evolved cases only
- Journal-supported cases only (no FAO reports or personal communications)
- Rarefied host breadth
- *Tetranychus cinnabarinus* merged with and separate from *T. urticae*

## Validation

- Temporal hold-out: fit on cases reported up to 2010; test on species whose first case came after 2010
- Leave-one-species-out cross-validation

## Name standard

Stage 1 decisions are in `data/manual/name_map_stage1.csv`. The full naming standard (Spider Mites Web for Tetranychidae; current catalogues for Tenuipalpidae and Eriophyidae) is fixed at stage 2.

## Change log

| Date | Change | Reason |
| --- | --- | --- |
| 13 November 2026 | Case definition changed from a resistance-ratio rule to database cases as recorded, with two new sensitivity runs (field-evolved only; journal-supported only) | Database records do not report resistance ratios, and many early sources (FAO reports, personal communications) cannot be checked; approved by V. K. Chaudhari |
