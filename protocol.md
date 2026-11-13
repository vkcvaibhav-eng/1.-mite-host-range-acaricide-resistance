# Protocol

Status: draft. To be completed and committed with a date at the end of stage 1, before any modelling.

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

One species, one active ingredient, one country, supported by a bioassay with a resistance ratio of 10 or more against a susceptible reference, or by a field control failure confirmed in a bioassay.

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

## Models and sensitivity runs

To be completed at end of stage 1. Sensitivity runs planned: without *Tetranychus urticae*; Tetranychidae only; taxonomy-only tree; resistance ratio ≥ 100; rarefied host breadth; *T. cinnabarinus* merged and separate.

## Change log

| Date | Change | Reason |
| --- | --- | --- |
