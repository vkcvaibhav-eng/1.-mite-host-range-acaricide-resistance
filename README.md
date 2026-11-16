# Host range, phylogeny and acaricide resistance in phytophagous mites

**Working title:** Do host range and phylogeny predict acaricide resistance in phytophagous mites? A comparative test of the pre-adaptation hypothesis

V. K. Chaudhari, AINP on Agricultural Acarology, Department of Entomology, N. M. College of Agriculture, Navsari Agricultural University, Navsari, Gujarat, India

## Aim

A comparative study, built only from public data, that tests whether host range and evolutionary relatedness predict which phytophagous mites (Tetranychidae, Tenuipalpidae, Eriophyidae) evolve acaricide resistance.

| Code | Hypothesis |
| --- | --- |
| H1 | Resistance is clustered on the mite phylogeny |
| H2 | Species with broader host range resist more IRAC mode-of-action groups |
| H3 | Broader host range shortens the time from acaricide launch to first resistance report |
| H4 | The host-range effect in mites runs in the same direction as in insects |

The full plan is in [`docs/project_plan.md`](docs/project_plan.md). Hypotheses, variables, models and exclusion rules are fixed in [`protocol.md`](protocol.md) before any modelling.

## Data sources

- Resistance cases: Arthropod Pesticide Resistance Database and published papers
- Host plants: Spider Mites Web and published host catalogues
- Sequences: GenBank and BOLD (COI, 18S, 28S)
- Life-history traits: literature
- Insect comparison: supplements of Hardy et al. (2018) and Crossley et al. (2021)

## Layout

```
protocol.md            hypotheses, variables, models, exclusion rules
environment.yml        conda programs (IQ-TREE, MAFFT, trimAl, BEAST 2, R)
Makefile               one command rebuilds tables and figures
data/raw/              downloads, never edited
data/manual/           hand decisions read by scripts (name maps)
data/clean/            curated tables, alignments, trees
scripts/               numbered scripts, one per stage
cloud/                 Google Cloud spot machine scripts
results/tables/        every number used in the manuscript
results/figures/
logs/                  decisions and excluded records
docs/                  project plan and notes
```

## Stages

1. Feasibility and protocol
2. Resistance dataset
3. Host range dataset
4. Sequences and alignment
5. Phylogeny
6. Comparative analyses
7. Comparison with insects
8. Risk ranking and validation
9. Manuscript

## Rules

- Files in `data/raw/` are never edited.
- Every excluded record is logged in `logs/` with its reason.
- Every number in the manuscript comes from a committed script writing to `results/tables/`.

## Reproduce

```bash
conda env create -f environment.yml
conda activate mite-resistance
make all
```

## Status

Stage 1 (feasibility) run on 13 November 2026 with the full resistance-database species list and BOLD counts; approved (go) on 13 November 2026. Stage 2 (resistance dataset) approved 16 November 2026: `data/clean/resistance_cases.csv`, `data/clean/resistance_by_species.csv`, introduction years in `data/clean/ai_introduction_years.csv`; record flow in `results/tables/stage2_record_flow.csv`; paper-level verification in `logs/verification_stage2.md` and `results/tables/stage2_verification_summary.csv`; corrections in `data/manual/case_corrections_stage2.csv`. Results: `results/tables/stage1_summary.csv`, `data/clean/feasibility_counts.csv`; decisions: `logs/decisions.md`.
