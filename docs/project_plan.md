# Host Range, Phylogeny and Acaricide Resistance in Phytophagous Mites — Project Plan

Nov 13, 2026 · V. K. Chaudhari

Copy of the project plan kept in the Claude Project "1. Mite Resistance And Phylogeny". If the two differ, update this file.

A comparative study testing whether host range and evolutionary relatedness predict which phytophagous mites evolve acaricide resistance, built entirely from public data.

## Highlights

- First phylogenetically corrected test of the pre-adaptation hypothesis in mites; earlier tests covered insects only.
- Sequence data are already sufficient: GenBank holds COI for about 183 named Tetranychidae, 46 Tenuipalpidae and 264 Eriophyidae species.
- Combines four public sources: resistance records, host-plant records, DNA sequences and life-history traits.
- Tests whether resistance is clustered on the mite tree, and whether broad host range predicts wider or faster resistance.
- Sets mites beside the published insect results.
- Tree searches and Bayesian runs on a Google Cloud spot machine at an estimated ₹2,000–3,000; no field or laboratory work.
- Ends with a ranked list of mite species at high risk of resistance, usable in Indian resistance management.

## Feasibility check: sequence availability

GenBank queried on 13 November 2026 for records of 150–3,000 bp, named species only.

| Family | Species with COI | With 18S | With 28S | With COI plus 18S or 28S | COI records |
| --- | --- | --- | --- | --- | --- |
| Tetranychidae | 183 | 103 | 123 | 116 | 3,597 |
| Tenuipalpidae | 46 | 4 | 25 | 23 | 971 |
| Eriophyidae | 264 | 167 | 256 | 149 | 3,606 |

Main pests well covered in COI: *Tetranychus urticae* (630 records), *Brevipalpus yothersi* (339), *Panonychus ulmi* (217), *Tetranychus evansi* (168), *Panonychus citri* (104), *Aculops lycopersici* (103), *Tetranychus truncatus* (99). Counts will shift after name standardisation. BOLD records and the count of species with resistance cases are added in stage 1.

## Background and research gap

- Pre-adaptation hypothesis: herbivores meet synthetic pesticides with detoxification systems evolved against host-plant chemicals. Supported in 615 North American arthropod pests, without separating feeding habit from taxonomy (Rosenheim et al.).
- Hardy et al. (2018), Evolutionary Applications: phylogenetically corrected test, US insect pests only.
- Crossley et al. (2021), Evolutionary Applications 14: 290–296: >17,000 herbivore–insecticide combinations; faster resistance in insects with broad, but not the broadest, diets.
- Jin et al. (2019), BMC Evolutionary Biology 19: phylogenetic signal in spider mite abundance and range, positively related to host range; resistance not analysed.

**Gap.** No phylogenetically corrected test of whether host range predicts acaricide resistance across phytophagous mites.

**Backbone.** RNA-Seq phylogeny of 72 Tetranychinae species (Matsuda et al., 2018, PLoS ONE 13: e0203136). COI alone does not resolve deeper Tetranychidae nodes (Matsuda et al., 2014, PLoS ONE 9: e108672).

## Objectives

1. Compile a curated global dataset of acaricide resistance cases in Tetranychidae, Tenuipalpidae and Eriophyidae.
2. Compile host range for the same species (host plants and host plant families).
3. Build a species-level phylogeny from public sequences.
4. Measure phylogenetic signal in resistance traits.
5. Test whether host range, voltinism and pest status predict resistance after phylogenetic correction.
6. Compare mite results with published insect results and rank mite species by resistance risk.

Hypotheses H1–H4 and variables: see `protocol.md`.

## Data sources

| Data | Source | Access |
| --- | --- | --- |
| Resistance cases | Arthropod Pesticide Resistance Database (pesticideresistance.org), plus papers | Website search open; full export by request to Dr David Mota-Sanchez |
| Host plants | Spider Mites Web (Migeon and Dorkeld), plus published records | Open |
| Sequences | GenBank and BOLD: COI, ITS2, 28S, mitogenomes | Open |
| Voltinism, body size, pest status | Literature | Open |
| Acaricide launch years | IRAC classification and registration records | Open |
| Insect comparison | Hardy et al. (2018), Crossley et al. (2021) supplements | Open |

Checkpoint at end of stage 1: about 30 species with at least one resistance case, inside a tree of 80 or more species.

## Stages

### Stage 1. Feasibility and protocol
- Export Acari records from the Arthropod Pesticide Resistance Database, family by family.
- Match names against GenBank and BOLD; species-by-marker availability table.
- Write and date `protocol.md` before modelling.
- Output: `feasibility_counts.csv`, `protocol.md`. Check: go/no-go on species count.

### Stage 2. Resistance dataset
- Database export plus documented search of Scopus, CAB Abstracts and Google Scholar (strings and dates logged).
- Case: one species, one active ingredient, one country; bioassay RR ≥ 10 vs susceptible reference, or field failure confirmed by bioassay.
- Fields: species as published, standardised name, active ingredient, IRAC group, country, crop, collection year, publication year, RR, bioassay method, reference with DOI.
- First registration year of each active ingredient from IRAC documents and pesticide compendia.
- Output: `resistance_cases.csv`, `resistance_by_species.csv`, record-flow diagram. Check: 10% of rows re-checked against sources.

### Stage 3. Host range dataset
- Spider Mites Web for Tetranychidae; catalogues for Tenuipalpidae and Eriophyidae.
- Plant names standardised against World Flora Online; laboratory-only records excluded.
- Metrics: host species, genera, families; host phylogenetic diversity; rarefied breadth.
- Output: `host_records.csv`, `host_breadth.csv`. Check: 20 broadest and 20 narrowest species.

### Stage 4. Sequences and alignment
- COI, 18S, 28S from GenBank and BOLD with voucher, country, host.
- Filters: length; COI stop codons (invertebrate mitochondrial code); genus-nesting misidentification flag; duplicates.
- One representative per species and marker, same voucher across markers where possible.
- MAFFT (codon-aware for COI); trimAl for ribosomal genes; partition by gene and COI codon position; Trombidiformes outgroups.
- Output: `concat.fasta`, `partitions.nex`, `accessions.csv`, excluded-record log.

### Stage 5. Phylogeny
- Separate Tetranychoidea and Eriophyoidea trees joined at a fixed backbone.
- IQ-TREE: ModelFinder partitions, ten searches, 1,000 UFBoot, SH-aLRT; Tetranychinae constrained to Matsuda et al. (2018).
- BEAST 2 relaxed clock, relative dating (root = 1), four chains, ESS > 200; 1,000 posterior trees kept.
- Output: `mite_ml.treefile`, `mite_dated.tre`, `posterior_1000.trees`.

### Stage 6. Comparative analyses
- H1: D statistic (caper); Pagel's λ and Blomberg's K (phytools), 1,000 permutations.
- H2: MCMCglmm or brms, Poisson / negative binomial, quadratic host breadth.
- H3: coxme with species and active ingredient random effects and phylogenetic covariance; right-censoring.
- Research effort and geographic range in every model; run over 100 posterior trees and pooled.
- Sensitivity: without *T. urticae*; Tetranychidae only; taxonomy-only tree; RR ≥ 100; rarefied host breadth; *T. cinnabarinus* merged/separate.

### Stage 7. Comparison with insects
- Refit H2 and H3 structures on Hardy et al. (2018) and Crossley et al. (2021) data; joint model with group × host-breadth interaction.

### Stage 8. Risk ranking and validation
- Predicted risk for every species; temporal hold-out (fit to 2010, test after 2010); leave-one-species-out; Indian crop species extracted.

### Stage 9. Manuscript
- Fig. 1 tree with traits; Fig. 2 H2 effect sizes; Fig. 3 time-to-resistance curves; Fig. 4 mites vs insects; Fig. 5 risk ranking.
- Supplement: dataset, record-flow diagram, model tables, accessions, code archive with DOI.

## Cloud computing

c2d-standard-16 spot machine in asia-south1, 50 GB disk, one Cloud Storage bucket, hourly sync. Estimated 80–130 machine hours; ₹1,100–1,900 on spot. Budget alert ₹8,000. One-hour pilot before full runs; machine stopped when idle and deleted after trees are archived. Free Trial (US$300, 90 days) started only when the alignment is ready.

## Target journals

Evolutionary Applications; Journal of Pest Science; Pest Management Science; Experimental and Applied Acarology.

## Risks

| Risk | Fallback |
| --- | --- |
| Fewer than 30 resistant species | Keep zeros; mite–insect comparison leads |
| Database export refused | Compile from website search and papers |
| Reporting bias | Research effort in every model; run without *T. urticae* |
| Deep Tetranychoidea–Eriophyoidea divergence | Separate trees; Tetranychidae-only analysis |
| COI saturation | 18S, 28S and phylogenomic backbone |
| Spot machine reclaimed | Checkpoints, hourly sync, restart script |
| *T. urticae* red/green form naming | Fixed naming standard; both treatments as sensitivity |
