#!/usr/bin/env bash
# Stage 5c: maximum-likelihood tree of Eriophyoidea (11 taxa), no constraint.
# Partition scheme and models from ModelFinder in the pilot run (05a); ten independent searches,
# each with 1,000 ultrafast bootstraps and 1,000 SH-aLRT replicates; outgroups Phytoptus and Nalepella.
set -euo pipefail
A=data/clean/alignments
mkdir -p work/stage5 results/trees
iqtree2 -s $A/concat_eriophyoidea.fasta -p results/trees/pilot/pilot_eriophyoidea.best_scheme.nex \
  --runs 10 -B 1000 -alrt 1000 -T 2 --seed 2026 --prefix work/stage5/final_erio \
  -o Phytoptus_avellanae,Nalepella_tsugifoliae "$@"
cp work/stage5/final_erio.treefile results/trees/eriophyoidea_ml.treefile
cp work/stage5/final_erio.contree  results/trees/eriophyoidea.contree
cp work/stage5/final_erio.iqtree   results/trees/eriophyoidea_ml.iqtree
cp work/stage5/final_erio.runtrees results/trees/eriophyoidea.runtrees
cp work/stage5/final_erio.log      results/trees/eriophyoidea.log
