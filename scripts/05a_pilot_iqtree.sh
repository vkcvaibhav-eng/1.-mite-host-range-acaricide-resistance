#!/usr/bin/env bash
# Stage 5a: pilot maximum-likelihood trees, no topological constraint.
# IQ-TREE 2.3.6; ModelFinder with partition merging; 1,000 ultrafast bootstraps; 1,000 SH-aLRT; seed 2026.
# Runs on 2 cores in about 15 minutes (Tetranychoidea) and under 1 minute (Eriophyoidea).
set -euo pipefail
A=data/clean/alignments
mkdir -p work/stage5 results/trees/pilot
iqtree2 -s $A/concat_eriophyoidea.fasta -p $A/partitions_eriophyoidea.nex -m MFP+MERGE -B 1000 -alrt 1000 \
  -T 2 --seed 2026 --prefix work/stage5/pilot_erio -o Phytoptus_avellanae,Nalepella_tsugifoliae
# no -o here: IQ-TREE 2.3.6 aborts in setRootNode when two outgroups are given for this matrix; the tree is
# rooted afterwards on Leptotrombidium deliense + Cheyletus malaccensis
iqtree2 -s $A/concat_tetranychoidea.fasta -p $A/partitions_tetranychoidea.nex -m MFP+MERGE -B 1000 -alrt 1000 \
  -T 2 --seed 2026 --prefix work/stage5/pilot_tetra
cp work/stage5/pilot_tetra.contree results/trees/pilot/pilot_tetranychoidea.contree
cp work/stage5/pilot_tetra.iqtree  results/trees/pilot/pilot_tetranychoidea.iqtree
cp work/stage5/pilot_erio.contree  results/trees/pilot/pilot_eriophyoidea.contree
cp work/stage5/pilot_erio.iqtree   results/trees/pilot/pilot_eriophyoidea.iqtree
cp work/stage5/pilot_tetra.treefile results/trees/pilot/pilot_tetranychoidea.treefile
cp work/stage5/pilot_erio.treefile  results/trees/pilot/pilot_eriophyoidea.treefile
cp work/stage5/pilot_tetra.best_scheme.nex results/trees/pilot/pilot_tetranychoidea.best_scheme.nex
cp work/stage5/pilot_erio.best_scheme.nex  results/trees/pilot/pilot_eriophyoidea.best_scheme.nex
cp work/stage5/pilot_tetra.best_model.nex results/trees/pilot/pilot_tetranychoidea.best_model.nex
cp work/stage5/pilot_erio.best_model.nex  results/trees/pilot/pilot_eriophyoidea.best_model.nex
