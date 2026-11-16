#!/usr/bin/env bash
# Stage 5b: maximum-likelihood tree of Tetranychoidea with the Tetranychinae backbone of Matsuda et al. (2018).
# The constraint (data/manual/constraint_matsuda2018.nwk) is pruned to the taxa of the alignment first.
# Partition scheme and models are those chosen by ModelFinder in the pilot run (05a).
# Ten independent searches, each with 1,000 ultrafast bootstraps and 1,000 SH-aLRT replicates; IQ-TREE keeps the best.
set -euo pipefail
A=data/clean/alignments
mkdir -p work/stage5 results/trees
python3 -W ignore - <<'PY'
import dendropy
t = dendropy.Tree.get(path="data/manual/constraint_matsuda2018.nwk", schema="newick", preserve_underscores=True)
aln = {l[1:].strip() for l in open("data/clean/alignments/concat_tetranychoidea.fasta") if l.startswith(">")}
t.retain_taxa_with_labels([l.taxon.label for l in t.leaf_node_iter() if l.taxon.label in aln])
t.suppress_unifurcations()
open("work/stage5/constraint_pruned.nwk", "w").write(
    t.as_string(schema="newick", suppress_rooting=True, suppress_edge_lengths=True, unquoted_underscores=True))
PY
iqtree2 -s $A/concat_tetranychoidea.fasta -p results/trees/pilot/pilot_tetranychoidea.best_scheme.nex -g work/stage5/constraint_pruned.nwk \
  --runs 10 -B 1000 -alrt 1000 -T 2 --seed 2026 --prefix work/stage5/constrained_tetra "$@"
cp work/stage5/constrained_tetra.treefile results/trees/tetranychoidea_constrained_ml.treefile
cp work/stage5/constrained_tetra.contree  results/trees/tetranychoidea_constrained.contree
cp work/stage5/constrained_tetra.iqtree   results/trees/tetranychoidea_constrained.iqtree
cp work/stage5/constrained_tetra.runtrees results/trees/tetranychoidea_constrained.runtrees
cp work/stage5/constrained_tetra.log      results/trees/tetranychoidea_constrained.log
cp work/stage5/constraint_pruned.nwk      results/trees/constraint_used.nwk
