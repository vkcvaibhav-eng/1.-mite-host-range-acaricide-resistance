#!/usr/bin/env python3
"""Stage 5g: ultrametric starting tree for BEAST 2 from the constrained maximum-likelihood tree.

The tree is rooted between Tetranychoidea and the two outgroups, node heights are set by the mean path length
from each node to its tips, every branch is kept at least 1e-4 long, and the tree is scaled so that its root height
equals the mean root-to-tip distance of the maximum-likelihood tree. The starting tree satisfies all monophyly
constraints because the constraints are taken from the same tree (05f). It only sets where the chains start.

Usage: python -W ignore scripts/05g_start_tree.py
Output: data/clean/beast_start_tree_tetranychoidea.nwk
"""
import re
from pathlib import Path

import dendropy

ROOT = Path(__file__).resolve().parents[1]
OUTGROUP = ["Leptotrombidium_deliense", "Cheyletus_malaccensis"]
text = re.sub(r"^\[[^\]]*\]", "",
              (ROOT / "results" / "trees" / "tetranychoidea_constrained_ml.treefile").read_text().strip())
t = dendropy.Tree.get(data=text, schema="newick", preserve_underscores=True, rooting="force-rooted")
for nd in t.preorder_node_iter():
    nd.label = None
ingroup = [l.taxon.label for l in t.leaf_node_iter() if l.taxon.label not in OUTGROUP]
t.reroot_at_edge(t.find_node_with_taxon_label(OUTGROUP[0]).edge, update_bipartitions=False)
mrca = t.mrca(taxon_labels=ingroup)
assert {l.taxon.label for l in mrca.leaf_iter()} == set(ingroup), "ingroup is not a clade in the ML tree"
half = (mrca.edge.length or 0.0) / 2
t.reroot_at_edge(mrca.edge, length1=half, length2=half, update_bipartitions=False)
t.suppress_unifurcations()
assert len(t.seed_node.child_nodes()) == 2

dist = {t.seed_node: 0.0}
for nd in t.preorder_node_iter():
    if nd.parent_node is not None:
        dist[nd] = dist[nd.parent_node] + (nd.edge.length or 0.0)
tips = list(t.leaf_node_iter())
target = sum(dist[l] for l in tips) / len(tips)

height = {}
for nd in t.postorder_node_iter():
    if nd.is_leaf():
        height[nd] = 0.0
    else:
        leaves = list(nd.leaf_iter())
        mpl = sum(dist[l] - dist[nd] for l in leaves) / len(leaves)
        height[nd] = max(mpl, max(height[c] for c in nd.child_node_iter()) + 1e-4)
scale = target / height[t.seed_node]
for nd in t.preorder_node_iter():
    if nd.parent_node is not None:
        nd.edge.length = (height[nd.parent_node] - height[nd]) * scale
t.seed_node.edge.length = None
out = ROOT / "data" / "clean" / "beast_start_tree_tetranychoidea.nwk"
out.write_text(t.as_string(schema="newick", suppress_rooting=True, unquoted_underscores=True,
                           suppress_internal_node_labels=True).strip() + "\n")
print(f"{out.relative_to(ROOT)}: {len(tips)} tips, root height {target:.3f}")
