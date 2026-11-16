#!/usr/bin/env python3
"""Stage 5f: monophyly constraints for BEAST 2, taken from the constrained maximum-likelihood tree.

IQ-TREE accepts a backbone constraint on a subset of taxa and places the other taxa freely. BEAST 2 has no such
option: a monophyly constraint must list every member of the clade. Each backbone clade of Matsuda et al. (2018)
is therefore written with the membership it has in the constrained maximum-likelihood tree, that is, the taxa of
the alignment that are not in Matsuda et al. are assigned to the backbone clade in which IQ-TREE placed them.
Inside and between the constrained clades BEAST samples the topology freely.

Every backbone clade found in the tree is enforced (MIN_SUPPORT = 0). The ultrafast bootstrap value of each clade
is written to the log: all bootstrap trees obey the backbone, so this number measures how firmly the other taxa
are assigned to the clade. The uncertainty of these assignments is not carried by the BEAST trees; it is covered
by the sensitivity analysis on the dated bootstrap trees, in which the unassigned taxa move freely.

The backbone of Matsuda et al. is rooted on Bryobiinae by assumption (they had no other outgroup), so
Bryobiinae is not required to be one clade here: it is enforced only if the rooted tree shows it as one.

Usage: python -W ignore scripts/05f_beast_constraints.py [treefile]
Output: data/clean/beast_constraints_tetranychoidea.csv, logs/beast_constraints_tetranychoidea.csv
"""
import csv
import re
import sys
from pathlib import Path

import dendropy
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TREES = ROOT / "results" / "trees"
MIN_SUPPORT = 0.0
OUTGROUP = ["Leptotrombidium_deliense", "Cheyletus_malaccensis"]
treefile = Path(sys.argv[1]) if len(sys.argv) > 1 else TREES / "tetranychoidea_constrained_ml.treefile"


def support(node):
    m = re.match(r"^([\d.]+)/([\d.]+)$", node.label or "")
    return (float(m.group(1)), float(m.group(2))) if m else (None, None)


text = treefile.read_text().strip().splitlines()[0]
text = re.sub(r"^\[[^\]]*\]", "", text)
ml = dendropy.Tree.get(data=text, schema="newick", preserve_underscores=True, rooting="force-rooted")
ml.reroot_at_edge(ml.find_node_with_taxon_label(OUTGROUP[0]).edge, update_bipartitions=False)
mrca = ml.mrca(taxon_labels=OUTGROUP)
if {l.taxon.label for l in mrca.leaf_iter()} == set(OUTGROUP) and mrca is not ml.seed_node:
    ml.reroot_at_edge(mrca.edge, update_bipartitions=False)
tips = [l.taxon.label for l in ml.leaf_node_iter()]

family = {}
for r in csv.DictReader(open(ROOT / "data" / "clean" / "analysis_taxa.csv")):
    family[r["taxon"].replace(" ", "_")] = r["family"]
for r in csv.DictReader(open(ROOT / "data" / "manual" / "outgroups.csv")):
    family[r["taxon"].replace(" ", "_")] = r["family"]

bb = dendropy.Tree.get(path=str(TREES / "constraint_used.nwk"), schema="newick", preserve_underscores=True,
                       rooting="force-rooted")
bb_taxa = {l.taxon.label for l in bb.leaf_node_iter()}

rows, seen = [], {}


def clade(name, source, wanted, backbone_side=None):
    """Clade of the ML tree at the common ancestor of the wanted taxa."""
    node = ml.mrca(taxon_labels=sorted(wanted))
    members = frozenset(l.taxon.label for l in node.leaf_iter())
    if backbone_side is not None:
        ok = (members & bb_taxa) == set(backbone_side)
    else:
        ok = members == frozenset(wanted)
    a, b = support(node)
    row = {"clade": name, "source": source, "n_taxa": len(members),
           "n_backbone_taxa": len(members & bb_taxa), "n_assigned_taxa": len(members - bb_taxa),
           "clade_found_in_ml_tree": ok, "sh_alrt": a, "ufboot": b, "taxa": " ".join(sorted(members))}
    if members in seen:
        row["same_as"] = seen[members]
    else:
        seen[members] = name
    rows.append(row)


clade("ingroup", "all Tetranychoidea (roots the tree on the outgroups)", [t for t in tips if t not in OUTGROUP])
for fam in ["Tetranychidae", "Tenuipalpidae"]:
    clade(fam, "family", [t for t in tips if family.get(t) == fam])
i = 0
for nd in bb.postorder_internal_node_iter():
    if nd is bb.seed_node:
        continue
    side = [l.taxon.label for l in nd.leaf_iter()]
    i += 1
    genera = sorted({s.split("_")[0] for s in side})
    label = f"backbone{i:02d}_" + "_".join(g[:5] for g in genera[:3]) + f"_{len(side)}"
    clade(label, "Matsuda et al. (2018) backbone", side, backbone_side=side)

df = pd.DataFrame(rows)
if "same_as" not in df:
    df["same_as"] = ""
df["same_as"] = df["same_as"].fillna("")
df["enforced"] = (df["clade_found_in_ml_tree"] & (df["same_as"] == "") & (df["n_taxa"] < len(tips))
                  & ((df["source"] != "Matsuda et al. (2018) backbone") | (df["ufboot"].fillna(0) >= MIN_SUPPORT)))
# the ingroup and families are standard taxonomy; they are enforced whatever their support
cols = ["clade", "source", "enforced", "n_taxa", "n_backbone_taxa", "n_assigned_taxa", "clade_found_in_ml_tree",
        "sh_alrt", "ufboot", "same_as", "taxa"]
(ROOT / "logs").mkdir(exist_ok=True)
df[cols].to_csv(ROOT / "logs" / "beast_constraints_tetranychoidea.csv", index=False)
df.loc[df["enforced"], ["clade", "source", "n_taxa", "ufboot", "taxa"]].to_csv(
    ROOT / "data" / "clean" / "beast_constraints_tetranychoidea.csv", index=False)
print(df[cols[:-1]].to_string(index=False))
print(f"\nenforced: {int(df['enforced'].sum())} of {len(df)} clades; backbone clades not enforced: "
      f"{int(((df['source'] == 'Matsuda et al. (2018) backbone') & ~df['enforced'] & (df['same_as'] == '')).sum())}")
