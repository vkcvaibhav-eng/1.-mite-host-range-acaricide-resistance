#!/usr/bin/env python3
"""Stage 5d: check and summarise the maximum-likelihood trees.

- confirms that the constrained Tetranychoidea tree contains every split of the constraint that was used
- compares the constrained and the unconstrained (pilot) tree on the same data (approximately unbiased test)
- roots the trees on their outgroups and reports support and the monophyly of families and genera
- draws the trees for visual checking

Writes results/tables/stage5_ml_summary.csv, results/tables/stage5_genus_monophyly.csv,
results/trees/*_rooted.tre and results/figures/stage5_*.png. Run: python -W ignore scripts/05d_ml_summary.py
"""
import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import dendropy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TREES = ROOT / "results" / "trees"
WORK = ROOT / "work" / "stage5"
ALN = ROOT / "data" / "clean" / "alignments"
FITTED = TREES / "pilot" / "pilot_tetranychoidea.best_model.nex"   # scheme with the parameter values fitted in the pilot
OUTGROUPS = {
    "tetranychoidea": ["Leptotrombidium_deliense", "Cheyletus_malaccensis"],
    "eriophyoidea": ["Phytoptus_avellanae", "Nalepella_tsugifoliae"],
}
rows = []


def add(tree, item, value, note=""):
    rows.append({"tree": tree, "item": item, "value": value, "note": note})


def read(path, taxa=None):
    return dendropy.Tree.get(path=str(path), schema="newick", preserve_underscores=True,
                             taxon_namespace=taxa, rooting="force-unrooted")


def splits(tree, keep):
    """Non-trivial bipartitions of a tree restricted to the taxa in keep, as frozensets of the smaller side."""
    t = tree.clone(depth=1)
    t.retain_taxa_with_labels(keep)
    t.suppress_unifurcations()
    allt = frozenset(l.taxon.label for l in t.leaf_node_iter())
    out = set()
    for nd in t.internal_nodes():
        side = frozenset(l.taxon.label for l in nd.leaf_iter())
        if 1 < len(side) < len(allt) - 1:
            out.add(min(side, allt - side, key=lambda s: (len(s), sorted(s))))
    return out


def root_on(tree, outgroup):
    """Root on the branch leading to the outgroup taxa (their common ancestor if monophyletic, else the first)."""
    tree.is_rooted = True
    present = [o for o in outgroup if tree.find_node_with_taxon_label(o) is not None]
    tree.reroot_at_edge(tree.find_node_with_taxon_label(present[0]).edge, update_bipartitions=False)
    mono = True
    if len(present) > 1:
        mrca = tree.mrca(taxon_labels=present)
        mono = {l.taxon.label for l in mrca.leaf_iter()} == set(present)
        if mono and mrca is not tree.seed_node:
            tree.reroot_at_edge(mrca.edge, update_bipartitions=False)
    return mono


def support(node):
    """IQ-TREE writes 'SH-aLRT/UFBoot' as the node label."""
    lab = node.label or ""
    m = re.match(r"^([\d.]+)/([\d.]+)$", lab)
    return (float(m.group(1)), float(m.group(2))) if m else (None, None)


def clade_support(tree, members):
    """(is monophyletic, SH-aLRT, UFBoot) for a set of tip labels in a rooted tree."""
    members = set(members)
    if len(members) < 2:
        return None, None, None
    mrca = tree.mrca(taxon_labels=list(members))
    mono = {l.taxon.label for l in mrca.leaf_iter()} == members
    a, b = support(mrca) if mono else (None, None)
    return mono, a, b


def draw(tree, path, title, marked=(), width=11.0):
    """Rectangular phylogram with UFBoot values of 70 or more; tips in marked are printed in bold."""
    tips = [l for l in tree.leaf_node_iter()]
    y = {l: i for i, l in enumerate(tips)}
    x = {}
    for nd in tree.preorder_node_iter():
        x[nd] = 0.0 if nd.parent_node is None else x[nd.parent_node] + (nd.edge.length or 0.0)
    for nd in tree.postorder_node_iter():
        if not nd.is_leaf():
            ys = [y[c] for c in nd.child_node_iter()]
            y[nd] = (min(ys) + max(ys)) / 2
    fig, ax = plt.subplots(figsize=(width, max(4.0, 0.16 * len(tips) + 1)))
    for nd in tree.preorder_node_iter():
        if nd.parent_node is not None:
            ax.plot([x[nd.parent_node], x[nd]], [y[nd], y[nd]], color="black", lw=0.7)
        if not nd.is_leaf():
            ys = [y[c] for c in nd.child_node_iter()]
            ax.plot([x[nd], x[nd]], [min(ys), max(ys)], color="black", lw=0.7)
            ufb = support(nd)[1]
            if ufb is not None and ufb >= 70:
                ax.text(x[nd], y[nd], f"{ufb:.0f}", fontsize=4.5, ha="right", va="bottom", color="#555555")
    xmax = max(x.values()) or 1.0
    for l in tips:
        name = l.taxon.label.replace("_", " ")
        ax.text(x[l] + 0.005 * xmax, y[l], name, fontsize=6, va="center", style="italic",
                weight="bold" if l.taxon.label in marked else "normal")
    ax.set_ylim(len(tips), -1)
    ax.set_xlim(0, xmax * 1.35)
    ax.axis("off")
    ax.set_title(title, fontsize=9, loc="left")
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


taxa_tab = pd.read_csv(ROOT / "data" / "clean" / "analysis_taxa.csv")
resistant = {t.replace(" ", "_") for t in taxa_tab.loc[taxa_tab["resistant"], "taxon"]}
family = {t.replace(" ", "_"): f for t, f in zip(taxa_tab["taxon"], taxa_tab["family"])}
for r in csv.DictReader(open(ROOT / "data" / "manual" / "outgroups.csv")):
    family[r["taxon"].replace(" ", "_")] = r["family"]
(ROOT / "results" / "figures").mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- Tetranychoidea
con = read(TREES / "tetranychoidea_constrained_ml.treefile")
tips = [l.taxon.label for l in con.leaf_node_iter()]
used = read(TREES / "constraint_used.nwk")
used_taxa = [l.taxon.label for l in used.leaf_node_iter()]
want = splits(used, used_taxa)
have = splits(con, used_taxa)
full = read(ROOT / "data" / "manual" / "constraint_matsuda2018.nwk")
add("tetranychoidea", "taxa in tree", len(tips))
add("tetranychoidea", "taxa in the constraint of Matsuda et al. (2018) as read from Figure 1",
    sum(1 for _ in full.leaf_node_iter()))
add("tetranychoidea", "constraint taxa present in the alignment", len(used_taxa))
add("tetranychoidea", "constraint splits enforced", len(want))
add("tetranychoidea", "constraint splits found in the constrained tree", len(want & have),
    "must equal the number enforced")
assert want <= have, "constrained tree violates the constraint"

unc = read(TREES / "pilot" / "pilot_tetranychoidea.treefile")
add("tetranychoidea", "constraint splits already present in the unconstrained tree",
    len(want & splits(unc, used_taxa)), f"of {len(want)}")

# approximately unbiased test, constrained against unconstrained tree. Both trees are scored with the model
# parameters fitted on the unconstrained tree in the pilot (read from best_model.nex) and their own branch lengths.
# Re-estimating the free-rate parameters in every run gives likelihoods that differ by tens of units for the same
# tree, so the values printed by the searches are not comparable with each other.
WORK.mkdir(parents=True, exist_ok=True)
UNC = TREES / "pilot" / "pilot_tetranychoidea.treefile"


def score(trees_text, prefix):
    """log-likelihood, deltaL and p-AU of each tree in trees_text under the fitted pilot parameters."""
    f = WORK / f"{prefix}.nwk"
    f.write_text(trees_text)
    subprocess.run(["iqtree2", "-s", str(ALN / "concat_tetranychoidea.fasta"), "-p", str(FITTED), "-te", str(UNC),
                    "-z", str(f), "-zb", "10000", "-au", "-T", "2", "--seed", "2026",
                    "--prefix", str(WORK / prefix), "-redo", "-quiet"], check=True)
    res = []
    for line in (WORK / f"{prefix}.iqtree").read_text().splitlines():
        m = re.match(r"^\s*(\d+)\s+(-[\d.]+)\s+([\d.]+)\s+(.*[+-])\s*$", line)
        if m:
            res.append((float(m.group(2)), float(m.group(3)), re.findall(r"[\d.eE+-]+(?=\s*[+-])", m.group(4))[-1]))
    return res


au = score(UNC.read_text().strip() + "\n" + (TREES / "tetranychoidea_constrained_ml.treefile").read_text().strip() + "\n",
           "au_test")
assert len(au) == 2
for (logl, dl, pau), nm in zip(au, ["unconstrained (pilot)", "constrained"]):
    add("tetranychoidea", f"log-likelihood under common parameters, {nm} tree", f"{logl:.2f}")
    add("tetranychoidea", f"approximately unbiased test p, {nm} tree", pau,
        "10,000 RELL replicates; p below 0.05 means the tree is rejected by the three-gene data")
add("tetranychoidea", "log-likelihood cost of the constraint", f"{au[0][0] - au[1][0]:.2f}")

# the ten independent searches
log = (TREES / "tetranychoidea_constrained.log").read_text()
scores = [float(x) for x in re.findall(r"BEST SCORE FOUND : (-[\d.]+)", log)]
seeds = re.findall(r"START RUN NUMBER \d+ \(seed: (\d+)\)", log)
run_lines = [re.sub(r"^\[[^\]]*\]", "", l) for l in
             (TREES / "tetranychoidea_constrained.runtrees").read_text().splitlines() if l.strip()]
common = score("\n".join(run_lines) + "\n", "runs_common")
runs = [dendropy.Tree.get(data=l, schema="newick", preserve_underscores=True,
                          taxon_namespace=con.taxon_namespace, rooting="force-unrooted") for l in run_lines]
best_splits = splits(con, tips)
shared = [len(best_splits & splits(r, tips)) for r in runs]
cl = [c[0] for c in common]
pd.DataFrame({"search": range(1, len(scores) + 1), "seed": seeds, "log_likelihood_reported_by_search": scores,
              "log_likelihood_common_parameters": cl, "splits_shared_with_best_tree": shared,
              "splits_in_best_tree": len(best_splits)}).to_csv(
    ROOT / "results" / "tables" / "stage5_ml_searches.csv", index=False)
best = cl.index(max(cl)) + 1
add("tetranychoidea", "independent constrained searches", len(scores))
add("tetranychoidea", "search kept by IQ-TREE (highest reported log-likelihood)", scores.index(max(scores)) + 1)
add("tetranychoidea", "search with the highest log-likelihood under common parameters", best,
    "must be the search kept")
assert best == scores.index(max(scores)) + 1
add("tetranychoidea", "log-likelihood under common parameters, best search", f"{max(cl):.2f}")
add("tetranychoidea", "log-likelihood under common parameters, worst search", f"{min(cl):.2f}")
others = sorted(x for i, x in enumerate(shared) if i != best - 1)
add("tetranychoidea", "splits of the best tree found in the other nine searches, median", others[len(others) // 2],
    f"of {len(best_splits)}; range {others[0]}-{others[-1]}")

og_mono = root_on(con, OUTGROUPS["tetranychoidea"])
add("tetranychoidea", "outgroups form one group at the root", og_mono,
    "if False the tree is rooted on Leptotrombidium deliense alone")
ing = [t for t in tips if t not in OUTGROUPS["tetranychoidea"]]
mono, a, b = clade_support(con, ing)
add("tetranychoidea", "Tetranychoidea (all taxa except the two outgroups) monophyletic", mono,
    f"SH-aLRT {a}, ultrafast bootstrap {b}" if mono else "")
ufb = [support(n)[1] for n in con.internal_nodes() if support(n)[1] is not None]
add("tetranychoidea", "internal nodes with support values", len(ufb))
add("tetranychoidea", "nodes with ultrafast bootstrap 95 or more", sum(u >= 95 for u in ufb))
add("tetranychoidea", "nodes with ultrafast bootstrap below 70", sum(u < 70 for u in ufb))
for fam in ["Tetranychidae", "Tenuipalpidae"]:
    mono, a, b = clade_support(con, [t for t in tips if family.get(t) == fam])
    add("tetranychoidea", f"{fam} monophyletic", mono, f"SH-aLRT {a}, ultrafast bootstrap {b}" if mono else "")

genera = defaultdict(list)
for t in tips:
    if family.get(t) in ("Tetranychidae", "Tenuipalpidae"):
        genera[t.split("_")[0]].append(t)
grow = []
for g, members in sorted(genera.items()):
    mono, a, b = clade_support(con, members)
    grow.append({"genus": g, "taxa_in_tree": len(members),
                 "monophyletic": "" if mono is None else mono,
                 "sh_alrt": "" if a is None else a, "ufboot": "" if b is None else b})
pd.DataFrame(grow).to_csv(ROOT / "results" / "tables" / "stage5_genus_monophyly.csv", index=False)

con.ladderize(ascending=True)
con.write(path=str(TREES / "tetranychoidea_constrained_rooted.tre"), schema="newick", unquoted_underscores=True)
draw(con, ROOT / "results" / "figures" / "stage5_ml_tree_tetranychoidea.png",
     "Tetranychoidea, constrained maximum-likelihood tree (bold: taxa with resistance records; "
     "numbers: ultrafast bootstrap of 70 or more)", marked=resistant)

# constraint as read from the figure, taxa of the alignment in bold
full.is_rooted = True
for nd in full.preorder_node_iter():
    nd.edge.length = 1.0 if nd.parent_node is not None else 0.0
depth = {}
for nd in full.postorder_node_iter():
    depth[nd] = 0 if nd.is_leaf() else 1 + max(depth[c] for c in nd.child_node_iter())
for nd in full.preorder_node_iter():
    if nd.parent_node is not None:
        nd.edge.length = float(depth[nd.parent_node] - depth[nd])
draw(full, ROOT / "results" / "figures" / "stage5_constraint_matsuda2018.png",
     "Backbone read from Matsuda et al. (2018), Figure 1 (bold: taxa present in our alignment; "
     "unresolved groups are left free)", marked=set(used_taxa), width=9.0)

# ---------------------------------------------------------------- Eriophyoidea
ef = TREES / "eriophyoidea_ml.treefile"
if ef.exists():
    er = read(ef)
    etips = [l.taxon.label for l in er.leaf_node_iter()]
    emono = root_on(er, OUTGROUPS["eriophyoidea"])
    add("eriophyoidea", "taxa in tree", len(etips))
    add("eriophyoidea", "outgroups form one group at the root", emono,
        "if False the tree is rooted on Phytoptus avellanae alone")
    mono, a, b = clade_support(er, [t for t in etips if family.get(t) == "Eriophyidae"])
    add("eriophyoidea", "Eriophyidae monophyletic", mono, f"SH-aLRT {a}, ultrafast bootstrap {b}" if mono else "")
    escores = re.findall(r"BEST SCORE FOUND : (-[\d.]+)", (TREES / "eriophyoidea.log").read_text())
    add("eriophyoidea", "independent searches", len(escores),
        f"log-likelihood {min(map(float, escores)):.2f} to {max(map(float, escores)):.2f}")
    eufb = [support(n)[1] for n in er.internal_nodes() if support(n)[1] is not None]
    add("eriophyoidea", "internal nodes with support values", len(eufb))
    add("eriophyoidea", "nodes with ultrafast bootstrap 95 or more", sum(u >= 95 for u in eufb))
    add("eriophyoidea", "nodes with ultrafast bootstrap below 70", sum(u < 70 for u in eufb))
    m = re.search(r"Log-likelihood of the tree: (-[\d.]+)", (TREES / "eriophyoidea_ml.iqtree").read_text())
    if m:
        add("eriophyoidea", "log-likelihood", m.group(1))
    lens = [e.length for e in er.postorder_edge_iter() if e.length]
    add("eriophyoidea", "tree length (substitutions per site)", f"{sum(lens):.2f}")
    er.ladderize(ascending=True)
    er.write(path=str(TREES / "eriophyoidea_rooted.tre"), schema="newick", unquoted_underscores=True)
    draw(er, ROOT / "results" / "figures" / "stage5_ml_tree_eriophyoidea.png",
         "Eriophyoidea, maximum-likelihood tree (bold: taxa with resistance records)", marked=resistant, width=8.0)

lens = [e.length for e in con.postorder_edge_iter() if e.length]
add("tetranychoidea", "tree length (substitutions per site)", f"{sum(lens):.2f}")
pd.DataFrame(rows).to_csv(ROOT / "results" / "tables" / "stage5_ml_summary.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False))
