"""Stage 3c: host breadth per mite species.

Inputs:  data/clean/host_records.csv
         data/raw/GBOTB_smith_brown_2018.tre  (Smith and Brown 2018, American Journal of Botany 105: 302-314;
                                               GBOTB tree, release v0.1, md5 ce280a81ab824dc7cc83f395a7c7b54e)
Output:  data/clean/host_breadth.csv
         results/tables/stage3_host_breadth_extremes.csv (20 broadest and 20 narrowest by host genera)
         logs/host_genera_not_in_tree.csv

Metrics per mite species
  n_host_records        host-plant records (plant name x mite) in Spider Mites Web
  n_references          distinct references behind those records
  n_host_species        distinct WFO accepted species
  n_host_genera, n_host_families
  n_crop_genera         host genera in the FAO crop list; on_crop = n_crop_genera > 0
  rarefied_genera_10    expected number of host genera in 10 host-reference records (Hurlbert);
                        missing when a species has fewer than 10 records
  host_pd_genus         Faith's phylogenetic diversity of host genera (sum of branch lengths, in million
                        years, of the subtree linking one tip per genus, root excluded); missing when fewer
                        than 2 genera are in the tree
"""
import logging
from math import comb
from pathlib import Path

import dendropy
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
K = 10
logging.basicConfig(filename=ROOT / "logs/host_breadth.log", level=logging.INFO, filemode="w",
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()


def hurlbert(counts, k):
    n = sum(counts)
    if n < k:
        return np.nan
    denom = comb(n, k)
    return float(sum(1 - comb(n - c, k) / denom for c in counts))


def main():
    h = pd.read_csv(ROOT / "data/clean/host_records.csv", dtype=str, keep_default_na=False)
    h["crop_genus"] = h["crop_genus"].eq("True")
    h = h[h["plant"] != ""]  # species whose only record was "Unknown" keep no host rows
    # check against the counts printed on each Spider Mites Web page ("Host plants (N plants, ...)")
    raw = pd.read_csv(ROOT / "data/raw/smw_host_records_2026-11-16.csv", dtype=str, keep_default_na=False)
    stated = raw.drop_duplicates("smw_id").set_index("smw_id")["host_header"].str.extract(r"\((\d+) plant")[0].astype(int)
    got = raw.groupby("smw_id").size()
    bad = (stated != got.reindex(stated.index)).sum()
    log.info("species whose parsed host count differs from the page count: %d", int(bad))
    if bad:
        raise ValueError("parsed host counts differ from the Spider Mites Web page counts")

    # one row per host x reference for rarefaction; records with no reference count once
    pairs = h.assign(ref=h["ref_ids"].str.split()).explode("ref")
    pairs["ref"] = pairs["ref"].fillna("none")

    tree = dendropy.Tree.get(path=ROOT / "data/raw/GBOTB_smith_brown_2018.tre", schema="newick",
                             preserve_underscores=True, rooting="force-rooted")
    tip_genus = {}
    for leaf in tree.leaf_node_iter():
        g = leaf.taxon.label.split("_")[0]
        tip_genus.setdefault(g, leaf.taxon.label)
    host_genera = sorted(set(h.loc[h["host_genus"] != "", "host_genus"]))
    in_tree = [g for g in host_genera if g in tip_genus]
    missing = [g for g in host_genera if g not in tip_genus]
    pd.DataFrame({"host_genus": missing}).to_csv(ROOT / "logs/host_genera_not_in_tree.csv", index=False)
    log.info("host genera %d; in tree %d; missing %d", len(host_genera), len(in_tree), len(missing))
    keep = {tip_genus[g] for g in in_tree}
    tree.retain_taxa_with_labels(keep)
    tree.suppress_unifurcations()
    label_of = {g: tip_genus[g] for g in in_tree}

    def faith_pd(genera):
        labels = {label_of[g] for g in genera if g in label_of}
        if len(labels) < 2:
            return np.nan, len(labels)
        leaves = [n for n in tree.leaf_node_iter() if n.taxon.label in labels]
        # edges on paths from each leaf up to the most recent common ancestor
        mrca = tree.mrca(taxon_labels=list(labels))
        edges = set()
        for lf in leaves:
            node = lf
            while node is not mrca:
                edges.add(node)
                node = node.parent_node
        return float(sum(n.edge.length or 0 for n in edges)), len(labels)

    rows = []
    for sp, g in h.groupby("species"):
        genera = set(g.loc[g["host_genus"] != "", "host_genus"])
        accepted = set(g.loc[g["accepted_name"] != "", "accepted_name"])
        pr = pairs[pairs["species"] == sp]
        counts = pr[pr["host_genus"] != ""].groupby("host_genus").size().tolist()
        pd_val, n_tree = faith_pd(genera)
        rows.append(dict(
            species=sp, smw_id=g["smw_id"].iloc[0], smw_valid_name=g["smw_valid_name"].iloc[0],
            n_host_records=len(g), n_references=len({r for x in g["ref_ids"] for r in x.split()}),
            n_host_species=len(accepted), n_host_genera=len(genera),
            n_host_families=g.loc[g["host_family"] != "", "host_family"].nunique(),
            n_crop_genera=g.loc[g["crop_genus"], "host_genus"].nunique(),
            rarefied_genera_10=hurlbert(counts, K), host_pd_genus=pd_val, n_genera_in_tree=n_tree))
    b = pd.DataFrame(rows)
    nm = pd.read_csv(ROOT / "data/manual/smw_name_map.csv", dtype=str)
    b = nm[["species", "smw_id", "smw_valid_name"]].merge(b.drop(columns=["smw_id", "smw_valid_name"]), on="species", how="left")
    b = b.fillna({c: 0 for c in ["n_host_records", "n_references", "n_host_species", "n_host_genera",
                                 "n_host_families", "n_crop_genera", "n_genera_in_tree"]})
    b["on_crop"] = b["n_crop_genera"] > 0
    b.to_csv(ROOT / "data/clean/host_breadth.csv", index=False)

    # check table: 20 broadest and 20 narrowest by host genera (distinct Spider Mites Web taxa only)
    u = b.drop_duplicates("smw_id").sort_values(["n_host_genera", "n_host_records"], ascending=False)
    ext = pd.concat([u.head(20).assign(group="broadest"), u.tail(20).assign(group="narrowest")])
    cols = ["group", "smw_valid_name", "n_host_records", "n_references", "n_host_species", "n_host_genera",
            "n_host_families", "n_crop_genera", "rarefied_genera_10", "host_pd_genus"]
    ext[cols].to_csv(ROOT / "results/tables/stage3_host_breadth_extremes.csv", index=False)
    cand = pd.read_csv(ROOT / "data/clean/candidate_species.csv")
    res = pd.read_csv(ROOT / "data/clean/resistance_by_species.csv")
    tet = set(cand.loc[cand["family"] == "Tetranychidae", "species"])
    b["has_coi"] = b["species"].isin(set(cand["species"]))
    b["has_resistance_record"] = b["species"].isin(set(res["species"]))
    zero = b[b["species"].isin(tet) & b["on_crop"] & ~b["has_resistance_record"]]
    summ = pd.DataFrame([
        ("spider mite names with host data", len(b)),
        ("distinct Spider Mites Web taxa", b["smw_id"].nunique()),
        ("host-plant records (distinct taxa)", int(raw.shape[0])),
        ("names with at least one crop host", int(b["on_crop"].sum())),
        ("resistant spider mite names", int(b["has_resistance_record"].sum())),
        ("zero-resistance species: COI, crop host, no resistance record", len(zero)),
        ("host genera", len(host_genera)),
        ("host genera not in the plant tree", len(missing)),
    ], columns=["item", "n"])
    summ.to_csv(ROOT / "results/tables/stage3_summary.csv", index=False)
    b["zero_resistance_species"] = b["species"].isin(zero["species"])
    b.to_csv(ROOT / "data/clean/host_breadth.csv", index=False)
    print(summ.to_string(index=False))
    log.info("species %d; distinct SMW taxa %d; on crop %d", len(b), b["smw_id"].nunique(), int(b["on_crop"].sum()))
    print(ext[cols].to_string(index=False))


if __name__ == "__main__":
    main()
