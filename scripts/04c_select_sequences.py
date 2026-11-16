"""Stage 4c: filter GenBank records and choose one representative per taxon and marker.

Input:   data/clean/sequence_records.csv
Outputs: data/clean/accessions.csv                 chosen record per taxon x marker, with coverage and flags
         logs/sequences_excluded.csv               every record not used, with the reason
         results/tables/stage4_marker_coverage.csv taxa per marker by group
         intermediate files go to work/stage4 (not tracked in git)

Filters (in order)
  1 title: rRNA searches also return protein-coding and COI records, and COI searches return rRNA and ITS
    records of the same studies; these are dropped
    records that GenBank marks UNVERIFIED are dropped as well
  2 length: COI 300-1650 bp; 18S 500-2600 bp; 28S 400-5500 bp
  3 ambiguity: at most 1% non-ACGT bases
  4 COI reading frame: at least one frame with no internal stop codon (invertebrate mitochondrial code)
  5 duplicates: identical sequences within a taxon collapsed
  6 placement: each record is placed on a reference of the marker (Tetranychus urticae; Aceria tosichella
    for the rRNA genes of Eriophyoidea); records more than 45% different from the reference are dropped;
    (MAFFT --addfragments --keeplength) to find which part of the gene it covers
  7 conspecific check: records more than THRESH from the taxon's medoid are flagged and not chosen
Choice: among remaining records the one covering the gene columns shared by most taxa; ties go to a record
whose voucher or isolate also occurs in another marker, then to the longer record, then to the accession.
After the choice, a representative whose nearest neighbour belongs to another genus and lies closer than\nNEST is flagged for review. For COI a second record is added where it covers another part of the gene.
"""
import logging
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work/stage4"
WORK.mkdir(parents=True, exist_ok=True)
LEN = {"COI": (300, 1650), "18S": (500, 2600), "28S": (400, 5500)}
THRESH = {"COI": 0.08, "18S": 0.03, "28S": 0.05}
MIN_OVERLAP = 200
EXTRA_MIN = 250
NEST = {"COI": 0.03, "18S": 0.005, "28S": 0.01}  # nearest neighbour in another genus closer than this: review
REF_TAXON = "Tetranychus urticae"
REF_ERIO = "Aceria tosichella"          # rRNA reference for Eriophyoidea, too divergent to place on Tetranychus
ERIOPHYOIDEA = ["Eriophyidae", "Phytoptidae"]
REF_MAX = {"18S": 2600, "28S": 1500}    # reference limited to the commonly sequenced region
REF_DIST = 0.45                         # p-distance to the reference above which a placement is not homologous

logging.basicConfig(filename=ROOT / "logs/sequences_select.log", level=logging.INFO, filemode="w",
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()


def coi_frame(seq):
    """Return the in-frame, correctly oriented sequence with no internal stop, or None."""
    best = None
    for strand, s in (("+", seq), ("-", str(Seq(seq).reverse_complement()))):
        for frame in range(3):
            sub = s[frame:]
            sub = sub[: len(sub) - len(sub) % 3]
            prot = str(Seq(sub).translate(table=5))
            stops = prot[:-1].count("*")
            if stops == 0 and (best is None or len(sub) > len(best)):
                best = sub[:-3] if prot.endswith("*") else sub
    return best


def place(marker, ref_seq, recs):
    """Place records on the reference; return dict accession -> aligned string (reference length)."""
    ref = WORK / f"{marker}_ref.fasta"
    frag = WORK / f"{marker}_frag.fasta"
    ref.write_text(f">REF\n{ref_seq}\n")
    with open(frag, "w") as fh:
        for r in recs.itertuples():
            fh.write(f">{r.accession}\n{r.sequence}\n")
    cmd = ["mafft", "--quiet", "--thread", "4", "--keeplength", "--addfragments", str(frag)]
    if marker != "COI":
        cmd.insert(1, "--adjustdirection")
    out = subprocess.run(cmd + [str(ref)], capture_output=True, text=True, check=True).stdout
    aln, name = {}, None
    for line in out.splitlines():
        if line.startswith(">"):
            name = line[1:].replace("_R_", "")
            aln[name] = []
        else:
            aln[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in aln.items() if k != "REF"}


def to_array(aln, order):
    lut = np.zeros(256, dtype=np.uint8)
    for i, c in enumerate("ACGT", start=1):
        lut[ord(c)] = i
    return np.vstack([lut[np.frombuffer(aln[a].encode(), dtype=np.uint8)] for a in order])


def pdist_row(arr, i):
    both = (arr != 0) & (arr[i] != 0)
    n = both.sum(axis=1)
    diff = ((arr != arr[i]) & both).sum(axis=1)
    d = np.full(len(arr), np.nan)
    ok = n >= MIN_OVERLAP
    d[ok] = diff[ok] / n[ok]
    return d


def main():
    d = pd.read_csv(ROOT / "data/clean/sequence_records.csv", keep_default_na=False)
    d["reason"] = ""
    d.loc[d["marker_length"] == 0, "reason"] = d["extraction"]
    bad_title = d["marker"].ne("COI") & d["title"].str.contains(
        r"ribosomal protein|methyltransferase|cox1|cytochrome|mRNA|\(LOC\d+\)", case=False, regex=True)
    d.loc[bad_title & (d["reason"] == ""), "reason"] = "not an rRNA gene record"
    for m, pat in (("18S", r"18S|small subunit|SSU"), ("28S", r"28S|large subunit|LSU")):
        wrong = d["marker"].eq(m) & d["extraction"].eq("whole record") & ~d["title"].str.contains(pat, case=False, regex=True)
        d.loc[wrong & (d["reason"] == ""), "reason"] = f"not a {m} record"
    wgs = d["title"].str.contains("whole genome shotgun|scaffold|chromosome", case=False, regex=True)
    d.loc[wgs & (d["reason"] == ""), "reason"] = "genome assembly sequence, not a gene record"
    unver = d["title"].str.startswith("UNVERIFIED")
    d.loc[unver & (d["reason"] == ""), "reason"] = "marked UNVERIFIED by GenBank"
    # the COI search also returns rRNA and ITS records of the same studies (matched through other fields)
    not_coi = d["marker"].eq("COI") & d["extraction"].eq("whole record") & ~d["title"].str.contains(
        r"cytochrome|oxidase|\bCOI\b|COX1|\bCO1\b|mitochondri", case=False, regex=True)
    d.loc[not_coi & (d["reason"] == ""), "reason"] = "not a COI record"
    for m, (lo, hi) in LEN.items():
        sel = (d["marker"] == m) & (d["reason"] == "") & ~d["marker_length"].between(lo, hi)
        d.loc[sel, "reason"] = f"length outside {lo}-{hi} bp"
    amb = d["sequence"].map(lambda s: len(re.findall(r"[^ACGT]", s)) / max(len(s), 1))
    d.loc[(amb > 0.01) & (d["reason"] == ""), "reason"] = "more than 1% ambiguous bases"
    coi = (d["marker"] == "COI") & (d["reason"] == "")
    framed = d.loc[coi, "sequence"].map(coi_frame)
    d.loc[framed[framed.isna()].index, "reason"] = "COI: stop codons in every reading frame"
    d.loc[framed.dropna().index, "sequence"] = framed.dropna()
    d.loc[coi & (d["reason"] == "") & (d["sequence"].str.len() < LEN["COI"][0]), "reason"] = "length outside 300-1650 bp"
    good = d[d["reason"] == ""]
    dup = good.sort_values("accession").duplicated(["taxon", "marker", "sequence"])
    d.loc[dup[dup].index, "reason"] = "identical to another record of the same taxon"
    good = d[d["reason"] == ""].copy()
    good["tag"] = (good["voucher"] + "|" + good["isolate"] + "|" + good["strain"]).str.strip("|")

    # vouchers or isolates seen in more than one marker within a taxon
    tags = good[good["tag"] != ""].groupby(["taxon", "tag"])["marker"].nunique()
    shared = set(tags[tags > 1].index)

    chosen, placed_all = [], {}
    erio = good["family"].isin(ERIOPHYOIDEA)
    jobs = [("COI", good[good["marker"] == "COI"], REF_TAXON)]
    for m in ("18S", "28S"):
        jobs.append((m, good[(good["marker"] == m) & ~erio], REF_TAXON))
        jobs.append((m, good[(good["marker"] == m) & erio], REF_ERIO))
    for m, g, ref_taxon in jobs:
        g = g.reset_index(drop=True)
        refs = g[g["taxon"] == ref_taxon].sort_values(["marker_length", "accession"], ascending=[False, True])
        if m != "COI":
            refs = refs[refs["marker_length"] <= REF_MAX[m]]
        ref_row = refs.iloc[0]
        log.info("%s reference %s %s (%d bp); records to place %d", m, ref_taxon, ref_row.accession,
                 len(ref_row.sequence), len(g))
        aln = place(m, ref_row.sequence, g)
        g = g[g["accession"].isin(aln)].reset_index(drop=True)
        arr = to_array(aln, g["accession"].tolist())
        # records that do not match the reference region (another part of the gene, or another gene)
        ref_arr = to_array({"r": ref_row.sequence}, ["r"])[0]
        both = (arr != 0) & (ref_arr != 0)
        d_ref = ((arr != ref_arr) & both).sum(axis=1) / np.maximum(both.sum(axis=1), 1)
        off = d_ref > REF_DIST
        d.loc[d["accession"].isin(g.loc[off, "accession"]) & (d["marker"] == m) & (d["reason"] == ""),
              "reason"] = "does not match the reference region of the marker"
        g, arr = g[~off].reset_index(drop=True), arr[~off]
        g["covered"] = (arr != 0).sum(axis=1)
        g["cov_start"] = [int(np.argmax(r != 0)) + 1 if r.any() else 0 for r in arr]
        g["cov_end"] = [int(len(r) - np.argmax(r[::-1] != 0)) if r.any() else 0 for r in arr]
        low = g["covered"] < MIN_OVERLAP
        d.loc[d["accession"].isin(g.loc[low, "accession"]) & (d["marker"] == m) & (d["reason"] == ""),
              "reason"] = "less than 200 bp could be placed on the reference"
        keep = ~low.values
        g, arr = g[keep].reset_index(drop=True), arr[keep]

        # conspecific check against the taxon medoid
        g["flag"] = ""
        for taxon, idx in g.groupby("taxon").indices.items():
            if len(idx) < 3:
                continue
            sub = arr[idx]
            dm = np.vstack([pdist_row(sub, i) for i in range(len(idx))])
            med = np.nanmedian(np.where(np.eye(len(idx), dtype=bool), np.nan, dm), axis=1)
            if np.all(np.isnan(med)):
                continue
            medoid = int(np.nanargmin(med))
            far = dm[medoid] > THRESH[m]
            g.loc[idx[far], "flag"] = "divergent from conspecific records"
        div = g[g["flag"] != ""]
        d.loc[d["accession"].isin(div["accession"]) & (d["marker"] == m) & (d["reason"] == ""),
              "reason"] = f"more than {THRESH[m]:.0%} from the taxon's medoid sequence"

        # column weights: share of taxa with any usable record covering the column
        usable = g["flag"] == ""
        col_w = np.zeros(arr.shape[1])
        taxa_idx = g[usable].groupby("taxon").indices
        for taxon, idx in taxa_idx.items():
            col_w += (arr[g[usable].index[idx]] != 0).any(axis=0)
        col_w /= max(len(taxa_idx), 1)
        g["score"] = ((arr != 0) * col_w).sum(axis=1)
        g["shared_tag"] = [(t, tag) in shared for t, tag in zip(g["taxon"], g["tag"])]
        cand = g[usable].sort_values(["taxon", "score", "shared_tag", "marker_length", "accession"],
                                     ascending=[True, False, False, False, True])
        rep = cand.groupby("taxon").head(1).copy()
        rep["n_records_usable"] = rep["taxon"].map(cand.groupby("taxon").size())

        # nearest neighbour of each representative among all usable records of other taxa
        genus = g["taxon"].str.split().str[0].values
        notes = []
        for i in rep.index:
            dist = pdist_row(arr, i)
            other = (g["taxon"].values != g.at[i, "taxon"]) & usable.values & ~np.isnan(dist)
            if not other.any():
                notes.append(("", np.nan, ""))
                continue
            j = np.where(other)[0][np.argmin(dist[other])]
            note = ""
            same_genus = other & (genus == genus[i])
            if genus[j] != genus[i] and (not same_genus.any() or dist[j] < np.nanmin(dist[same_genus])) \
                    and dist[j] < NEST[m]:
                note = "nearest neighbour in another genus"
            elif m == "COI" and dist[j] < 0.01:
                note = "COI within 1% of another taxon"
            notes.append((g.at[j, "taxon"] + " " + g.at[j, "accession"], round(float(dist[j]), 4), note))
        rep["nearest_other_taxon"], rep["nearest_distance"], rep["review_flag"] = zip(*notes)
        rep["segment"] = 1
        chosen.append(rep)
        if m == "COI":
            # COI records fall in two regions of the gene (5' barcode region and 3' region). Where a second
            # usable record of the taxon adds at least EXTRA_MIN reference columns, it is kept as segment 2
            # and merged with the first after alignment (first record has priority where both have data).
            extra = []
            for i in rep.index:
                same = np.where((g["taxon"].values == g.at[i, "taxon"]) & usable.values)[0]
                same = same[same != i]
                if not len(same):
                    continue
                new_cols = ((arr[same] != 0) & (arr[i] == 0)).sum(axis=1)
                dist = pdist_row(arr, i)[same]
                overlap = ((arr[same] != 0) & (arr[i] != 0)).sum(axis=1)
                ok = (new_cols >= EXTRA_MIN) & ~((overlap >= MIN_OVERLAP) & (dist > 0.04))
                if ok.any():
                    j = same[ok][np.argmax(new_cols[ok])]
                    row = g.loc[[j]].copy()
                    row["n_records_usable"] = rep.at[i, "n_records_usable"]
                    # a second segment cannot always be compared with the first, so it is checked against
                    # other taxa instead: rejected when its nearest neighbour sits in another genus within NEST, or in any other taxon within 1%
                    dj = pdist_row(arr, j)
                    oth = (g["taxon"].values != g.at[j, "taxon"]) & usable.values & ~np.isnan(dj)
                    if oth.any():
                        k = np.where(oth)[0][np.argmin(dj[oth])]
                        if (genus[k] != genus[j] and dj[k] < NEST[m]) or dj[k] < 0.01:
                            log.info("COI segment 2 rejected: %s %s nearest %s at %.3f", g.at[j, "taxon"],
                                     g.at[j, "accession"], g.at[k, "taxon"], dj[k])
                            continue
                        row["nearest_other_taxon"] = g.at[k, "taxon"] + " " + g.at[k, "accession"]
                        row["nearest_distance"] = round(float(dj[k]), 4)
                    else:
                        row["nearest_other_taxon"], row["nearest_distance"] = "", np.nan
                    row["review_flag"] = "" if overlap[same == j][0] >= MIN_OVERLAP else "segment 2 not overlapping segment 1; checked against other taxa only"
                    row["segment"] = 2
                    extra.append(row)
            if extra:
                chosen.append(pd.concat(extra))
                rep = pd.concat([rep] + extra)
                log.info("COI: second segments added for %d taxa", len(extra))
        d.loc[d["accession"].isin(set(g["accession"]) - set(rep["accession"])) & (d["marker"] == m) &
              (d["reason"] == ""), "reason"] = "usable, another record of the taxon chosen"
        log.info("%s: usable records %d; taxa with a representative %d; review flags %d", m, int(usable.sum()),
                 len(rep), int((rep["review_flag"] != "").sum()))

    acc = pd.concat(chosen, ignore_index=True)
    cols = ["taxon", "family", "role", "marker", "segment", "accession", "organism", "marker_length", "covered", "cov_start",
            "cov_end", "extraction", "voucher", "isolate", "strain", "country", "host", "collection_date", "title",
            "n_records_found", "n_records_usable", "shared_tag", "nearest_other_taxon", "nearest_distance",
            "review_flag", "sequence"]
    acc[cols].sort_values(["marker", "family", "taxon", "segment"]).to_csv(ROOT / "data/clean/accessions.csv", index=False)
    ex = d[(d["reason"] != "")]
    ex[["taxon", "marker", "accession", "organism", "record_length", "marker_length", "title", "reason"]].to_csv(
        ROOT / "logs/sequences_excluded.csv", index=False)

    taxa = pd.read_csv(ROOT / "data/clean/analysis_taxa.csv")
    taxa = taxa[taxa["in_analysis"]]
    og = pd.read_csv(ROOT / "data/manual/outgroups.csv")
    allt = pd.concat([taxa.assign(group=np.where(taxa["resistant"], "resistant", "zero-resistance"))[["taxon", "family", "group"]],
                      og.assign(group="outgroup or anchor")[["taxon", "family", "group"]]], ignore_index=True)
    for m in ("COI", "18S", "28S"):
        allt[m] = allt["taxon"].isin(acc.loc[acc["marker"] == m, "taxon"])
    allt["n_markers"] = allt[["COI", "18S", "28S"]].sum(axis=1)
    allt.to_csv(ROOT / "data/clean/taxon_marker_matrix.csv", index=False)
    cov = allt.groupby(["family", "group"]).agg(taxa=("taxon", "size"), COI=("COI", "sum"), r18S=("18S", "sum"),
                                                r28S=("28S", "sum"), no_marker=("n_markers", lambda s: int((s == 0).sum())),
                                                all_three=("n_markers", lambda s: int((s == 3).sum()))).reset_index()
    cov.to_csv(ROOT / "results/tables/stage4_marker_coverage.csv", index=False)
    print(cov.to_string(index=False))
    print(ex["reason"].value_counts().to_string())
    print(acc.loc[acc["review_flag"] != "", ["taxon", "marker", "accession", "nearest_other_taxon", "nearest_distance"]].to_string(index=False))


if __name__ == "__main__":
    main()
