"""Stage 4a: analysis taxa. One row per valid taxon after merging synonyms (logs/decisions.md 47, 49).

Inputs:  data/clean/host_breadth.csv, data/clean/resistance_by_species.csv, data/clean/candidate_species.csv
Output:  data/clean/analysis_taxa.csv
A taxon is resistant if any of its names has a resistance record; a zero-resistance taxon is a
Tetranychidae taxon with COI in GenBank, a crop host and no resistance record under any name.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
b = pd.read_csv(ROOT / "data/clean/host_breadth.csv")
res = pd.read_csv(ROOT / "data/clean/resistance_by_species.csv")
cand = pd.read_csv(ROOT / "data/clean/candidate_species.csv")
fam = dict(zip(cand["species"], cand["family"])) | dict(zip(res["species"], res["family"]))

rows = []
for valid, g in b.groupby("smw_valid_name"):
    names = sorted(set(g["species"]) | {valid})
    family = next((fam[n] for n in names if n in fam), "Tetranychidae")
    resistant = bool(g["has_resistance_record"].any())
    on_crop = bool(g["on_crop"].any())
    has_coi = bool(g["has_coi"].any())
    keep = resistant or (family == "Tetranychidae" and has_coi and on_crop)
    rows.append(dict(taxon=valid, family=family, search_names="; ".join(names), resistant=resistant,
                     on_crop=on_crop, has_coi_genbank=has_coi, in_analysis=keep,
                     n_cases=int(res.loc[res["species"].isin(names), "n_cases"].sum()),
                     host_source=g["host_source"].iloc[0]))
t = pd.DataFrame(rows).sort_values(["family", "taxon"])
t.to_csv(ROOT / "data/clean/analysis_taxa.csv", index=False)
a = t[t["in_analysis"]]
print(a.groupby(["family", "resistant"]).size().to_string())
print("taxa in analysis:", len(a), "| resistant:", int(a["resistant"].sum()), "| resistant without COI:",
      a.loc[a["resistant"] & ~a["has_coi_genbank"], "taxon"].tolist())
