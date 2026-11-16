"""Stage 3b: clean Spider Mites Web host records and standardise plant names against World Flora Online.

Inputs
  data/raw/smw_host_records_2026-11-16.csv   host records saved from Spider Mites Web (checksum below)
  data/manual/smw_name_map.csv               mite names -> Spider Mites Web species id
  data/raw/tenuipalpidae_db_hosts_2026-11-16.csv  Brevipalpus host records from the Tenuipalpidae Database
                                             (Castro et al. 2020), saved from the user's Chrome
  data/manual/crop_genera.csv                crop genera from the FAO Indicative Crop Classification
  WFO backbone (World Flora Online Plant List June 2026, doi:10.5281/zenodo.20782718,
  _DwC_backbone_R.zip, md5 0e4486945cd9f7af548ca87eb9a870ed), path given by WFO_CLASSIFICATION
Output
  data/clean/host_records.csv  one row per mite species x host record, with WFO accepted names
  logs/host_names_unmatched.csv
"""
import logging
import os
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW_FILE = ROOT / "data/raw/smw_host_records_2026-11-16.csv"
EXPECTED_ROWS, EXPECTED_ID_SUM = 7716, 5967422
WFO = Path(os.environ.get("WFO_CLASSIFICATION", "/tmp/wfo/classification.csv"))

logging.basicConfig(filename=ROOT / "logs/host_records.log", level=logging.INFO, filemode="w",
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()


def clean_name(s):
    """Return (genus, epithet or '') from a host name as listed."""
    s = re.sub(r"\s+", " ", s.replace("×", " ").replace(" x ", " ")).strip()
    s = re.sub(r"\b(cv|var|subsp|ssp|f)\..*$", "", s).strip()
    parts = s.split(" ")
    genus = parts[0].capitalize() if parts else ""
    ep = parts[1].lower() if len(parts) > 1 else ""
    if ep in ("sp", "sp.", "spp", "spp.", "sps", "hybrid", "hybr.", "") or not re.fullmatch(r"[a-z\-]+", ep):
        ep = ""
    return genus, ep


def main():
    raw = pd.read_csv(RAW_FILE, dtype=str, keep_default_na=False)
    if len(raw) != EXPECTED_ROWS or raw["smw_id"].astype(int).sum() != EXPECTED_ID_SUM:
        raise ValueError("Spider Mites Web host file does not match its checksum")
    raw["n_refs"] = raw["ref_ids"].str.split().str.len().fillna(0).astype(int)
    # Brevipalpus records from the Tenuipalpidae Database, put in the same shape
    tp = pd.read_csv(ROOT / "data/raw/tenuipalpidae_db_hosts_2026-11-16.csv", dtype=str, keep_default_na=False)
    if len(tp) != 669:
        raise ValueError("Tenuipalpidae Database file does not have the expected 669 rows")
    fam = tp["plant_family"].str.capitalize().replace({"Leguminosae": "Fabaceae", "Compositae": "Asteraceae",
                                                       "Gramineae": "Poaceae", "Palmae": "Arecaceae",
                                                       "Cruciferae": "Brassicaceae", "Labiatae": "Lamiaceae",
                                                       "Umbelliferae": "Apiaceae", "Guttiferae": "Clusiaceae"})
    tp = pd.DataFrame({"smw_id": "TPDB:" + tp["mite_species"], "species_title": tp["mite_species"],
                       "host_header": "Tenuipalpidae Database", "plant_family": fam, "plant": tp["plant"],
                       "ref_ids": tp["references"].str.replace("; ", "|", regex=False)})
    tp["n_refs"] = tp["ref_ids"].str.split("|").str.len()
    raw = pd.concat([raw, tp], ignore_index=True)
    unknown = raw["plant"].str.lower().isin(["unknown", ""]) | raw["plant_family"].str.lower().eq("unknown")
    log.info("records %d; unknown host %d", len(raw), int(unknown.sum()))
    rec = raw[~unknown].copy()
    rec[["genus_listed", "epithet_listed"]] = rec["plant"].apply(lambda s: pd.Series(clean_name(s)))

    cols = ["taxonID", "scientificName", "taxonRank", "family", "genus", "specificEpithet",
            "taxonomicStatus", "acceptedNameUsageID", "majorGroup"]
    # read the large backbone in chunks: keep species and genus rows of the host genera,
    # then add the accepted names those rows point to (which may sit in other genera)
    genera = set(rec["genus_listed"])
    read = dict(sep="\t", usecols=cols, dtype=str, keep_default_na=False, quoting=3, encoding="utf-8",
                on_bad_lines="skip", chunksize=200_000)
    parts = []
    for ch in pd.read_csv(WFO, **read):
        ch = ch[ch["taxonRank"].isin(["species", "genus"])]
        ch["scientificName"] = ch["scientificName"].str.strip('"')
        gcol = ch["genus"].where(ch["taxonRank"] == "species", ch["scientificName"])
        parts.append(ch[gcol.isin(genera)])
    w = pd.concat(parts, ignore_index=True)
    need = set(w["acceptedNameUsageID"]) - set(w["taxonID"]) - {""}
    extra = []
    for ch in pd.read_csv(WFO, **read):
        ch = ch[ch["taxonID"].isin(need)]
        ch["scientificName"] = ch["scientificName"].str.strip('"')
        extra.append(ch)
    w = pd.concat([w] + extra, ignore_index=True)
    log.info("WFO rows (species and genus rank): %d", len(w))
    by_id = w.set_index("taxonID")
    status_rank = {"Accepted": 0, "Synonym": 1, "Unchecked": 2}
    w["srank"] = w["taxonomicStatus"].map(status_rank).fillna(3)

    def resolve(row):
        """Follow a synonym to its accepted name; return (id, name, genus, family, status)."""
        if row["taxonomicStatus"] == "Synonym" and row["acceptedNameUsageID"] in by_id.index:
            a = by_id.loc[row["acceptedNameUsageID"]]
            if isinstance(a, pd.DataFrame):
                a = a.iloc[0]
            return row["acceptedNameUsageID"], a["scientificName"], a["genus"], a["family"], "synonym resolved"
        return row["taxonID"], row["scientificName"], row["genus"], row["family"], row["taxonomicStatus"].lower()

    sp = w[w["taxonRank"] == "species"].sort_values("srank")
    # hybrid names in WFO carry a multiplication sign ("Citrus ×aurantium"); index them without it
    sp = sp.assign(key=sp["scientificName"].str.replace("×", "", regex=False).str.replace(r"\s+", " ", regex=True).str.strip())
    sp_index = {k: g for k, g in sp.groupby("key")}
    gen = w[w["taxonRank"] == "genus"].sort_values("srank")
    gen_index = {k: g for k, g in gen.groupby("scientificName")}

    out = []
    for name in sorted(set(zip(rec["genus_listed"], rec["epithet_listed"], rec["plant_family"]))):
        g, e, fam_smw = name
        res = None
        if e:
            cand = sp_index.get(f"{g} {e}")
            if cand is not None:
                if len(cand) > 1 and (cand["family"] == fam_smw).any():
                    cand = cand[cand["family"] == fam_smw]
                res = resolve(cand.iloc[0]) + ("species",)
        if res is None:
            cand = gen_index.get(g)
            if cand is not None:
                if len(cand) > 1 and (cand["family"] == fam_smw).any():
                    cand = cand[cand["family"] == fam_smw]
                r = resolve(cand.iloc[0])
                res = (r[0], None, r[2] or r[1], r[3], r[4],
                       "genus only" if not e else "genus only (species not found)")
        if res is None:
            res = (None, None, None, None, "unmatched", "unmatched")
        wid, acc, accg, accf, status, level = res
        if level == "species":
            accg = acc.split()[0] if acc else accg
        out.append(dict(genus_listed=g, epithet_listed=e, plant_family=fam_smw, wfo_id=wid,
                        accepted_name=acc if level == "species" else None, accepted_genus=accg,
                        accepted_family=accf, wfo_status=status, match_level=level))
    names = pd.DataFrame(out)
    rec = rec.merge(names, on=["genus_listed", "epithet_listed", "plant_family"], how="left")
    # family fallback: use the Spider Mites Web family when WFO has no match
    rec["host_family"] = rec["accepted_family"].where(rec["accepted_family"].notna() & (rec["accepted_family"] != ""),
                                                      rec["plant_family"])
    rec["host_genus"] = rec["accepted_genus"].where(rec["accepted_genus"].notna(), rec["genus_listed"])
    # records naming only a family or a higher group carry no genus
    higher = rec["genus_listed"].str.endswith("aceae") | rec["genus_listed"].isin(["Bryophytes"])
    rec.loc[higher, ["host_genus", "match_level"]] = [None, "family only"]

    crops = pd.read_csv(ROOT / "data/manual/crop_genera.csv")
    rec["crop_genus"] = rec["host_genus"].isin(crops["genus"]) | rec["genus_listed"].isin(crops["genus"])

    nm = pd.read_csv(ROOT / "data/manual/smw_name_map.csv", dtype=str)
    bre = sorted(x[5:] for x in raw["smw_id"].unique() if x.startswith("TPDB:"))
    nm = pd.concat([nm, pd.DataFrame({"species": bre, "smw_id": ["TPDB:" + b for b in bre], "smw_valid_name": bre,
                                      "match": "Tenuipalpidae Database"})], ignore_index=True)
    rec = nm.merge(rec, on="smw_id", how="left")
    keep = ["species", "smw_id", "smw_valid_name", "plant", "plant_family", "genus_listed", "epithet_listed",
            "wfo_id", "accepted_name", "host_genus", "host_family", "wfo_status", "match_level", "crop_genus",
            "n_refs", "ref_ids"]
    rec[keep].to_csv(ROOT / "data/clean/host_records.csv", index=False)

    un = names[names["match_level"].str.startswith("unmatched") | names["match_level"].str.contains("not found")]
    un.to_csv(ROOT / "logs/host_names_unmatched.csv", index=False)
    log.info("distinct host names %d; match levels: %s", len(names), names["match_level"].value_counts().to_dict())
    print(names["match_level"].value_counts().to_string())
    print("records written:", len(rec))


if __name__ == "__main__":
    main()
