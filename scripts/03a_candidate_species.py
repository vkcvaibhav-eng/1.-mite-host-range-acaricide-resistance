"""Stage 3a: candidate species with COI sequences in GenBank for the three mite families.

For each family, find nuccore records matching COI terms under the family's taxonomy subtree,
fetch their organism names, and list formally named species (no "sp.", "cf.", "aff.").
Output: data/raw/ncbi_coi_species_<date>.csv (one row per species, with record counts)
        data/clean/candidate_species.csv (named species, flag for resistance record)
"""
import datetime as dt
import logging
import re
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
FAMILIES = {"Tetranychidae": 32262, "Tenuipalpidae": 50007, "Eriophyidae": 70320}
COI = ('(COI[All Fields] OR COX1[All Fields] OR "cytochrome oxidase subunit I"[All Fields] '
       'OR "cytochrome c oxidase subunit 1"[All Fields] OR "cytochrome c oxidase subunit I"[All Fields])')
PAUSE = 0.4
TODAY = dt.date.today().isoformat()

logging.basicConfig(filename=ROOT / "logs" / "candidate_species.log", level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()


def get(url, params):
    for attempt in range(8):
        try:
            r = requests.get(url, params=params, timeout=60)
            if r.status_code == 200:
                time.sleep(PAUSE)
                return r
        except requests.RequestException:
            pass
        time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"failed: {url} {params}")


def organisms(txid):
    term = f"txid{txid}[Organism:exp] AND {COI}"
    r = get(E + "esearch.fcgi", {"db": "nuccore", "term": term, "usehistory": "y", "retmode": "json", "retmax": 0}).json()
    res = r["esearchresult"]
    n, web, key = int(res["count"]), res["webenv"], res["querykey"]
    rows = []
    for start in range(0, n, 400):
        js = get(E + "esummary.fcgi", {"db": "nuccore", "WebEnv": web, "query_key": key, "retstart": start,
                                       "retmax": 400, "retmode": "json"}).json()["result"]
        for uid in js.get("uids", []):
            d = js[uid]
            rows.append((d.get("organism", ""), d.get("taxid"), d.get("title", "")))
    log.info("txid %s: %d records, %d fetched", txid, n, len(rows))
    return pd.DataFrame(rows, columns=["organism", "taxid", "title"])


def named(sp):
    parts = sp.split()
    return (len(parts) == 2 and parts[1].islower() and not re.search(r"\b(sp|cf|aff|nr)\b\.?", sp)
            and re.fullmatch(r"[A-Z][a-z]+ [a-z\-]+", sp) is not None)


def main():
    out = []
    for fam, txid in FAMILIES.items():
        df = organisms(txid)
        df["binomial"] = df["organism"].str.split().str[:2].str.join(" ")
        g = df.groupby("binomial").agg(n_coi_records=("title", "size"), taxids=("taxid", lambda s: ";".join(
            sorted({str(x) for x in s})))).reset_index()
        g["family"] = fam
        out.append(g)
    allsp = pd.concat(out, ignore_index=True)
    allsp.to_csv(ROOT / "data" / "raw" / f"ncbi_coi_species_{TODAY}.csv", index=False)

    cand = allsp[allsp["binomial"].map(named)].copy()
    res = pd.read_csv(ROOT / "data" / "clean" / "resistance_by_species.csv")
    cand["has_resistance_record"] = cand["binomial"].isin(res["species"])
    cand = cand.rename(columns={"binomial": "species"}).sort_values(["family", "species"])
    cand.to_csv(ROOT / "data" / "clean" / "candidate_species.csv", index=False)
    log.info("named species with COI: %d; with resistance record: %d", len(cand), int(cand["has_resistance_record"].sum()))
    print(cand.groupby("family").agg(species=("species", "size"), resistant=("has_resistance_record", "sum")))


if __name__ == "__main__":
    main()
