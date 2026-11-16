"""Stage 4b: download COI, 18S and 28S records from GenBank for the analysis taxa and outgroups.

Inputs:  data/clean/analysis_taxa.csv, data/manual/outgroups.csv (optional)
Outputs: data/raw/genbank_<date>/<marker>.gb.gz   GenBank flat files as downloaded
         data/clean/sequence_records.csv           one row per record x marker with metadata and the
                                                   extracted marker sequence
         logs/sequences_fetch.log

For each taxon and marker all matching nucleotide records are listed; when more than MAX_PER are found,
the MAX_PER/2 longest and a seeded random MAX_PER/2 of the rest are downloaded. The marker region is cut
from annotated features where present (CDS for COI, rRNA for 18S and 28S); otherwise the whole record is
used when it is shorter than 3,000 bp.
"""
import datetime as dt
import gzip
import io
import logging
import random
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TODAY = sys.argv[1] if len(sys.argv) > 1 else dt.date.today().isoformat()
OUT = ROOT / "data/raw" / f"genbank_{TODAY}"
SEED, MAX_PER, PAUSE = 2026, 60, 0.4
COI_PRODUCT = re.compile(r"cytochrome( c)? oxidase subunit (i|1)\b")

TERMS = {
    "COI": '(COI[All Fields] OR COX1[All Fields] OR CO1[All Fields] OR "cytochrome oxidase subunit I"[All Fields] '
           'OR "cytochrome c oxidase subunit I"[All Fields] OR "cytochrome c oxidase subunit 1"[All Fields] '
           'OR (mitochondrion[Title] AND "complete genome"[Title]))',
    "18S": '("18S"[All Fields] OR "small subunit ribosomal"[All Fields])',
    "28S": '("28S"[All Fields] OR "large subunit ribosomal"[All Fields])',
}

logging.basicConfig(filename=ROOT / "logs/sequences_fetch.log", level=logging.INFO, filemode="w",
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()


def get(endpoint, params, post=False):
    for attempt in range(8):
        try:
            r = (requests.post if post else requests.get)(E + endpoint, data=params if post else None,
                                                          params=None if post else params, timeout=120)
            if r.status_code == 200:
                time.sleep(PAUSE)
                return r
        except requests.RequestException:
            pass
        time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"failed {endpoint}")


def marker_region(rec, marker):
    """Return (sequence string, how) for the marker inside a GenBank record, or (None, reason)."""
    def q(f, key):
        return " ".join(f.qualifiers.get(key, [])).lower()
    hits = []
    for f in rec.features:
        text = q(f, "gene") + " | " + q(f, "product")
        if marker == "COI" and f.type == "CDS":
            gene = q(f, "gene").replace(" ", "").replace("-", "")
            if gene in ("coi", "cox1", "co1", "coxi") or COI_PRODUCT.search(q(f, "product")):
                hits.append(f)
        elif marker in ("18S", "28S") and f.type in ("rRNA", "misc_RNA", "misc_feature"):
            key = ("18s", "small subunit") if marker == "18S" else ("28s", "large subunit")
            prod = q(f, "product") + " " + q(f, "note")
            other = ("28s", "large subunit", "5.8s", "internal transcribed") if marker == "18S" else \
                ("18s", "small subunit", "5.8s", "internal transcribed")
            if any(k in prod for k in key) and not any(o in prod for o in other):
                hits.append(f)
    if hits:
        f = max(hits, key=lambda x: len(x))
        try:
            return str(f.extract(rec.seq)).upper(), "feature"
        except Exception:
            return None, "feature could not be extracted"
    desc = rec.description.lower()
    if marker != "COI" and ("internal transcribed" in desc or "5.8s" in desc or " its" in desc):
        return None, "ITS-region record without a separate rRNA feature"
    if len(rec.seq) <= 3000:
        return str(rec.seq).upper(), "whole record"
    return None, "no annotated marker in a long record"


def main():
    random.seed(SEED)
    OUT.mkdir(parents=True, exist_ok=True)
    taxa = pd.read_csv(ROOT / "data/clean/analysis_taxa.csv")
    taxa = taxa[taxa["in_analysis"]][["taxon", "family", "search_names"]].assign(role="ingroup")
    og = ROOT / "data/manual/outgroups.csv"
    if og.exists():
        o = pd.read_csv(og)
        taxa = pd.concat([taxa, pd.DataFrame({"taxon": o["taxon"], "family": o["family"],
                                              "search_names": o["taxon"], "role": "outgroup"})], ignore_index=True)
    rows = []
    for marker, term in TERMS.items():
        handle = gzip.open(OUT / f"{marker}.gb.gz", "wt")
        for t in taxa.itertuples(index=False):
            ids = []
            for name in t.search_names.split("; "):
                q = f'"{name}"[Organism] AND {term} AND 150:20000[SLEN]'
                res = get("esearch.fcgi", {"db": "nuccore", "term": q, "retmax": 5000, "retmode": "json"}).json()
                ids += res["esearchresult"].get("idlist", [])
            ids = list(dict.fromkeys(ids))
            n_found = len(ids)
            if n_found > MAX_PER:
                lens = {}
                for i in range(0, len(ids), 300):
                    js = get("esummary.fcgi", {"db": "nuccore", "id": ",".join(ids[i:i + 300]), "retmode": "json"},
                             post=True).json()["result"]
                    lens.update({u: int(js[u].get("slen", 0)) for u in js.get("uids", [])})
                ranked = sorted(ids, key=lambda u: -lens.get(u, 0))
                long_ids = ranked[:MAX_PER // 2]
                ids = long_ids + random.sample(ranked[MAX_PER // 2:], MAX_PER // 2)
            log.info("%s %s: found %d, fetching %d", marker, t.taxon, n_found, len(ids))
            if not ids:
                continue
            txt = get("efetch.fcgi", {"db": "nuccore", "id": ",".join(ids), "rettype": "gbwithparts",
                                      "retmode": "text"}, post=True).text
            handle.write(txt)
            for rec in SeqIO.parse(io.StringIO(txt), "genbank"):
                src = next((f for f in rec.features if f.type == "source"), None)
                qual = (lambda k: "; ".join(src.qualifiers.get(k, [])) if src else "")
                seq, how = marker_region(rec, marker)
                rows.append(dict(
                    taxon=t.taxon, family=t.family, role=t.role, marker=marker, accession=rec.id,
                    organism=rec.annotations.get("organism", ""), record_length=len(rec.seq),
                    title=rec.description, voucher=qual("specimen_voucher"), isolate=qual("isolate"),
                    strain=qual("strain"), country=qual("geo_loc_name") or qual("country"), host=qual("host"),
                    collection_date=qual("collection_date"), n_records_found=n_found,
                    extraction=how, marker_length=len(seq) if seq else 0, sequence=seq or ""))
        handle.close()
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "data/clean/sequence_records.csv", index=False)
    log.info("records %d", len(df))
    print(df.groupby("marker").agg(records=("accession", "size"), taxa=("taxon", "nunique")))


if __name__ == "__main__":
    main()
