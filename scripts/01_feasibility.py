"""Stage 1 feasibility: resistant phytophagous mite species and their GenBank sequence coverage.

Inputs
  data/raw/aprd_acari_ag_species_2026-11-13.csv      full APRD species table for order Acari, group AG
  data/raw/aprd_species_ai_cases_2026-11-13.csv      cases per species and active ingredient (APRD species pages)
  data/raw/bold_summary_2026-11-13.csv               BOLD specimen and COI counts per species
  data/raw/literature_resistant_species_2026-11-13.csv  resistant species found in papers, not in APRD
  data/manual/name_map_stage1.csv                    name decisions (synonyms, spellings, flags)
  (data/raw/aprd_ai_pages_acari_2026-11-13.csv, the first 48-page transcription, is superseded and not read)

Outputs
  data/clean/stage1_species_list.csv     one row per standardised species name
  data/clean/feasibility_counts.csv      species list plus NCBI taxonomy check and record counts per marker
  results/tables/stage1_summary.csv      counts by family, used in the checkpoint report
  logs/feasibility_excluded.csv          every APRD row left out, with the reason
  logs/feasibility_queries.tsv           every NCBI query and its count
  logs/feasibility.log                   run information and zero-hit names

Run: python scripts/01_feasibility.py
"""

from __future__ import annotations

import datetime as dt
import logging
import platform
import sys
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW_APRD = ROOT / "data/raw/aprd_acari_ag_species_2026-11-13.csv"
RAW_APRD_AI = ROOT / "data/raw/aprd_species_ai_cases_2026-11-13.csv"
RAW_BOLD = ROOT / "data/raw/bold_summary_2026-11-13.csv"
RAW_LIT = ROOT / "data/raw/literature_resistant_species_2026-11-13.csv"
NAME_MAP = ROOT / "data/manual/name_map_stage1.csv"
OUT_SPECIES = ROOT / "data/clean/stage1_species_list.csv"
OUT_COUNTS = ROOT / "data/clean/feasibility_counts.csv"
OUT_SUMMARY = ROOT / "results/tables/stage1_summary.csv"
LOG_EXCLUDED = ROOT / "logs/feasibility_excluded.csv"
LOG_QUERIES = ROOT / "logs/feasibility_queries.tsv"
LOG_RUN = ROOT / "logs/feasibility.log"

IN_SCOPE_FAMILIES = {"tetranychidae", "tenuipalpidae", "eriophyidae"}

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL = "mite-host-range-acaricide-resistance"
PAUSE = 1.0  # NCBI allows 3 requests per second without a key; slower to avoid shared-IP throttling

# Length window matches the plan's family-level counts (150-3,000 bp).
LEN = "150:3000[SLEN]"
MARKERS = {
    "COI": (
        '(COI[Gene] OR COX1[Gene] OR CO1[Gene] OR COXI[Gene] '
        'OR "cytochrome oxidase subunit I"[Title] OR "cytochrome oxidase subunit 1"[Title] '
        'OR "cytochrome c oxidase subunit I"[Title] OR "cytochrome c oxidase subunit 1"[Title]) AND ' + LEN
    ),
    "18S": '("18S ribosomal RNA"[Title] OR "18S rRNA"[Title] OR "small subunit ribosomal RNA"[Title]) AND ' + LEN,
    "28S": '("28S ribosomal RNA"[Title] OR "28S rRNA"[Title] OR "large subunit ribosomal RNA"[Title]) AND ' + LEN,
    "ITS2": '("internal transcribed spacer 2"[Title] OR ITS2[Title]) AND ' + LEN,
    "mitogenome": '"complete genome"[Title] AND mitochondrion[filter] AND 10000:30000[SLEN]',
}

log = logging.getLogger("feasibility")


def setup_logging() -> None:
    LOG_RUN.parent.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(LOG_RUN, mode="w", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    log.addHandler(handler)
    log.addHandler(logging.StreamHandler(sys.stdout))
    log.setLevel(logging.INFO)
    log.info("run started %s", dt.datetime.now().isoformat(timespec="seconds"))
    log.info("python %s, pandas %s, requests %s, %s", platform.python_version(), pd.__version__,
             requests.__version__, platform.platform())


def binomial(name: str) -> str:
    name = " ".join(name.strip().split())
    return name[:1].upper() + name[1:].lower()


# ---------------------------------------------------------------- species list

def build_species_list() -> pd.DataFrame:
    aprd = pd.read_csv(RAW_APRD, dtype=str).fillna("")
    for col in ("species_as_listed", "family_as_listed"):
        aprd[col] = aprd[col].str.strip().str.lower()

    keep = aprd["family_as_listed"].isin(IN_SCOPE_FAMILIES)
    excluded = aprd.loc[~keep].copy()
    excluded["reason"] = "family outside Tetranychidae, Tenuipalpidae, Eriophyidae"
    hal = excluded["species_as_listed"].eq("halotydeus destructor")
    excluded.loc[hal, "reason"] = ("Penthaleidae (database lists family as pentatomidae); "
                                   "phytophagous but outside the three study families")
    LOG_EXCLUDED.parent.mkdir(parents=True, exist_ok=True)
    excluded.to_csv(LOG_EXCLUDED, index=False)
    log.info("APRD rows read: %d; kept: %d; excluded: %d (%d species names)", len(aprd), keep.sum(),
             len(excluded), excluded["species_as_listed"].nunique())

    aprd = aprd.loc[keep]
    names = pd.read_csv(NAME_MAP, dtype=str).fillna("")
    name_map = dict(zip(names["name_as_listed"].str.lower(), names["query_name"]))
    flags = dict(zip(names["name_as_listed"].str.lower(), names["confidence"]))

    aprd = aprd.assign(query_name=aprd["species_as_listed"].map(lambda s: name_map.get(s, binomial(s))))
    ai = pd.read_csv(RAW_APRD_AI, dtype={"arId": str, "ai_id": str})
    ai = ai.merge(aprd[["arId", "query_name"]], on="arId", how="inner")
    ai_per = (ai.groupby("query_name")
              .agg(aprd_active_ingredients=("ai_id", "nunique"),
                   aprd_compounds=("ai_name", lambda s: "; ".join(sorted(set(s)))))
              .reset_index())
    aprd["n_cases"] = aprd["n_cases"].astype(int)
    aprd["n_citations"] = pd.to_numeric(aprd["n_citations"], errors="coerce").fillna(0).astype(int)
    per_species = (aprd.groupby("query_name")
                   .agg(family=("family_as_listed", "first"),
                        names_as_listed=("species_as_listed", lambda s: "; ".join(sorted(set(s)))),
                        aprd_ids=("arId", lambda s: "; ".join(sorted(set(s), key=int))),
                        aprd_cases=("n_cases", "sum"),
                        aprd_citations=("n_citations", "sum"),
                        aprd_countries=("countries", lambda s: len({c for x in s for c in x.split("; ") if c})))
                   .reset_index()
                   .merge(ai_per, on="query_name", how="left"))
    per_species["source"] = "APRD"

    lit = pd.read_csv(RAW_LIT, dtype=str).fillna("")
    lit = lit.assign(query_name=lit["species"].map(binomial), names_as_listed=lit["species"],
                     aprd_cases=0, aprd_active_ingredients=0, source="literature")
    lit = lit.loc[~lit["query_name"].isin(per_species["query_name"]),
                  ["query_name", "family", "names_as_listed", "aprd_cases", "aprd_active_ingredients", "source"]]
    log.info("literature species already in APRD and dropped from the literature list: %s",
             sorted(set(pd.read_csv(RAW_LIT)["species"].map(binomial)) & set(per_species["query_name"])))

    species = pd.concat([per_species, lit], ignore_index=True)
    species["family"] = species["family"].str.capitalize()
    species["name_flag"] = species["names_as_listed"].map(
        lambda s: "; ".join(sorted({flags[n] for n in s.split("; ") if flags.get(n) in ("check", "decision")})))
    species = species.sort_values(["family", "query_name"]).reset_index(drop=True)
    OUT_SPECIES.parent.mkdir(parents=True, exist_ok=True)
    species.to_csv(OUT_SPECIES, index=False)
    log.info("standardised species: %d (APRD %d, literature %d)", len(species),
             (species["source"] == "APRD").sum(), (species["source"] == "literature").sum())
    return species


# ---------------------------------------------------------------- NCBI queries

class Entrez:
    def __init__(self) -> None:
        self.session = requests.Session()
        self.rows: list[dict] = []

    def _get(self, endpoint: str, params: dict) -> dict:
        params = {**params, "retmode": "json", "tool": TOOL}
        for attempt in range(8):
            time.sleep(PAUSE)
            r = self.session.get(EUTILS + endpoint, params=params, timeout=60)
            if r.status_code == 200:
                return r.json()
            log.warning("HTTP %s on %s (attempt %d)", r.status_code, endpoint, attempt + 1)
            time.sleep(5 * (attempt + 1))
        r.raise_for_status()
        return {}

    def count(self, db: str, term: str, label: str, species: str) -> int:
        res = self._get("esearch.fcgi", {"db": db, "term": term, "retmax": 0})
        n = int(res["esearchresult"]["count"])
        self.rows.append({"species": species, "label": label, "db": db, "term": term, "count": n})
        return n

    def taxon(self, name: str) -> tuple[str, str, str]:
        """Return (taxid, NCBI scientific name, match type)."""
        for field, match in (("Scientific Name", "exact"), ("All Names", "synonym or other name")):
            res = self._get("esearch.fcgi", {"db": "taxonomy", "term": f'"{name}"[{field}]', "retmax": 5})
            ids = res["esearchresult"]["idlist"]
            self.rows.append({"species": name, "label": f"taxonomy {field}", "db": "taxonomy",
                              "term": f'"{name}"[{field}]', "count": len(ids)})
            if ids:
                summ = self._get("esummary.fcgi", {"db": "taxonomy", "id": ids[0]})
                sci = summ["result"][ids[0]].get("scientificname", "")
                return ids[0], sci, match
        return "", "", "not in NCBI taxonomy"


def count_sequences(species: pd.DataFrame) -> pd.DataFrame:
    ez = Entrez()
    records = []
    for name in species["query_name"]:
        taxid, sci, match = ez.taxon(name)
        row = {"query_name": name, "ncbi_taxid": taxid, "ncbi_name": sci, "ncbi_match": match}
        organism = f"txid{taxid}[Organism:exp]" if taxid else f'"{name}"[Organism]'
        row["nuccore_all"] = ez.count("nuccore", organism, "all nucleotide", name)
        for marker, term in MARKERS.items():
            row[marker] = ez.count("nuccore", f"{organism} AND {term}", marker, name)
        records.append(row)
        log.info("%-32s taxid=%-8s COI=%-5d 18S=%-4d 28S=%-4d ITS2=%-4d mito=%d", name, taxid or "-",
                 row["COI"], row["18S"], row["28S"], row["ITS2"], row["mitogenome"])
        if not taxid:
            log.warning("zero taxonomy hits: %s", name)
        elif sci != name:
            log.warning("NCBI uses a different name: %s -> %s", name, sci)
        if row["nuccore_all"] == 0:
            log.warning("zero nucleotide records: %s", name)

    pd.DataFrame(ez.rows).to_csv(LOG_QUERIES, sep="\t", index=False)
    counts = species.merge(pd.DataFrame(records), on="query_name", how="left")
    bold = pd.read_csv(RAW_BOLD).rename(columns={"species": "query_name"})
    counts = counts.merge(bold, on="query_name", how="left")
    missing = counts.loc[counts["bold_specimens"].isna(), "query_name"].tolist()
    if missing:
        log.warning("no BOLD row for: %s", missing)
    counts["has_COI"] = (counts["COI"] > 0) | (counts["bold_COI_all_regions"].fillna(0) > 0)
    counts["has_any_marker"] = (counts[["COI", "18S", "28S", "ITS2"]].sum(axis=1) > 0) | counts["has_COI"]
    counts["COI_plus_rDNA"] = counts["has_COI"] & ((counts["18S"] > 0) | (counts["28S"] > 0))
    counts["accessed"] = dt.date.today().isoformat()
    counts.to_csv(OUT_COUNTS, index=False)
    return counts


def summarise(counts: pd.DataFrame) -> pd.DataFrame:
    def block(df: pd.DataFrame, label: str) -> dict:
        return {"group": label,
                "species": len(df),
                "from_APRD": int((df["source"] == "APRD").sum()),
                "from_literature": int((df["source"] == "literature").sum()),
                "with_any_marker": int(df["has_any_marker"].sum()),
                "APRD_cases": int(df["aprd_cases"].fillna(0).sum()),
                "with_2plus_APRD_cases": int((df["aprd_cases"].fillna(0) >= 2).sum()),
                "with_COI_GenBank_or_BOLD": int(df["has_COI"].sum()),
                "with_COI_GenBank": int((df["COI"] > 0).sum()),
                "with_COI_plus_18S_or_28S": int(df["COI_plus_rDNA"].sum()),
                "with_mitogenome": int((df["mitogenome"] > 0).sum()),
                "name_needs_check": int(df["name_flag"].str.contains("check").sum())}

    rows = [block(g, fam) for fam, g in counts.groupby("family")]
    rows.append(block(counts, "All three families"))
    merged = counts.loc[counts["query_name"] != "Tetranychus cinnabarinus"]
    rows.append(block(merged, "All, T. cinnabarinus merged with T. urticae"))
    summary = pd.DataFrame(rows)
    OUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUT_SUMMARY, index=False)
    return summary


def main() -> None:
    setup_logging()
    species = build_species_list()
    counts = count_sequences(species)
    summary = summarise(counts)
    log.info("summary\n%s", summary.to_string(index=False))
    log.info("run finished %s", dt.datetime.now().isoformat(timespec="seconds"))


if __name__ == "__main__":
    main()
