"""Stage 2: case-level resistance dataset and per-species summary.

Inputs
  data/raw/aprd_cases_part1_2026-11-13.txt      APRD cases (space-separated, as transferred)
  data/raw/aprd_cases_part2_2026-11-13.txt      APRD cases (pipe-separated)
  data/raw/aprd_references_2026-11-13.txt       reference list, index|citation
  data/raw/aprd_acari_ag_species_2026-11-13.csv species table (arId, names)
  data/raw/aprd_species_ai_cases_2026-11-13.csv active-ingredient names
  data/raw/literature_resistant_species_2026-11-13.csv
  data/manual/name_map_stage1.csv               name decisions
  data/manual/irac_map.csv                      IRAC group per active ingredient
  data/manual/literature_cases.csv              cases from papers for species absent from APRD

Outputs
  data/clean/resistance_cases.csv               one row per APRD case, with inclusion flags
  data/clean/resistance_by_species.csv          one row per standardised species
  data/manual/verification_sample_stage2.csv    random 10% of included cases for hand checking
  results/tables/stage2_record_flow.csv         counts for the record-flow diagram
  logs/resistance.log

Run: python scripts/02_resistance.py
"""

from __future__ import annotations

import datetime as dt
import logging
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw"
MANUAL = ROOT / "data/manual"
CLEAN = ROOT / "data/clean"
TABLES = ROOT / "results/tables"
LOG = ROOT / "logs/resistance.log"

EXPECTED_ROWS = 1027          # count reported by the browser export
EXPECTED_ID_SUM = 7748965     # sum of numeric case IDs reported by the browser export
SEED = 2026

TYPE_CODES = {
    "F": "field-evolved",
    "U": "unspecified",
    "L": "lab-selected only",
    "LF": "lab-selected after field-evolved resistance",
    "L+F": "lab-selected and field-evolved",
    "L+LF": "lab-selected and lab-selected after field resistance",
}
PART1_TYPES = {
    "F": "F", "U": "U", "Lab selected*": "L",
    "Lab selected after field-evolved resistance**": "LF",
    "Lab selected*Field Evolved Resistance": "L+F",
    "Lab selected*Lab selected after field-evolved resistance**": "L+LF",
}
NOT_GROUPS = {"UN", "NC", "MIX"}  # not counted as IRAC mode-of-action groups

log = logging.getLogger("resistance")


def setup_logging() -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(LOG, mode="w", encoding="utf-8")
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    log.addHandler(fh)
    log.addHandler(logging.StreamHandler(sys.stdout))
    log.setLevel(logging.INFO)
    log.info("run started %s; pandas %s; numpy %s", dt.datetime.now().isoformat(timespec="seconds"),
             pd.__version__, np.__version__)


def binomial(name: str) -> str:
    name = " ".join(name.strip().split())
    return name[:1].upper() + name[1:].lower()


def read_cases() -> pd.DataFrame:
    pat = re.compile(r"^(\d+) (\d+) (G\d+) (\d{4}) (" + "|".join(
        re.escape(k) for k in sorted(PART1_TYPES, key=len, reverse=True)) + r") (.*) (\d+)$")
    rows = []
    for line in (RAW / "aprd_cases_part1_2026-11-13.txt").read_text(encoding="utf-8").splitlines():
        m = pat.match(line)
        if not m:
            raise ValueError(f"unparsed line in part 1: {line}")
        g = list(m.groups())
        g[4] = PART1_TYPES[g[4]]
        rows.append(g)
    for line in (RAW / "aprd_cases_part2_2026-11-13.txt").read_text(encoding="utf-8").splitlines():
        p = line.split("|")
        if len(p) != 7:
            raise ValueError(f"unparsed line in part 2: {line}")
        rows.append(p)
    df = pd.DataFrame(rows, columns=["arId", "ai_id", "case_id", "year_of_report", "type_code",
                                     "location", "ref_idx"])
    id_sum = df["case_id"].str[1:].astype(int).sum()
    if len(df) != EXPECTED_ROWS or id_sum != EXPECTED_ID_SUM or df["case_id"].duplicated().any():
        raise ValueError(f"checksum failed: rows {len(df)}, id sum {id_sum}")
    log.info("cases read: %d; checksum matches the browser export", len(df))
    df["year_of_report"] = df["year_of_report"].astype(int)
    df["ref_idx"] = df["ref_idx"].astype(int)
    df["location"] = df["location"].str.replace(r"\s*--\s*", " -- ", regex=True).str.strip(" -")
    df["country"] = df["location"].str.split(" -- ").str[0].str.strip()
    return df


def ref_key(citation: str) -> str:
    """Lower-case title with punctuation removed, used to spot one paper entered twice."""
    m = re.search(r"\(\d{4}\)\.\s*(.+?)\.\s", citation + " ")
    title = m.group(1) if m else citation
    first_author = re.split(r"[ ,.]", citation.strip(), maxsplit=1)[0].lower()
    return first_author + ":" + re.sub(r"[^a-z0-9]", "", title.lower())[:80]


def main() -> None:
    setup_logging()
    cases = read_cases()

    refs = pd.DataFrame([l.split("|", 1) for l in
                         (RAW / "aprd_references_2026-11-13.txt").read_text(encoding="utf-8").splitlines()],
                        columns=["ref_idx", "reference"])
    refs["ref_idx"] = refs["ref_idx"].astype(int)
    refs["ref_key"] = refs["reference"].map(ref_key)
    cases = cases.merge(refs, on="ref_idx", how="left")

    sp = pd.read_csv(RAW / "aprd_acari_ag_species_2026-11-13.csv", dtype={"arId": str})
    names = pd.read_csv(MANUAL / "name_map_stage1.csv", dtype=str).fillna("")
    name_map = dict(zip(names["name_as_listed"].str.lower(), names["query_name"]))
    sp["species_as_listed"] = sp["species_as_listed"].str.lower()
    sp["species"] = sp["species_as_listed"].map(lambda s: name_map.get(s, binomial(s)))
    sp["family"] = sp["family_as_listed"].str.capitalize()
    cases = cases.merge(sp[["arId", "species_as_listed", "species", "family"]], on="arId", how="left")
    if cases["species"].isna().any():
        raise ValueError("cases with unknown arId")

    ai = (pd.read_csv(RAW / "aprd_species_ai_cases_2026-11-13.csv", dtype={"ai_id": str})
          [["ai_id", "ai_name"]].drop_duplicates("ai_id"))
    irac = pd.read_csv(MANUAL / "irac_map.csv", dtype={"ai_id": str})
    cases = cases.merge(ai, on="ai_id", how="left").merge(
        irac[["ai_id", "irac_group", "chemical_class"]], on="ai_id", how="left")
    missing = cases.loc[cases["irac_group"].isna(), "ai_id"].unique()
    if len(missing):
        raise ValueError(f"active ingredients without IRAC mapping: {missing}")

    cases["resistance_type"] = cases["type_code"].map(TYPE_CODES)

    # possible duplicates: same species, compound, country and paper (paper entered twice in the database)
    key = ["species", "ai_id", "country", "year_of_report", "ref_key"]
    grp = cases.groupby(key)["ref_idx"].transform("nunique")
    cases["possible_duplicate"] = grp > 1

    # inclusion for the main analysis
    cases["include_main"] = True
    cases["exclusion_reason"] = ""
    lab = cases["type_code"].eq("L")
    cases.loc[lab, ["include_main", "exclusion_reason"]] = [False, "lab-selected population only"]
    # Where one paper is entered twice in the database, keep the entry with more location detail
    # (ties: keep the lower reference index) and exclude the other.
    for _, g in cases.loc[cases["possible_duplicate"] & cases["include_main"]].groupby(key):
        if g["ref_idx"].nunique() > 1:
            detail = g.groupby("ref_idx")["location"].apply(lambda s: s.str.count(" -- ").mean())
            keep_ref = detail.sort_values(ascending=False, kind="stable").index[0] \
                if detail.nunique() > 1 else detail.index.min()
            idx = g.index[g["ref_idx"] != keep_ref]
            cases.loc[idx, ["include_main", "exclusion_reason"]] = [False, "same paper entered twice in database"]

    cols = ["case_id", "species", "species_as_listed", "family", "arId", "ai_id", "ai_name", "irac_group",
            "chemical_class", "year_of_report", "resistance_type", "country", "location", "ref_idx",
            "reference", "possible_duplicate", "include_main", "exclusion_reason"]
    cases = cases[cols].sort_values(["family", "species", "ai_name", "year_of_report", "case_id"])
    cases.insert(len(cols), "source", "Arthropod Pesticide Resistance Database")
    CLEAN.mkdir(parents=True, exist_ok=True)
    cases.to_csv(CLEAN / "resistance_cases.csv", index=False)

    # literature cases for species absent from the database (logs/literature_search_stage2.md)
    lit_cases = pd.read_csv(MANUAL / "literature_cases.csv", dtype=str)
    ai_by_name = irac.set_index(irac["ai_name"].str.lower())
    lit_cases["ai_id"] = lit_cases["ai_name"].str.lower().map(ai_by_name["ai_id"])
    lit_cases["irac_group"] = lit_cases["ai_name"].str.lower().map(ai_by_name["irac_group"])
    lit_cases["chemical_class"] = lit_cases["ai_name"].str.lower().map(ai_by_name["chemical_class"])
    if lit_cases["ai_id"].isna().any():
        raise ValueError("literature compound without IRAC mapping")
    lit_cases = lit_cases.assign(
        species_as_listed=lit_cases["species"], family=lit_cases["family"].str.capitalize(), arId="",
        year_of_report=lit_cases["year_of_report"].astype(int),
        country=lit_cases["country"], ref_idx=-1,
        possible_duplicate=False,
        include_main=lit_cases["meets_case_definition"].eq("yes"),
        exclusion_reason=lit_cases["meets_case_definition"].map(
            {"yes": "", "no": "below case definition", "pending": "resistance ratio not yet verified"}),
        source="literature")
    cases = pd.concat([cases, lit_cases[cases.columns.drop([]).intersection(lit_cases.columns)]],
                      ignore_index=True)
    cases.to_csv(CLEAN / "resistance_cases.csv", index=False)
    log.info("literature cases added: %d (%d meet the case definition)", len(lit_cases),
             int(lit_cases["include_main"].sum()))

    main_cases = cases[cases["include_main"]]

    def groups(s: pd.Series) -> int:
        return s[~s.isin(NOT_GROUPS)].nunique()

    by_sp = (main_cases.groupby(["species", "family"])
             .agg(n_cases=("case_id", "size"),
                  n_active_ingredients=("ai_id", "nunique"),
                  n_irac_groups=("irac_group", groups),
                  irac_groups=("irac_group", lambda s: "; ".join(sorted(set(s) - NOT_GROUPS))),
                  n_chemical_classes=("chemical_class", "nunique"),
                  n_countries=("country", "nunique"),
                  first_report_year=("year_of_report", "min"),
                  last_report_year=("year_of_report", "max"),
                  n_references=("ref_idx", "nunique"),
                  n_field_evolved=("resistance_type", lambda s: int((s == "field-evolved").sum())))
             .reset_index())
    excluded_only = sorted(set(cases["species"]) - set(by_sp["species"]))
    if excluded_only:
        log.warning("species whose every case is excluded: %s", excluded_only)

    src = cases.groupby("species")["source"].agg(lambda s: "; ".join(sorted(set(s))))
    by_sp["source_note"] = by_sp["species"].map(src)
    lit = pd.DataFrame(columns=["species"])  # all literature species now enter through literature_cases.csv
    for c in ["n_cases", "n_active_ingredients", "n_irac_groups", "n_chemical_classes", "n_countries",
              "n_references", "n_field_evolved"]:
        by_sp[c] = by_sp[c].fillna(0).astype(int)
    by_sp = by_sp.sort_values(["family", "species"])
    by_sp.to_csv(CLEAN / "resistance_by_species.csv", index=False)

    flow = pd.DataFrame([
        ("APRD Acari cases in the three families", int((cases["source"] != "literature").sum())),
        ("excluded: lab-selected population only", int((cases["exclusion_reason"] == "lab-selected population only").sum())),
        ("excluded: same paper entered twice", int((cases["exclusion_reason"] == "same paper entered twice in database").sum())),
        ("cases in main analysis", len(main_cases)),
        ("species with at least one main-analysis case", main_cases["species"].nunique()),
        ("literature cases screened", int((cases["source"] == "literature").sum())),
        ("literature cases meeting the case definition", int(((cases["source"] == "literature") & cases["include_main"]).sum())),
        ("APRD references cited by main-analysis cases", main_cases.loc[main_cases["ref_idx"] >= 0, "ref_idx"].nunique()),
        ("main-analysis cases from FAO or personal communications", int(main_cases["reference"].str.contains(
            r"FAO|Personal communication|Communication to", case=False, regex=True).sum())),
    ], columns=["step", "n"])
    TABLES.mkdir(parents=True, exist_ok=True)
    flow.to_csv(TABLES / "stage2_record_flow.csv", index=False)

    sample = main_cases.sample(frac=0.10, random_state=SEED).sort_values(["species", "ai_name"])
    sample = sample[["case_id", "species", "ai_name", "year_of_report", "resistance_type", "location",
                     "reference"]].assign(matches_source="", note="")
    # The sample file holds hand-entered checks: write it only once, then only compare.
    sample_path = MANUAL / "verification_sample_stage2.csv"
    if sample_path.exists():
        kept = pd.read_csv(sample_path, dtype=str)
        if sorted(kept["case_id"]) != sorted(sample["case_id"].astype(str)):
            log.warning("verification sample differs from the cases now drawn; hand checks kept, sample not rewritten")
    else:
        sample.to_csv(sample_path, index=False)

    log.info("record flow\n%s", flow.to_string(index=False))
    log.info("verification sample: %d cases (seed %d)", len(sample), SEED)
    log.info("run finished %s", dt.datetime.now().isoformat(timespec="seconds"))


if __name__ == "__main__":
    main()
