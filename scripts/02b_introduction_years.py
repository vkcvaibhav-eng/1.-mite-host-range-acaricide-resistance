"""Stage 2b: first introduction year per active ingredient, and years from introduction to first resistance case.

Inputs
  data/raw/ppdb_introduction_dates_2026-11-13.txt   'Introduction & key dates' text from PPDB/BPDB/VSDB
  data/manual/irac_map.csv                          active-ingredient list
  data/manual/intro_year_overrides.csv              second-source years where PPDB is missing or contradicted
  data/clean/resistance_cases.csv                   from scripts/02_resistance.py

Outputs
  data/clean/ai_introduction_years.csv
  data/clean/resistance_time_to_first_case.csv     one row per species x active ingredient with a main-analysis case
  logs/introduction_years.log

Rule: intro_year is the earliest year attached to 'introduced', 'registered', 'registration', 'marketed',
'reported' or 'approval'. Years attached only to discovery, isolation, synthesis or description are kept in
intro_year_fallback and flagged. Decades ('1960s', 'Mid-1990s') take the first year of the decade.
"""

import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
KEEP = re.compile(r"introduc|intoduc|register|registration|market|reported|approval", re.I)
EARLY = re.compile(r"discover|isolat|synthes|describ", re.I)


def years(entry: str) -> list[int]:
    return [int(y) for y in re.findall(r"(1[89]\d\d|20\d\d)", entry)]


def parse(text: str) -> tuple[float, float, str]:
    strict, fallback = [], []
    for entry in text.split(";"):
        ys = years(entry)
        if not ys:
            continue
        if KEEP.search(entry):
            strict += ys
        elif EARLY.search(entry):
            fallback += ys
    s = min(strict) if strict else float("nan")
    f = min(fallback) if fallback else float("nan")
    basis = "introduction/registration" if strict else ("discovery only" if fallback else "no year")
    return s, f, basis


def main() -> None:
    log = open(ROOT / "logs/introduction_years.log", "w", encoding="utf-8")
    raw = pd.read_csv(ROOT / "data/raw/ppdb_introduction_dates_2026-11-13.txt", sep="|", header=None,
                      names=["ai_id", "ppdb_name", "ppdb_record", "introduction_text", "availability"], dtype=str)
    parsed = raw["introduction_text"].map(parse)
    raw["intro_year"] = [p[0] for p in parsed]
    raw["intro_year_fallback"] = [p[1] for p in parsed]
    raw["intro_basis"] = [p[2] for p in parsed]

    irac = pd.read_csv(ROOT / "data/manual/irac_map.csv", dtype={"ai_id": str})
    ai = irac[["ai_id", "ai_name", "irac_group"]].merge(raw, on="ai_id", how="left")
    ai["intro_basis"] = ai["intro_basis"].fillna("not in PPDB")
    ai.loc[ai["irac_group"].isin(["MIX"]) | ai["ai_name"].str.contains("Unspecified"), "intro_basis"] = \
        "not a single compound"
    ai["intro_year_used"] = ai["intro_year"].fillna(ai["intro_year_fallback"])
    # second-source corrections where the database year is missing or contradicted by the case records
    ov = pd.read_csv(ROOT / "data/manual/intro_year_overrides.csv", dtype={"ai_id": str})
    ai = ai.merge(ov[["ai_id", "intro_year", "source", "confidence"]].rename(columns={
        "intro_year": "override_year", "source": "override_source", "confidence": "override_confidence"}),
        on="ai_id", how="left")
    has = ai["override_year"].notna()
    ai.loc[has, "intro_year_used"] = ai.loc[has, "override_year"]
    ai.loc[has, "intro_basis"] = "second source (data/manual/intro_year_overrides.csv)"
    ai.to_csv(ROOT / "data/clean/ai_introduction_years.csv", index=False)
    print(ai["intro_basis"].value_counts().to_string(), file=log)

    cases = pd.read_csv(ROOT / "data/clean/resistance_cases.csv", dtype={"ai_id": str})
    cases = cases[cases["include_main"]]
    pairs = (cases.groupby(["species", "family", "ai_id", "ai_name", "irac_group"])
             .agg(first_report_year=("year_of_report", "min"), n_cases=("case_id", "size")).reset_index()
             .merge(ai[["ai_id", "intro_year_used", "intro_basis"]], on="ai_id", how="left"))
    pairs["years_to_first_case"] = pairs["first_report_year"] - pairs["intro_year_used"]
    pairs["conflict"] = pairs["years_to_first_case"] < 0
    pairs.to_csv(ROOT / "data/clean/resistance_time_to_first_case.csv", index=False)

    conf = pairs[pairs["conflict"]]
    ai_conf = conf.groupby(["ai_name", "intro_year_used"]).agg(
        n_pairs=("species", "size"), earliest_case=("first_report_year", "min")).reset_index()
    print("\nfirst case earlier than introduction year (source year probably wrong):", file=log)
    print(ai_conf.to_string(index=False), file=log)
    print(f"\npairs: {len(pairs)}; with a usable year: {pairs['intro_year_used'].notna().sum()}; "
          f"conflicts: {len(conf)}", file=log)
    log.close()
    print(open(ROOT / "logs/introduction_years.log").read())


if __name__ == "__main__":
    main()
