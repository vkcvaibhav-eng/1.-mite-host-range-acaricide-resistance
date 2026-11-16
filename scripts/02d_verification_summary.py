"""Summarise the paper-level check of the stage 2 verification sample.

Input:  data/manual/verification_sample_stage2.csv (matches_source filled by V. K. Chaudhari)
Output: results/tables/stage2_verification_summary.csv
Rule: a blank matches_source on a grey-literature case (FAO report, personal
communication, short citation) counts as "not traceable".
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
s = pd.read_csv(ROOT / "data" / "manual" / "verification_sample_stage2.csv", dtype=str, keep_default_na=False)

grey = s["claude_note"].str.startswith("grey literature")
verdict = s["matches_source"].where(s["matches_source"] != "", "not traceable").where(
    ~((s["matches_source"] == "") & ~grey), "not checked")

counts = verdict.value_counts()
traceable = int(counts.get("yes", 0) + counts.get("partly", 0) + counts.get("no", 0))
rows = [
    ("cases in sample", len(s)),
    ("yes", int(counts.get("yes", 0))),
    ("partly", int(counts.get("partly", 0))),
    ("no", int(counts.get("no", 0))),
    ("not traceable", int(counts.get("not traceable", 0))),
    ("not checked", int(counts.get("not checked", 0))),
    ("traceable cases (yes + partly + no)", traceable),
    ("share of traceable cases marked no (%)", round(100 * counts.get("no", 0) / traceable, 1) if traceable else None),
]
out = ROOT / "results" / "tables" / "stage2_verification_summary.csv"
pd.DataFrame(rows, columns=["item", "value"]).to_csv(out, index=False)
print(pd.DataFrame(rows, columns=["item", "value"]).to_string(index=False))
