"""Summarise the paper-level check of the stage 2 verification sample.

Input:  data/manual/verification_sample_stage2.csv (matches_source filled by V. K. Chaudhari)
Output: results/tables/stage2_verification_summary.csv
Grey-literature "partly" means only a related or citing source was found, not
the original record; it is reported separately from journal-backed cases.
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
# cases checked against the original paper: all journal-backed cases plus grey cases whose original was found
original_checked = int((~grey).sum()) + int(((grey) & (verdict == "yes")).sum())
g = pd.DataFrame({"grey": grey, "verdict": verdict})
grey_counts = g[g.grey].verdict.value_counts()
journal_counts = g[~g.grey].verdict.value_counts()
rows = [
    ("cases in sample", len(s)),
    ("yes", int(counts.get("yes", 0))),
    ("partly", int(counts.get("partly", 0))),
    ("no", int(counts.get("no", 0))),
    ("not traceable", int(counts.get("not traceable", 0))),
    ("not checked", int(counts.get("not checked", 0))),
    ("journal-backed cases", int((~grey).sum())),
    ("journal-backed: yes", int(journal_counts.get("yes", 0))),
    ("journal-backed: partly", int(journal_counts.get("partly", 0))),
    ("journal-backed: no", int(journal_counts.get("no", 0))),
    ("grey-literature cases", int(grey.sum())),
    ("grey: yes (original found)", int(grey_counts.get("yes", 0))),
    ("grey: partly (related or citing source only)", int(grey_counts.get("partly", 0))),
    ("grey: not traceable", int(grey_counts.get("not traceable", 0))),
    ("cases checked against the original source", original_checked),
    ("share of those marked no (%)", round(100 * counts.get("no", 0) / original_checked, 1) if original_checked else None),
]
out = ROOT / "results" / "tables" / "stage2_verification_summary.csv"
pd.DataFrame(rows, columns=["item", "value"]).to_csv(out, index=False)
print(pd.DataFrame(rows, columns=["item", "value"]).to_string(index=False))
