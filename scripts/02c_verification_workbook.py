"""Build the Excel workbook the user fills in for the paper-level check of the stage 2 sample.

Input:  data/manual/verification_sample_stage2.csv (with Claude's first-pass columns)
        data/clean/resistance_cases.csv (for the extra Pan et al. 2020 cases)
        data/manual/verification_dois.json (DOIs found on Crossref; optional)
        data/manual/verification_papers_stage2.csv (paper list; optional)
Output: data/manual/verification_sample_stage2_for_checking.xlsx
"""
import json
from pathlib import Path

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.worksheet.datavalidation import DataValidation

ROOT = Path(__file__).resolve().parents[1]
MANUAL = ROOT / "data" / "manual"
OUT = MANUAL / "verification_sample_stage2_for_checking.xlsx"

FONT = "Arial"
YELLOW = PatternFill("solid", start_color="FFF2CC")
HEAD = PatternFill("solid", start_color="D9E1F2")
EXTRA_IDS = ["G20604", "G20605", "G20606"]  # Pan et al. 2020 cyflumetofen cases outside the sample


def priority(row):
    if row["claude_paper_check"] in ("no", "not supported", "partly", "unclear"):
        return 1
    if row["claude_note"].startswith("paper not reachable"):
        return 1
    if row["claude_paper_check"] == "yes":
        return 2
    return 3


def main():
    s = pd.read_csv(MANUAL / "verification_sample_stage2.csv", dtype=str).fillna("")
    cases = pd.read_csv(ROOT / "data" / "clean" / "resistance_cases.csv", dtype=str).fillna("")
    extra = cases[cases["case_id"].isin(EXTRA_IDS)][
        ["case_id", "species", "ai_name", "year_of_report", "resistance_type", "location", "reference"]].copy()
    extra["claude_paper_check"] = "no"
    extra["claude_note"] = "Not in the random sample. Same paper as G20601; check for the decision on the Pan et al. 2020 cases"
    doi_file = MANUAL / "verification_dois.json"
    dois = json.loads(doi_file.read_text()) if doi_file.exists() else {}

    s["priority"] = s.apply(priority, axis=1)
    extra["priority"] = 1
    s["in_sample"] = "yes"
    extra["in_sample"] = "no"
    allc = pd.concat([s, extra], ignore_index=True)
    allc["doi_link"] = allc["reference"].map(lambda r: f"https://doi.org/{dois[r]}" if dois.get(r) else "")
    order = {"no": 0, "not supported": 1, "unclear": 2, "partly": 3, "not checked": 4, "yes": 5}
    allc["_o"] = allc["claude_paper_check"].map(order)
    allc = allc.sort_values(["priority", "_o", "reference", "case_id"]).drop(columns="_o")

    wb = Workbook()
    ws = wb.active
    ws.title = "How to check"
    lines = [
        ("Paper-level check of the stage 2 verification sample", True),
        ("", False),
        ("What this check answers: does the original paper really report the resistance case the database lists?", False),
        ("", False),
        ("For each case, open the paper and check these five points:", True),
        ("1. Species: the paper's mite is the listed species (allow old names, e.g. Tetranychus telarius = T. urticae).", False),
        ("2. Acaricide: the paper tested the listed active ingredient and reports resistance to it (not only to other compounds).", False),
        ("3. Place: the resistant population comes from the listed country.", False),
        ("4. Year: the listed year is the paper's publication year (or the collection year given in the paper).", False),
        ("5. Field population: the resistant mites were collected from the field or glasshouse, not only selected in the laboratory.", False),
        ("", False),
        ("Then fill in the yellow columns on the Cases sheet:", True),
        ("species_ok, acaricide_ok, place_ok, year_ok, field_ok: yes / no / not stated", False),
        ("resistance_ratio: the ratio or range the paper gives (e.g. 46, or 11-46); leave blank if none", False),
        ("matches_source: yes = all five points hold; partly = species and acaricide hold but another point does not; no = the paper does not report this resistance; not traceable = you cannot get the paper", False),
        ("your_note: page, table or sentence that supports your answer", False),
        ("", False),
        ("Order of work (column priority):", True),
        ("1 = check first (37 cases): my check found a problem, could not decide from the abstract, or could not open the paper. Includes three extra Pan et al. 2020 cases outside the sample.", False),
        ("2 = spot-check (21 cases): my abstract check says the paper supports the entry. Check as many as time allows.", False),
        ("3 = grey literature (45 cases): FAO reports, personal communications, short citations. Mark not traceable unless you hold the source.", False),
        ("", False),
        ("Example of a filled row:", True),
        ("G20230 | Tetranychus urticae | abamectin | species_ok yes | acaricide_ok yes | place_ok yes | year_ok yes | field_ok yes | resistance_ratio 11-46 | matches_source yes | your_note Table 2, all 12 Iranian populations", False),
        ("", False),
        ("Do not edit the grey columns; they come from the database and from my first check. Send the file back when done.", False),
        ("", False),
        ("Progress (counts of matches_source on the Cases sheet)", True),
    ]
    for i, (text, bold) in enumerate(lines, start=1):
        c = ws.cell(row=i, column=1, value=text)
        c.font = Font(name=FONT, bold=bold, size=12 if i == 1 else 10)
        c.alignment = Alignment(wrap_text=True, vertical="top")
    ws.column_dimensions["A"].width = 110
    r0 = len(lines) + 1
    n = len(allc) + 1
    for j, val in enumerate(["yes", "partly", "no", "not traceable"]):
        ws.cell(row=r0 + j, column=1, value=val).font = Font(name=FONT, size=10)
        c = ws.cell(row=r0 + j, column=2, value=f'=COUNTIF(Cases!$S$2:$S${n},A{r0 + j})')
        c.font = Font(name=FONT, size=10)
    ws.cell(row=r0 + 4, column=1, value="still to do").font = Font(name=FONT, size=10)
    ws.cell(row=r0 + 4, column=2, value=f'=COUNTBLANK(Cases!$S$2:$S${n})').font = Font(name=FONT, size=10)

    cs = wb.create_sheet("Cases")
    given = ["priority", "case_id", "in_sample", "species", "ai_name", "year_of_report", "resistance_type",
             "location", "reference", "doi_link", "claude_paper_check", "claude_note"]
    yours = ["species_ok", "acaricide_ok", "place_ok", "year_ok", "field_ok", "resistance_ratio",
             "matches_source", "your_note"]
    headers = given + yours
    for j, h in enumerate(headers, start=1):
        c = cs.cell(row=1, column=j, value=h)
        c.font = Font(name=FONT, bold=True, size=10)
        c.fill = YELLOW if h in yours else HEAD
        c.alignment = Alignment(wrap_text=True, vertical="top")
    for i, row in enumerate(allc.itertuples(index=False), start=2):
        d = row._asdict()
        for j, h in enumerate(headers, start=1):
            v = d.get(h, "")
            if h == "priority":
                v = int(v)
            c = cs.cell(row=i, column=j, value=v)
            c.font = Font(name=FONT, size=10)
            c.alignment = Alignment(wrap_text=h in ("reference", "claude_note", "your_note"), vertical="top")
            if h in yours:
                c.fill = YELLOW
            if h == "doi_link" and v:
                c.hyperlink = v
                c.font = Font(name=FONT, size=10, color="0563C1", underline="single")
    widths = {"priority": 8, "case_id": 9, "in_sample": 9, "species": 22, "ai_name": 16, "year_of_report": 8,
              "resistance_type": 14, "location": 16, "reference": 60, "doi_link": 30,
              "claude_paper_check": 12, "claude_note": 45, "resistance_ratio": 12, "matches_source": 14,
              "your_note": 40}
    for j, h in enumerate(headers, start=1):
        cs.column_dimensions[cs.cell(row=1, column=j).column_letter].width = widths.get(h, 11)
    cs.freeze_panes = "C2"
    cs.auto_filter.ref = f"A1:{cs.cell(row=1, column=len(headers)).column_letter}{n}"

    ok = DataValidation(type="list", formula1='"yes,no,not stated"', allow_blank=True)
    ms = DataValidation(type="list", formula1='"yes,partly,no,not traceable"', allow_blank=True)
    cs.add_data_validation(ok)
    cs.add_data_validation(ms)
    col = {h: cs.cell(row=1, column=j).column_letter for j, h in enumerate(headers, start=1)}
    for h in ["species_ok", "acaricide_ok", "place_ok", "year_ok", "field_ok"]:
        ok.add(f"{col[h]}2:{col[h]}{n}")
    ms.add(f"{col['matches_source']}2:{col['matches_source']}{n}")
    assert col["matches_source"] == "S", col["matches_source"]

    pp = MANUAL / "verification_papers_stage2.csv"
    if pp.exists():
        papers = pd.read_csv(pp, dtype=str).fillna("")
        ps = wb.create_sheet("Papers", 1)
        ph = ["n", "cases", "bibliography", "doi_link", "what_to_check"]
        for j, h in enumerate(ph, start=1):
            c = ps.cell(row=1, column=j, value=h)
            c.font = Font(name=FONT, bold=True, size=10)
            c.fill = HEAD
        for i, row in enumerate(papers[ph].itertuples(index=False), start=2):
            for j, v in enumerate(row, start=1):
                c = ps.cell(row=i, column=j, value=int(v) if j == 1 else v)
                c.font = Font(name=FONT, size=10)
                c.alignment = Alignment(wrap_text=True, vertical="top")
                if j == 4 and v:
                    c.hyperlink = v
                    c.font = Font(name=FONT, size=10, color="0563C1", underline="single")
        for col, w in zip("ABCDE", [5, 18, 80, 34, 50]):
            ps.column_dimensions[col].width = w
        ps.freeze_panes = "A2"

    wb.save(OUT)
    print(f"wrote {OUT.name}: {len(allc)} cases;", allc["priority"].value_counts().sort_index().to_dict())


if __name__ == "__main__":
    main()
