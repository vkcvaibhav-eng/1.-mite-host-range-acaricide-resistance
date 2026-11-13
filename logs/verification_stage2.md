# Verification of the stage 2 case table

## Check 1: database transcription (13 November 2026, done by Claude through the user's browser)

All 99 database cases in the random 10% sample (`data/manual/verification_sample_stage2.csv`, seed 2026) were looked up again on their case pages (`display.php?page=rrc&aiId=…&arId=…`). For each, species, active ingredient, year of report, country and reference were compared with `data/clean/resistance_cases.csv`.

Result: 99 of 99 match; 0 errors. The 100th sampled case (LIT008) is a literature case, checked against its paper's abstract when it was added.

This confirms the transfer from the database. It does not test whether the database itself reports each paper correctly.

## Check 2: database against the original papers

### 2a. First pass by Claude (13 November 2026)

Columns `claude_paper_check` and `claude_note` in the sample file. Sources: PubMed and Crossref abstracts, and full pages on ScienceDirect read through the user's browser after he signed in. Oxford University Press, BioOne and Springer pages were behind bot checks and were not opened.

| Result | Cases |
|---|---|
| yes, paper supports the database entry | 21 |
| partly (species and compound match; country or compound not named in the abstract) | 5 |
| unclear (abstract too general to tell) | 7 |
| not supported (paper is on other compounds) | 1 (G20612) |
| no (paper says the opposite) | 1 (G20601) |
| paper not reachable online | 20 |
| grey literature or short citation, not traceable | 45 |

Of the 35 cases checked, 2 do not match (G20601, G20612).

G20601: Pan et al. 2020, *Journal of Economic Entomology*, found no cyflumetofen resistance in the Liangping, Wanzhou, Daying and Anyue populations. The database lists cyflumetofen cases from this paper (G20601 in the sample; G20604, G20605 and G20606 outside it). Decision on these four cases pending with the user.

G20612: Adesanya et al. 2018 deals with mite growth inhibitors; fenpyroximate is not in the abstract. Needs the full text before any decision.

Abstract-level checks cannot confirm a country or a single population; "partly" and "unclear" are left for the full-text check.

### 2b. Check by V. K. Chaudhari (to do)

Column `matches_source`: yes, no, or not traceable. The error rate reported in the paper is the share of traceable cases that do not match.
