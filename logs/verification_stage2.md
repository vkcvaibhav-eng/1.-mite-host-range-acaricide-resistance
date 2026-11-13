# Verification of the stage 2 case table

## Check 1: database transcription (13 November 2026, done by Claude through the user's browser)

All 99 database cases in the random 10% sample (`data/manual/verification_sample_stage2.csv`, seed 2026) were looked up again on their case pages (`display.php?page=rrc&aiId=…&arId=…`). For each, species, active ingredient, year of report, country and reference were compared with `data/clean/resistance_cases.csv`.

Result: 99 of 99 match; 0 errors. The 100th sampled case (LIT008) is a literature case, checked against its paper's abstract when it was added.

This confirms the transfer from the database. It does not test whether the database itself reports each paper correctly.

## Check 2: database against the original papers (to do, by V. K. Chaudhari)

Column `matches_source` in the sample file. 55 of the 100 cases cite a journal paper; the other 45 cite FAO reports, personal communications or short citations that can no longer be traced, and should be marked "not traceable" rather than "no". The error rate reported in the paper is the share of traceable cases that do not match.
