# Instructions for Claude in this repository

Project: comparative test of whether host range and phylogeny predict acaricide resistance in phytophagous mites. Full plan: `docs/project_plan.md`.

- Work one stage per session. Stop at each stage checkpoint and wait for the user's approval before starting the next stage.
- Never edit files in `data/raw/`. Write curated output to `data/clean/`.
- Log every excluded record, name change and judgement call in `logs/` with the reason.
- Every number used in the manuscript must be written by a committed script into `results/tables/`; never type numbers by hand.
- Ask before any step that costs money (creating or starting a Google Cloud machine) and report the running cost.
- Ask before installing anything large.
- Fix random seeds. Record program versions.
- `protocol.md` is fixed before modelling; any later change is logged in its change log with date and reason.
- Write full journal names, never abbreviations.
