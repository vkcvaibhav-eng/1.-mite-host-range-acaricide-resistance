# Raw data

Files here are never edited after they are saved. Corrections go in `data/manual/` or in the scripts, with the reason in `logs/`.

## aprd_ai_pages_acari_2026-11-13.csv

Rows for Acari species read from 48 active-ingredient pages of the Arthropod Pesticide Resistance Database (Mota-Sanchez and Wise, Michigan State University), accessed 13 November 2026. Page address: `https://www.pesticideresistance.org/display.php?page=ai&aiId=<ai_id>`.

How it was collected: the database's search form cannot be reached from the analysis workspace, so each page was read with a web reader and every row with order "acari" was transcribed: species, family and host exactly as the database shows them (lower case and spelling kept, including "hydrangaea", "galenod" and the family "pentatomidae" given for *Halotydeus destructor*).

Limitations:
- Covers only the 48 compounds listed; the database holds more. Species resistant only to compounds not read here are missing.
- Gives presence of a species on a compound page, not the number of cases.
- Rows were transcribed through a page reader, not exported; the full export requested from Dr David Mota-Sanchez replaces this file in stage 2.

The 48 compounds were chosen as the dedicated acaricides plus the organophosphates, carbamates, pyrethroids and organochlorines listed on the *Tetranychus urticae* page (arId 536).

The genus menu of the database search page also lists *Eotetranychus* and *Aculops*; no species of either genus appeared on the 48 pages read.

Permission note from the database: permission must be granted if a significant part of a paper is based on its data.

## aprd_species_page_ids_2026-11-13.csv

Species-page IDs found so far, for direct links in stage 2.

## literature_resistant_species_2026-11-13.csv

Phytophagous mites with published resistance evidence that were not on the 48 database pages. All four are marked "to verify in stage 2" against the case definition in `protocol.md`.

## Added later on 13 November 2026 (collected through the user's Chrome browser)

### aprd_acari_ag_species_2026-11-13.csv

The complete species table from the database search with Order = acari (adding "acarina" gave only one more species, a phytoseiid). All rows with group AG are kept: arId, species and family as listed, number of cases, and, for the phytophagous species, the number of citations and the countries on each species page. This replaces `aprd_ai_pages_acari_2026-11-13.csv` as the source of the species list; the older file is kept unchanged as a record.

### aprd_species_ai_cases_2026-11-13.csv

For each phytophagous species page: every active ingredient listed and its case count. The per-ingredient counts add up to the species total for every species except *Tetranychus urticae* (570 against 572 on the search table).

### bold_summary_2026-11-13.csv

BOLD Data Portal summary counts per species name (`tax:species:<name>`): specimens, records with any COI region, and records with the COI-5P barcode region.
