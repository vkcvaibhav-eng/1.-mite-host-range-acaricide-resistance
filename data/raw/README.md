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

## Stage 2 files (13 November 2026, collected through the user's Chrome browser)

### aprd_cases_part1_2026-11-13.txt and aprd_cases_part2_2026-11-13.txt

Every resistance case on the 291 species × active-ingredient pages (`display.php?page=rrc`) for the phytophagous species in the three families: arId, ai_id, case ID, year of report, type of resistance, location, reference index. Part 1 keeps the browser's space-separated layout and full type labels; part 2 is pipe-separated with type codes (F field-evolved, U unspecified, L lab-selected only, LF lab-selected after field-evolved resistance, L+F, L+LF). Checksum from the browser: 1,027 rows, sum of numeric case IDs 7,748,965; `scripts/02_resistance.py` stops if either differs. Garbled characters in a few Turkish and Brazilian place names were written in plain ASCII (Cicekli, Sao Paulo) during transfer.

### aprd_references_2026-11-13.txt

The 283 references cited by those cases (index|citation as shown in the database). The same plain-ASCII fix was applied to three names (Doker, Cagatay, page range 996-1004).

### ppdb_introduction_dates_2026-11-13.txt

"Introduction & key dates" and "Availability status" for 102 active ingredients, read on 13 November 2026 through the user's browser from the University of Hertfordshire Pesticide Properties DataBase (PPDB); abamectin and milbemectin from its BioPesticides DataBase (BPDB) and famphur from its Veterinary Substances DataBase (VSDB). Columns: ai_id | name used for matching | record (database:report number) | introduction text | availability. Not found in any of the three: chlorfensulfide, genite, MNFA, selenium sulfide, sulphenone.
