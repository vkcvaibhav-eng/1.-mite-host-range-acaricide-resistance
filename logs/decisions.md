# Decisions log

Each entry: date, stage, decision, reason. Name decisions are also listed in `data/manual/name_map_stage1.csv`.

## 2026-11-13 · Stage 1

1. **Source of the resistant-species list.** The search form of the Arthropod Pesticide Resistance Database cannot be reached from the analysis workspace (the network proxy refuses the site). Instead, 48 active-ingredient pages were read through a web reader and the Acari rows were transcribed into `data/raw/aprd_ai_pages_acari_2026-11-13.csv`. The 48 compounds are those on the *Tetranychus urticae* page that are dedicated acaricides or broad-spectrum older chemistry. This is a minimum list; the full export replaces it in stage 2.
2. **Genera seen but not reached.** The database's genus menu lists *Eotetranychus* and *Aculops*, but neither appeared on the 48 pages. They are not added by guesswork; stage 2 resolves them from the full export.
3. **Literature additions.** Four species with published resistance evidence but absent from the 48 pages were added and marked "to verify": *Oligonychus coffeae*, *Aculops lycopersici*, *Tetranychus truncatus*, *Tetranychus evansi*.
4. **Scope.** Rows from Phytoseiidae, Ixodidae, Acaridae, Dermanyssidae, Macronyssidae, Varroidae and Stigmaeidae were excluded (not phytophagous or not in the three study families). *Halotydeus destructor* (Penthaleidae; the database gives "pentatomidae") is phytophagous but outside the three families and was excluded. All exclusions are in `logs/feasibility_excluded.csv`.
5. **Synonyms merged:** *T. bimaculatus* and *T. althaeae* into *T. urticae*; *T. atlanticus* into *T. turkestani*. Current combinations used: *Amphitetranychus viennensis*, *Aculops pelekassi*. Spelling corrected: *T. hydrangeae*.
6. ***T. cinnabarinus*** kept separate in the main list; the summary also reports counts with it merged into *T. urticae*, as the protocol requires both treatments.
7. **Names flagged for checking in stage 2:** *Brevipalpus phoenicis* (the citrus record may be *B. yothersi* after the 2015 split), *Tetranychus crataegi*, *Aculus malivagrans*.
8. **GenBank query design.** Counts use the NCBI taxonomy ID with descendants where the name exists in NCBI, otherwise the name as organism. Marker filters use gene symbols and title phrases with a 150–3,000 bp length window, matching the window used for the family counts in the project plan. Every query is logged in `logs/feasibility_queries.tsv`.
9. **BOLD not counted.** The BOLD portal is also refused by the workspace proxy. BOLD counts are deferred; they can be added if boldsystems.org is allowed for this workspace or the script is run on a laptop.

## 2026-11-13 · Stage 1, rerun with full data

10. **Full species list obtained.** Through the user's Chrome browser the database search returned the complete Acari table. Stage 1 now uses it instead of the 48-page transcription. New phytophagous species found: *Eotetranychus hicoriae* (1 case). *Aculops lycopersici* is in the database (1 case, methamidophos, Egypt), so it moved from the literature list to the database list.
11. ***Aculops pelekassi* and *Aculus pelekassi*** have separate database pages (3 cases in total); merged.
12. **BOLD counts added** through the browser on the same day. A species counts as having COI if either GenBank or BOLD has a record. BOLD adds COI for *Aculus cornutus* (1 record) and raises *T. pacificus* (4 in GenBank, 55 in BOLD) and *T. ludeni* (21 and 140).
13. ***T. urticae* case total**: 575 = 572 (*T. urticae*) + 2 (*T. bimaculatus*) + 1 (*T. althaeae*).
14. **NCBI throttling.** The second run met HTTP 429 responses; the pause between requests was raised to 1 s with longer back-off. Counts are unchanged where the two runs overlap.
