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

## 2026-11-13 · Stage 1 checkpoint

15. **Go.** The user approved moving to stage 2 with 34 resistant species (24 with COI), using the plan's fallback: susceptible species kept as zeros and research effort modelled.
16. **Protocol approved** and fixed as of this date.
17. **Database use.** The user considers citation and acknowledgement of the Arthropod Pesticide Resistance Database sufficient because the data are openly viewable. Open point: the database states that permission must be granted when a significant portion of a paper is based on its data. To be settled before submission.
18. **Permission requested.** On 13 November 2026 the user emailed Dr David Mota-Sanchez (Cc Dr John C. Wise) asking consent to use the database's Acari records and offering to share the manuscript before submission. Awaiting reply.

## 2026-11-13 · Stage 2

19. **Case records** for all 291 species × compound pairs were taken from the database through the user's browser (1,027 cases; transfer checked against a row count and case-ID checksum).
20. **IRAC groups** assigned by hand in `data/manual/irac_map.csv` from the current IRAC classification. Obsolete compounds that IRAC never classified (bridged diphenyls such as chlorobenzilate, chlorfenson, tetrasul; selenium sulfide; zineb; dienochlor; MNFA; genite; azobenzene) are coded NC; UN (dicofol, bromopropylate, benzoximate, chinomethionate) and NC are not counted as mode-of-action groups, but are kept in a separate count of chemical classes. The hexythiazox + organotin mixture is coded MIX. Bifenazate is placed in 20D although the database lists it as unknown.
21. **Lab-selected only** cases (22) are excluded from the main analysis; the population did not evolve resistance in the field.
22. **Duplicate entries.** Xue et al. (2020) is entered twice in the database (reference 162 at country level, 163 at site level). The country-level rows (9) are excluded. *T. cinnabarinus* cyflumetofen (references 117 and 118) is the same paper, but the extra row is already excluded as lab-selected.
23. **Case definition (decided 13 November 2026, option A).** Database cases are accepted as recorded for the main analysis; sensitivity runs use field-evolved cases only and journal-supported cases only. Logged in the protocol change log.
24. **Open: first registration year** of each active ingredient (needed for H3) is not yet compiled.
25. **Open: literature cases** for *Oligonychus coffeae*, *Tetranychus evansi* and *Tetranychus truncatus*, and the documented literature search, are not yet done.
