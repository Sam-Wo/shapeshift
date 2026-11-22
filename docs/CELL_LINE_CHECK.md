# ShapeShift cell-line check: Cell Painting in Tahoe-100M cell lines

Scripts: `cellcheck/obs_cell_lines.py` (finds which lines were screened) and `cellcheck/cellline_match.py` (compound matching).
Outputs: `cellcheck/cellline_matches.csv` (one row per Tahoe drug × dataset), `cellcheck/data/summary.json`, and `cellcheck/data/dataset_compounds_*.csv`.
Only small metadata files were downloaded, about 15 MB in `cellcheck/data/`.

Compounds were matched with overlap.py's own functions (`std_mol`, `std_key`, `graph_hash`, `norm_name`, `METALS`), so the handling is the same as in Step 1:
- **Tier 1:** the standardized, stereo-flat InChIKey is equal. Both Tahoe keys (`tahoe_key`, `pubchem_key`) are compared with both the dataset's own InChIKey and the `std_key` of the dataset's SMILES. The split is 1a (achiral) or 1b (stereo unresolved).
- **Tier 2:** InChIKey block 1 is equal.
- **Tier 2b:** same molecular formula and same heavy-atom graph. Here the formula comes from the standardized molecule on both sides; overlap.py took the JUMP formula from JUMP's own InChI.
- **Tier 3:** a metal was removed, or the match is a name bridge (`norm_name`). Tier 3 is used only when there is no structure match.

## 1. Are A549 and HepG2 in Tahoe-100M?

`cell_line_metadata` has 1000 rows and 102 lines, so it lists more lines than were screened. The 50 screened lines were taken from `metadata/obs_metadata.parquet` without downloading it (2.3 GB). Only the footer and the dictionary-encoded `cell_line`/`cell_name` columns of row groups 0, 50 and 95 were read, over HTTP range requests.
- The pandas metadata says `num_categories: 50` for both columns.
- Each of the three row groups contains all 50 (Cellosaurus ID, name) pairs.
- All 50 IDs are present in `cell_line_metadata`.
- The list is saved in `cellcheck/data/obs_cell_lines.json`.

| Query | In the screened 50? |
|---|---|
| **A549** (CVCL_0023, ACH-000681, lung) | **Yes**, an exact match. |
| **HepG2 / Hep G2 / HEPG2** (Hep-G2, CVCL_0027) | **No.** The parent line is not in Tahoe under any spelling, and it is also absent from the 102-line table. |
| **HepG2/C3A** (Hep-G2/C3A, **CVCL_1098**, ACH-001021, liver) | **Yes, screened.** It is a derivative of Hep-G2, not Hep-G2 itself. |

**HepG2/C3A is a derivative, not HepG2.** Cellosaurus gives Hep-G2 (CVCL_0027) as its parent (`HI` line) and lists it as a patented subline (US5290684, ATCC CRL-10741) with its own accession. It comes from the same donor, but it is a separate line and its expression will differ from the parent.

**U2OS is not among the 50 lines and is not in the 102-line table either.**

The 50 screened lines:

A-172, A-427, A498, A549, AN3 CA, AsPC-1, BT-474, C-33 A, C32, CFPAC-1, CHP-212, COLO 205, H4, HCT15, HEC-1-A, HepG2/C3A, HOP62, Hs 766T, HS-578T, HT-29, hTERT-HPNE, J82, KATO III, LoVo, LOX-IMVI, LS 180, MIA PaCa-2, NCI-H1573, NCI-H1792, NCI-H2030, NCI-H2122, NCI-H23, NCI-H2347, NCI-H460, NCI-H596, NCI-H661, Panc 03.27, PANC-1, RKO, RPMI-7951, SHP-77, SK-MEL-2, SNU-1, SNU-423, SW 1088, SW 1271, SW 900, SW1417, SW48, SW480.

## 2. Cell Painting datasets in these lines

### Sources

**LINCS, cpg0004 (2016_04_01_a549_48hr_batch1)**
- Platemaps, 28 maps for 140 plates: [broadinstitute/lincs-cell-painting](https://github.com/broadinstitute/lincs-cell-painting) at `metadata/platemaps/`. These were spot-checked against gallery `cpg0004-lincs/broad/workspace/metadata/platemaps/`: the file sizes match, and two platemaps plus `barcode_platemap.csv` are byte-identical.
- Structures come from the Repurposing Hub sample files of four releases (2017-03-27, 2018-05-16, 2018-09-07 and 2020-03-24; the repo's `metadata/moa/clue/`). Lookup order:
  1. exact `broad_id`, newest release first
  2. `deprecated_broad_id`
  3. the 13-character BRD prefix, used only when its structures give a single `std_key`
- The gallery's `external_metadata.tsv` supplies InChIKey14 and names for the rest.
- 1570 broad_samples were found. 1551 resolved to a structure, giving 1544 unique structures, and 19 have no structure:
  - 18 have neither structure nor name in any source.
  - 1 (curcumin) has an ambiguous prefix.
- The `mmoles_per_liter` values are in µM. Way et al. 2022, Cell Systems (PMC10246468), gives "0.04 μM to 10 μM" across six doses, and the repo's consensus code uses `primary_dose_mapping = [0.04, 0.12, 0.37, 1.11, 3.33, 10, 20]`. The 20 µM wells are only the bortezomib and MG-132 controls.
- Exposure is 48 h. Way et al. state "exposed them to compound treatment at various doses for 48 hours", and note that L1000 in A549 was done at 24 h.

**LINCS Batch2 pilot (2017_12_05_Batch2)**
- This is profiles-only in the gallery: `profiles/2017_12_05_Batch2` exists, but there is no images folder under that name.
- It covers A549, MCF7 and U2OS at 6, 24 and 48 h; Way et al. confirm "6, 24, and 48 hours". Only the A549 platemaps were used.
- 90 of its 427 samples have no structure, so it is reported separately.

**CPJUMP1, cpg0000-jump-pilot/source_4**
- Compound metadata and platemap: `JUMP-Target-1_compound_metadata.tsv` and `JUMP-Target-1_compound_platemap.tsv` from [jump-cellpainting/JUMP-Target](https://github.com/jump-cellpainting/JUMP-Target).
- Conditions: `benchmark/output/experiment-metadata.tsv` from [2024_Chandrasekaran_NatureMethods](https://github.com/jump-cellpainting/2024_Chandrasekaran_NatureMethods).
- **Dose is 5 µM, verified.** The paper's Methods (PMC11166567) say "The treatment compounds were assayed at 5 μM", and the JUMP-Target README recommends 5 µM.
- Time points are 24 h and 48 h in both A549 and U2OS. Main batch 2020_11_04_CPJUMP1:
  - A549: 4 plates at 24 h; 4 parental + 4 Cas9 plates at 48 h; extra density variants.
  - U2OS: 4 plates at 24 h; 4 plates at 48 h.

**EU-OPENSCREEN bioactives, cpg0036**
- Per-site platemaps and `external_metadata.tsv`: gallery `cpg0036-EU-OS-bioactives/{FMP,IMTM,MEDINA,USC}/workspace/metadata/`. The metadata file has 2464 ECBL compounds with SMILES and InChIKey, and is identical at all four sites (same md5).
- **Doses, verified in the iScience 2025 Methods (PMC12084007):**
  - 10 µM for every compound, at a single dose.
  - 5 µM for the nocodazole and tetrandrine controls. The platemaps also list 28 control wells per control at concentration 0, which looks like a metadata quirk.
  - 24 h exposure, 4 replicates.
- **HepG2 (the parent, CVCL_0027) was imaged at all four sites, FMP, IMTM, MEDINA and USC, with 28 plates each. U2OS was imaged at FMP only, with 28 plates.**
- The gallery holds images and metadata only. Aggregated profiles are on Zenodo, concept DOI 10.5281/zenodo.13309565. The file list checked was version 1.0.0 (record 13309566, `Aggregated_Profiles.zip`, 1.3 GB). A v3.0.0 also exists but its file list was not checked. Nothing was downloaded.
- The metadata has no compound names, only Probes & Drugs `pdid`, so EU-OS could be matched on structure only.

### Matched compounds

"Unique compounds" collapses stereoisomers and salts by the dataset's `std_key`. Step-2 overlap is counted on Tahoe drug name, using `step2/compounds.csv` `tahoe_drugs`. Tier 2 (block 1) gave no matches in any dataset. MoA counts use Tahoe `moa-fine`, excluding "unclear", and count unique tier 1/2 compounds.

| Dataset | Cell line (vs Tahoe) | Dataset compounds | Tahoe drugs matched (1a / 1b / 2b / 3) | Tier 1+2 drugs | Unique compounds (tier 1+2) | Also Step-1 JUMP tier 1/2 | In Step-2 set (drugs / JUMP compounds) | MoA classes with ≥3 compounds |
|---|---|---|---|---|---|---|---|---|
| **LINCS cpg0004** | A549, exact | 1570 samples / 1544 structures | 92 / 108 / 3 / 3 | **203** | **185** | 200 | **99 / 92** | 13 |
| LINCS Batch2 pilot (A549 plates) | A549, exact | 427 samples / 315 structures | 27 / 38 / 1 / 3 | 66 | 60 | 63 | 30 / 27 | 6 |
| **CPJUMP1 cpg0000** | A549, exact (also U2OS) | 306 (300 structures) | 12 / 9 / 2 / 0 | 23 | 21 | 23 | 10 / 9 | 0 |
| **EU-OS cpg0036** | HepG2 parent; Tahoe has the C3A derivative (also U2OS) | 2464 + 2 controls (2463 structures) | 89 / 97 / 0 / 0 | **186** | **172** | 183 | **94 / 90** | 14 |

**Cross-dataset overlap:**
- A549 total, LINCS ∪ CPJUMP1: 211 Tahoe drugs. CPJUMP1 adds 8 drugs that LINCS lacks.
- Adding the Batch2 pilot brings the A549 total to 221 drugs, 110 of them in Step 2.
- LINCS ∩ EU-OS: 121 drugs.

**Notes on the matches:**
- All 23 CPJUMP1 matches are also Step-1 JUMP matches, as expected, because the JUMP-Target compounds are in cpg0016. Four of the 23 are CPJUMP1 positive controls: fedratinib, filgotinib, 5-azacytidine and lonafarnib.
- A few drugs match here but were unmatched in Step 1: bortezomib (LINCS, EU-OS), pentoxifylline (LINCS), pasireotide (LINCS), digitoxin (EU-OS) and trimetrexate (EU-OS).
- Tier 3 matches are all LINCS name bridges, and each pairs a salt or ester form with the bare compound:
  - Ixazomib citrate → ixazomib. This is the prodrug paired with its active form, so the match is plausible.
  - Methylprednisolone succinate → methylprednisolone. This is an ester prodrug.
  - **Estrone sulfate → estrone.** This is a different molecule; treat it as no match. It is still counted in the tier-3 column above and still appears in the CSV, so LINCS and Batch2 really have 2 usable tier-3 hits each.
- **Tier 1b caveat:** some stereo-flat matches are a different stereoisomer, or a racemate, of the Tahoe drug. Examples:
  - dexmedetomidine ↔ medetomidine
  - levobupivacaine ↔ bupivacaine
  - (R)-verapamil ↔ verapamil
  - citalopram ↔ escitalopram
  - quinidine ↔ quinidine and quinine samples

  The `name_agrees` column flags these; most other name differences are synonyms, such as KD025 = belumosudil. (S)-crizotinib ↔ crizotinib (LINCS, Batch2) is probably an enantiomer pair: the Hub sample is likely (R)-crizotinib. That was not checked.

**MoA classes with ≥3 unique compounds** (counts in parentheses):
- **LINCS:**
  - DNA synthesis/repair (15)
  - COX (8)
  - EGFR/ERBB (6)
  - adrenoceptor agonist (5)
  - JAK/STAT (4)
  - mTOR (4)
  - multi-TK (4)
  - DNMT (3)
  - AR antagonist (3)
  - microtubule (3)
  - RAF (3)
  - CDK (3)
  - MEK (3)
- **EU-OS:**
  - DNA synthesis/repair (16)
  - adrenoceptor agonist (6)
  - mTOR (6)
  - other TK (5)
  - EGFR/ERBB (5)
  - other MAPK (4)
  - microtubule (4)
  - MEK (4)
  - multi-TK (4)
  - CDK (3)
  - COX (3)
  - RAF (3)
  - PI3K/AKT (3)
  - Sonic (3)
- **Batch2:**
  - DNA synthesis/repair (10)
  - multi-TK (5)
  - EGFR (4)
  - CDK (3)
  - JAK/STAT (3)
  - mTOR (3)
- **CPJUMP1:** none.

### Doses vs Tahoe (0.05 / 0.5 / 5 µM, 24 h)

Tahoe's 24 h exposure is from the preprint Methods: "Cells were exposed to drug media for 24 hours". All 379 Tahoe drugs use 0.05, 0.5 and 5 µM.

| Dataset | Doses tested (µM) | Exposure | Nearest dose to each Tahoe dose (fold) | Time match |
|---|---|---|---|---|
| LINCS | 0.04, 0.12, 0.37, 1.11, 3.33, 10 (1550 of 1570 samples; 20 use a shifted or lower series, e.g. 0.004–1). Controls at 20. | 48 h | 0.05→0.041 (0.82×), 0.5→0.37 (0.74×), 5→3.33 (0.67×, i.e. 1.5× below) or 10 (2× above). | No (48 vs 24 h) |
| LINCS Batch2 | mostly 0.1, 1 and 10; some 6-point series | 6, 24, 48 h | 2× above each Tahoe dose | Yes (24 h plates exist) |
| CPJUMP1 | 5 | 24 h and 48 h | 5→5 (1×), top Tahoe dose only | Yes (24 h) |
| EU-OS | 10 (controls 5) | 24 h | 5→10 (2×); nothing near 0.05 or 0.5 | Yes |

The per-drug doses are in `cellline_matches.csv`, in the columns `dataset_doses_uM`, `tahoe_dose_nearest`, `tahoe_doses_within_1p5x` and `tahoe_doses_within_2x`. The 1.5× window has a 2% tolerance so that the nominal 3.333 µM (3.33 in the files) counts as within 1.5×. With that tolerance, all 203 LINCS matches have all three Tahoe doses within 1.5×.

## 3. Recommendation

**LINCS cpg0004 (A549) best supports a cell-line-matched morphology↔transcriptome comparison.**

Why:
- A549 is the same line in both datasets (CVCL_0023).
- It has the most matches: 203 Tahoe drugs and 185 unique compounds, 99 of them in the Step-2 set.
- Its 6-point series puts a dose within 1.5× of each Tahoe dose (0.04, 0.37 and 3.33 µM; 5 µM is bracketed by 3.33 and 10). The comparison can therefore be dose-matched, which no other dataset allows.
- Profiles are already in the gallery.

Trade-offs:
- Exposure is 48 h against Tahoe's 24 h.
- The data are older: 2016, Cell Painting protocol v2, a single site.
- Hub samples are stereo-specific while matching is stereo-flat, so the tier-1b racemate and enantiomer pairs need review.

**CPJUMP1** is the cleanest exact-condition subset: A549, 5 µM and 24 h all equal Tahoe's top dose. However, it has only 23 drugs and no MoA class with ≥3 compounds. Its best use is as a validation subset and a bridge, because the same compounds were imaged in A549 and U2OS. That lets you measure directly how much the morphology changes between A549 and U2OS for compounds also in JUMP.

**EU-OS** has good coverage (186 drugs, 14 MoA classes) at Tahoe's 24 h. Its limits:
- It is not an exact cell-line match: it imaged the HepG2 parent, while Tahoe screened the C3A derivative.
- It has a single 10 µM dose, 2× Tahoe's top dose.
- Its profiles have to be fetched from Zenodo.

Its four-site HepG2 design is useful for estimating site and batch noise.

**Compared with JUMP U2OS (Step 1, 357 tier 1/2 drugs, 10 µM and 0.625 µM):** JUMP covers about 1.75× more Tahoe drugs, but U2OS is not one of Tahoe's 50 lines, so every JUMP comparison is across cell lines. JUMP also carries multi-source batch structure.

The practical plan:
1. Use LINCS A549 as the primary matched-line comparison.
2. Add CPJUMP1 (A549 24 h) as the exact-condition check.
3. Use the roughly 121 drugs shared by LINCS and EU-OS, and the CPJUMP1 A549 vs U2OS pairs, to quantify how much cell-line mismatch costs before relying on JUMP U2OS for its extra drugs.

## Caveats and items not verified

- **LINCS unresolved samples:** 19 samples (18 with no structure or name, plus curcumin) could hide a few extra matches.
- **LINCS Batch2 images:** images were not located under the Batch2 name (`images/2018_02_23_LKCP_DBG` may be related). Only profiles are confirmed.
- **EU-OS controls:** the 28 control wells with concentration 0 are unexplained. Nocodazole and tetrandrine are not Tahoe drugs, so they do not affect the counts.
- **Unparseable SMILES:** two CPJUMP1 SMILES (briciclib, fludarabine phosphate) have broken phosphate valence and could not be parsed. Neither is a Tahoe drug.
- **Tier-2b formula:** it is computed from the standardized molecule on both sides. overlap.py used JUMP's own InChI formula on the JUMP side.
- **No PubChem bridge:** overlap.py's tier 3b PubChem parent and name bridge was not re-run here. EU-OS has no names, so it has no name bridge.
