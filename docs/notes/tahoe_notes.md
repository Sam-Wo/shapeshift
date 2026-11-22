# Tahoe-100M transcriptomic side: notes

Script: `step2/tahoe_signatures.py`, plus the decoder `step2/tahoe_pqdecode.py`. Run it with
`py -3.13 -I -X utf8 step2/tahoe_signatures.py --stage all`, which runs pilot → fetch2 → index → build. Caches and
the index live in `step2/tahoe_cache/`. Final numbers are in `tahoe_cache/build_summary.json`,
`tahoe_cache/index_summary.json` and `step2/tahoe_download_summary.json`.

Scope, per the coordinator's three-tier revision: the 324 Tahoe names in `step2/tahoe_drugs_needed.txt`, at all three doses
(0.05, 0.5 and 5 µM). Everything is keyed on `unit` × `dose_uM`, where `unit` comes from `step2/units.csv` (the parent InChIKey).
This replaces the earlier 186-drug, compound_id-keyed design; those outputs were never written.

## 1. How the DE was computed: confirmed vs not confirmed

**The dataset card (README) does not describe how the pseudobulk DE was computed.** The table was added on 2025-07-23
as `metadata/20250328.tahoe.pseudobulk.DE.long_format.parquet` and later split into 1,026 parts. There is no methods text.
HF discussion #7 asks exactly this question and was unanswered when checked. The points below come from outside the card (web search), or from our own data checks.

Confirmed:
- **Control (from Tahoe's public plugin, not the card).** The control is the vehicle (DMSO_TF) on the same plate, for the same cell line. Tahoe's public "Tara" plugin
  (github.com/tahoebio/tahoe-public-plugins, `references/tool-routing.md` and `evals/tara/*.md`) describes this table as
  "Differential expression, drug vs same-plate vehicle". It lists pooling controls across plates or lines as wrong.
- **One contrast per drug × dose × cell line × plate.** There are 65,218 contrasts. The index has no duplicate (drug, dose, plate, line) keys.
  Where sample_metadata lists the same drug twice on one plate, the DE has a single contrast for it.
- **Our own evidence on the control.** `n_cells_ctrl` is constant within every one of the 700 (line × plate) blocks, so all drugs on a plate
  share one control pseudobulk per line. That pseudobulk is the plate's DMSO_TF cells (2–3 wells per plate, sample_metadata).
- **Columns.** The columns are DESeq2-style: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj. NaN lfc appears only for genes with zero counts.

**Not confirmed:**
- **Replicates.** We could not establish what DESeq used as replicates. Tahoe's plugin eval requires stating a "pseudoreplicate caveat", which
  implies the cells of a well and/or of the controls were split into pseudo-replicates. The scheme itself is not public.
- **Tool and shrinkage.** The tool (pyDESeq2 is likely), any LFC shrinkage, and any minimum-cell filter applied upstream are all unknown.
  The release contains contrasts with 1–20 treated cells, so no such filter was applied before release.
  - The preprint (bioRxiv 10.1101/2025.02.20.639398) describes a different analysis (Vision scores, ≥50 cells per condition).
  - Search snippets of the Cell 2026 paper say "pyDESeq2". That paper is paywalled and I did not verify it.

**Cells sharing a well are not independent.** In Tahoe, all ~50 lines are pooled in one well per drug. A contrast is a pseudobulk of one
(or two) wells vs that plate's DMSO wells. Any replicates DESeq used inside a contrast are therefore pseudo-replicates of the same well. Its p-values and SEs
ignore well-to-well variance, so treat `padj` with caution. Section 7 shows that this well effect is large.

## 2. File layout discovered
- **Files and row groups.** There are 1,026 files of 3,986,180 rows each; files 0–99 have one extra row. Each file has 3,987 row groups of 1,000 rows and a ~5.75 MB footer.
- **Combination size and gene order.** Every combination has exactly **62,710 rows**, in one **fixed gene order** (Ensembl-ID order: TSPAN6, TNMD, DPM1, ...).
  This was verified on 23 pilot combinations from 10 files: gene_name per position was identical, and row labels matched exactly at combination boundaries.
- **Row-group boundaries are not combination boundaries.** 62,710 is not a multiple of 1,000, so a combination straddles row groups and sometimes files.
  The brief's "min == max per row group" is false for about 1 row group in 63. Rows are contiguous across files, so `combo = global_row // 62710`.
  The check 1026 × 3,986,180 + 100 = 65,218 × 62,710 holds exactly. Every one of the 4,090,662 row groups was checked: stats min==max ⇔ the row group lies inside a single combination (0 mismatches, 0 label conflicts).
- **Sort order.** Files are sorted by **cell line, then plate (numeric), then drug**, not by plate. Each line's ~20.5 files interleave all
  14 plates. Pruning files with sample_metadata therefore saved only 46 of 1,026 files (those wholly inside a 0.05 µM-only plate block).
- **Plates.** Plates 1/4/7/10 are 0.05 µM, 2/5/8/11 are 0.5 µM, 3/6/9/12/14 are 5 µM, and plate 13 is mixed. Plate 14 replicates plate 6 at 5 µM (93 drugs). The other plates carry
  one 0.05 µM well (Adagrasib, not among the 324). Contrasts per dose: 0.05 µM 20,843; 0.5 µM 19,979; 5 µM 24,396. There are 379 drugs and 50 lines.

## 3. Streaming: bytes and time
Total streamed: **14.94 GB**. No DE file or obs_metadata was downloaded in full.

| Phase | What | Bytes | Time | Requests | Retries |
|---|---|---|---|---|---|
| pilot | 10 footers + 23 full combinations | 0.09 GB | 20 s | | |
| fetch v1 | 979 footers (5.69 GB) + lfc chunks of the 186-drug × 2-dose scope (2.95 GB) | 8.63 GB | 2,812 s | 1.35 M | 847 (one ~7 min throttle, absorbed by backoff) |
| fetch2 | 46 missing footers (0.26 GB) + lfc chunks of every remaining combination (5.96 GB) | 6.22 GB | 3,599 s | 2.76 M | 0 |
| index + build | local | 0 | ~2 min | | |

- **Range requests.** fetch2 decodes each log2FoldChange column chunk (one range request per row group; multi-range requests return 416) with
  `tahoe_pqdecode.py`. It uses byte offsets cached in the row-group index, so no footer is read twice.
- **Decoder check.** Validation against pyarrow had 0 mismatches over 55k values. In fetch2, row groups straddling a v1 combination are cross-checked as well (0 mismatches).
- **Disk.** Peak disk use of `step2/tahoe_cache` was 1.52 GB, consisting of per-file fragments and a row-group index.
- **Excluded from the totals.** Interactive layout probes and two single-file test fetches (~23 MB), about 0.1 GB in all.

## 4. Gene subset
- **Pilot selection.** The 23 pilot combinations spanned 10 cell lines and plates 3/5/6/8/9/11/12/14.
- **Names excluded.** Unnamed ENSG*, MT-, rRNA, Y_RNA, miRNA, snoRNA and 7SK genes were dropped.
- **Ranking.** Genes needed a finite lfc in ≥90% of pilot combinations and were then ranked by median baseMean. The **top 5,000** were kept (median baseMean ≥ 41; mostly protein-coding, plus a few lncRNAs
  such as MALAT1 and NEAT1).
- **NaN filter.** Genes NaN in >20% of kept combinations would be dropped. None were (max NaN fraction 0.6%), so all **5,000 genes** are kept.
  Remaining NaN were set to 0. Gene columns appear in Ensembl-ID order, and the list is in `tahoe_cache/candidate_genes.txt`.

## 5. Filters
- **Cell threshold.** Each contrast needs **≥100 treated cells** (`MIN_CELLS = 100`, the same as `gate_t2.py`). Across all drugs this drops 5,122 of 65,218 contrasts.
  Among the 324 target drugs it drops 4,332 of 55,297: 1,448 at 0.05, 1,293 at 0.5 and 1,591 at 5 µM.
- **Lines removed entirely:** NCI-H596 (CVCL_1571, median 13 treated cells) and NCI-H2122 (CVCL_1531, median 21). This leaves **48 lines**.
- **Lines heavily reduced:**
  - SW 1088 (CVCL_1715, median 86.5) loses 57% of its contrasts.
  - NCI-H661 (CVCL_1577, median 92) loses 57%.
  - SW 1271 (median 148) loses 24%.
- **Per-line table:** `tahoe_cache/lines_cell_counts.csv`.
- **Treated-cell distribution (target contrasts).** 1% ≤ 8, 5% ≤ 36, 10% ≤ 160, median 1,213 and 90% ≤ 2,760. 1,805 contrasts have <20 cells.
- **Absent contrasts.** Expected target contrasts (sample_metadata × 50 lines) number 55,450; 55,297 are present, so 153 are absent, probably zero cells.
- **Names missing after filtering.** All 324 names are present at 0.05 and 5 µM. Three are missing at **0.5 µM: Encorafenib, Erdafitinib, SBI-0640756.**
  They are present in the DE on plate 2, but with a median of 3–12 treated cells. Their 5 µM wells have 400–2,200 cells, so these look like failed wells rather than cytotoxicity.

## 6. Signature construction
1. **Well level.** Each DE contrast (drug, dose, line, plate) is one well-level pseudobulk lfc vs that plate's DMSO_TF.
2. **Drug × dose × line.** Contrasts are averaged across plates (wells). `n_wells` is the number of contrasts, and `n_cells_trt` is summed.
   - The plate-replicate check uses the same drug, dose and line on different plates, with independent wells and controls. Median r = **0.285** over 6,039 pairs, mostly plates 6 and 14.
3. **Unit × dose × line.** Salt and free-base drugs of a unit are averaged, after their correlation is recorded. **Salt-duplicate r: median 0.156**
   over 3,440 pairs in 69 unit × doses; after plate-line centring it is 0.04. 50 of those 69 have the duplicates on different plates, so this is mostly a cross-well check, and most of these drugs are weak.
   - **Strong units agree well.** Trametinib vs its DMSO_TF solvate: r = 0.47 / 0.65 / 0.68 at 0.05 / 0.5 / 5 µM. Doxorubicin vs epirubicin (same plate): r = 0.28 / 0.48 / 0.75.
   - **Weak units are near 0.** Fluvoxamine is 0.12–0.14 and Berbamine is −0.16 to 0.20.
   - **Doxorubicin and epirubicin share one unit (units.csv).** Epirubicin is doxorubicin's 4′-epimer; the unit key ignores stereo.
4. **Consensus.** This is the mean over lines. `magnitude` is the mean L2 norm of the per-line signatures over the 5,000 genes.

## 7. Consistency, null, and the well effect (important)
- **Definitions.** `consistency` is the mean pairwise Pearson r across lines, with a bootstrap 95% CI over lines (1,000 resamples; pairs of the same resampled line are excluded).
  The median is 0.34 at 0.05 µM, 0.37 at 0.5 µM and 0.39 at 5 µM, with a 10–90% range of about 0.25–0.51.
- **Same-plate null.** The null is the mean r between different drugs, different lines, same plate and dose, computed over every drug on the plate. It averages **0.15**
  (0.11–0.23 by plate). For different drugs on the same line, sharing a control, it is 0.31–0.35. A "fake compound" null (one random different drug per line) has a p95 of 0.12–0.25.
  98–100% of units exceed it. **Do not read this as "nearly all drugs are reproducibly active".** Pairs from one drug on one plate share a **well**, and the plate-level null does not
  capture well effects.
- **Decomposition on the plate 6/14 replicates (5 µM, 93 drugs).** Mean r by pair type, raw and after subtracting the plate × dose × line mean over all drugs:

| Pair type | Raw | Plate-line-centred |
|---|---|---|
| Same drug, different line, same plate (same well) | 0.349 | 0.306 |
| Same drug, different line, different plate (no shared well or control) | 0.246 | 0.203 |
| Same drug, same line, different plate | 0.383 | 0.270 |
| Different drug, different line, same plate | 0.149 | −0.001 |
| Different drug, different line, different plate | 0.113 | 0.000 |
| Different drug, same line, same plate | 0.349 | −0.006 |

- **What that table shows.** There is a generic treated-vs-DMSO component (r ≈ 0.11 even across plates) plus a shared-control component; both are removed by centring.
  There is also a **well effect** (about 0.10 of r) that only cross-well pairs remove.
- **Cross-well consistency.** `consistency_crosswell_centred` covers units with ≥2 plates at a dose. The median is −0.005 at 0.05 µM (n=29), 0.03 at 0.5 µM (n=30) and
  **0.16 at 5 µM** (n=93). The centred same-well medians are 0.26, 0.31 and 0.32.
  - **At 5 µM (solid evidence: 93 units with plate 6/14 replicates):** about half of the within-well consistency reproduces across wells.
  - **At 0.05 and 0.5 µM:** the cross-well estimate rests on only about 30 units each. These have a second well, mostly on plate 13 or as salt variants on another plate.
    For those units, nearly all cross-line consistency within one well is a well effect. The claim is not demonstrated for the other ~270 units, which have a single well per dose.
- **Use for gating.** Use `consistency_crosswell_centred` where available, or `consistency_centred` minus the same-plate centred same-drug baseline (≈0.27–0.32, `tahoe_null_by_plate.csv`).
  Outputs keep the raw LFC vs DMSO as specified. `step2/tahoe_plate_line_mean.parquet` holds the plate × dose × line mean over all drugs, for anyone who wants to centre.

## 8. Activity (`tahoe_activity.csv`)
- **Score.** Robust z of the signature L2 norm against **all drugs** (not only the 324) on the same plate and dose: z = (norm − median) / (1.4826·MAD). A unit is active if z > 3.
  A (plate, dose[, line]) reference needs ≥5 distinct drugs; on plates whose only 0.05 µM drug is Adagrasib it is NaN, and if MAD = 0 it is NaN.
  - **A549 and HepG2/C3A.** Each well is scored against all drugs in that line on its plate and dose. The unit z is the mean over its wells (plates × salt drugs).
  - **Consensus.** For each drug on a plate, the cross-line mean signature is scored against those of all drugs on the plate. The unit z is the mean over its (plate, drug) entries.
  - **Why score per well.** An averaged signature has a smaller noise norm than the single-well reference, so z is computed per well before averaging.
- **Active fractions:**

| Context | 0.05 µM | 0.5 µM | 5 µM |
|---|---|---|---|
| A549 | 4.3% | 5.1% | 7.3% |
| HepG2/C3A | 3.0% | 4.0% | 5.1% |
| Consensus | 2.7% | 4.7% | 5.3% |
| Any context | 6.0% | 7.4% | 8.3% |

- **Secondary measure.** The secondary measure is `frac_abs_lfc_gt1`: the fraction of genes with |lfc| > 1 in the unit's averaged signature. Medians: A549 0.3–0.6%, HepG2/C3A 1.8–2.2%,
  consensus <0.1%. HepG2/C3A is noisier (median 820 treated cells vs 1,843 for A549).
- **Caveat: focused responses.** The norm-based call favours broad responses: proteasome, HDAC, CDK, translation inhibitors. It under-calls focused ones. MEK and mTOR inhibitors
  have z ≈ 2–3, below the threshold, yet have cross-well consistency of 0.2–0.3.
- **Caveat: cell count.** The norm tracks the treated-cell count. Within plate × dose × line, the median Spearman r(norm, n_cells_trt) is −0.57 (762 groups):
  LFC sampling noise grows as cells fall. Active calls have fewer cells per well (A549 median 1,392 vs 2,062 inactive; HepG2/C3A 520 vs 880).
  "Active" therefore partly means "fewer cells", which can be real cytotoxicity or a poor well. The z formula was kept as specified.
  `n_cells_trt` is included in `tahoe_activity.csv` so this can be checked or adjusted downstream.

## 9. Counts per dose
| Dose (µM) | Units | Per-line rows | Lines per unit (median, min) | A549 rows | HepG2/C3A rows |
|---|---|---|---|---|---|
| 0.05 | 300 | 13,998 | 47, 35 | 300 | 300 |
| 0.5 | 297 | 13,850 | 47, 36 | 297 | 297 |
| 5 | 300 | 13,908 | 47, 23 | 300 | 297 |

## 10. Sanity checks
- **Names.** 324/324 names were found at 0.05 and 5 µM, and 321/324 at 0.5 µM (section 5). Each unit × dose has ~47 of the 48 lines (min 23, for cytotoxic drugs at 5 µM).
- **Magnitude vs dose.**
  - **Overall: no trend.** Medians are 23.0 / 22.1 / 23.0 (centred: 18.5 / 18.6 / 18.6), and units with 5 µM > 0.5 µM are 52%.
    Most units are indistinguishable from the bulk of drugs, and the norm is dominated by noise and well effects.
  - **Active units (n=22): clear trend.** Median magnitudes are 25.2 / 34.8 / 49.0, and 68% have 5 > 0.5 and 0.5 > 0.05 µM.
  - **Examples:** Belinostat 25 / 38 / 64; Homoharringtonine 40 / 75 / 86.
- **Strong drugs.** Their consistency and activity are high:
  - Belinostat: 0.68 at 5 µM, z = 14.
  - Dinaciclib: 0.49–0.59, z = 17–27.
  - Homoharringtonine: 0.61, cross-well 0.58, z = 26.
  - Ixazomib citrate: 0.58, z = 9.
  - Bortezomib: 0.46–0.53, z = 9.8 / 1.6 / 3.4 across the three doses.
  - Harringtonine: 0.62.
  - The median consistency for proteasome/HDAC/mTOR/MEK/CDK/topoisomerase/HSP90/tubulin MOAs is 0.36 / 0.46 / 0.43, vs 0.34 / 0.37 / 0.39 for all units.
  - **Inert-looking drugs still score high:** Riluzole, Drospirenone and Sildenafil reach raw consistency of 0.56–0.59, which illustrates section 7.

## 11. Outputs
Sizes: 0.87 GB in total.

| File | Size | Contents |
|---|---|---|
| `tahoe_signatures.parquet` (root) | 19 MB | Consensus per unit × dose_uM: unit, dose_uM, n_lines, plates, magnitude, then 5,000 float32 gene columns. 897 rows. |
| `step2/tahoe_perline.parquet` | 765 MB | unit, dose_uM, cell_line (Cellosaurus ID), cell_name, tahoe_drugs, plates, n_wells, n_cells_trt, then genes. 41,756 rows over 48 lines; zstd + byte-stream-split, float32. |
| `step2/tahoe_A549.parquet` (CVCL_0023) and `step2/tahoe_HepG2C3A.parquet` (CVCL_1098) | 19 MB each | Subsets of the per-line table. |
| `step2/tahoe_consistency.csv` | | unit, dose_uM, n_lines, plates, magnitude, consistency, its CI, centred and cross-well variants, n_plates, the same-plate null and fake p95, salt_dup_r, n_wells, tahoe_drugs, tiers, moa. |
| `step2/tahoe_activity.csv` | | unit, dose_uM, context (A549 / HepG2C3A / consensus), n_wells, z_norm, active, norm, frac_abs_lfc_gt1, plates. |
| `step2/tahoe_null_by_plate.csv` | | Per plate × dose: the null pair means (raw and centred) and the fake-consistency distribution. |
| `step2/tahoe_plate_line_mean.parquet` | 24 MB | Centring vectors: plate × dose × line mean over all drugs. |
| `step2/tahoe_download_summary.json` | | Bytes, requests, time and peak cache disk. |

## 12. Deviations from the original brief
- **Scope.** The scope is now 324 names × 3 doses, keyed on unit × dose_uM; this is the coordinator's revision. Every combination was fetched, because the activity reference needs all drugs.
- **Footers.** About 5.9 GB of footers were read (979 + 46 files). The cell-line-major file order made sample_metadata pruning nearly useless. To avoid reading them a second time, fetch2 uses a small custom
  Parquet page decoder, validated against pyarrow.
- **Cell threshold.** The threshold is 100 treated cells (coordinator), not 20.
- **Data facts.** Row groups do **not** each hold one combination. Gene order **is** fixed.
- **Consistency.** It is reported raw as specified. Section 7 shows that well effects inflate it; at low doses this rests on about 30 units each. Centred and cross-well variants are added.
- **Replicates.** The coordinator asked for a statement of what DESeq used as replicates. That could not be established from public sources (section 1).

## 13. Per-plate signatures (noise-ceiling input; built from the cache only)
`--stage perplate` reads only the cached fragments and downloads nothing. It covers units × doses observed on **≥2 different plates**, after the 100-cell threshold.
Values are raw LFC vs DMSO, **not centred**. Centring is meant to be done with `tahoe_plate_line_mean.parquet`, which joins on (plate, dose_uM, cell_line) for every row.

| File | Size | Rows | Columns |
|---|---|---|---|
| `step2/tahoe_perplate.parquet` | 281 MB | 15,263 rows: 2,925 at 0.05 µM, 3,045 at 0.5 µM, 9,293 at 5 µM | unit, dose_uM, plate (single id, string), cell_line (Cellosaurus), n_cells_trt, then the same 5,000 float32 genes |
| `step2/tahoe_perplate_consensus.parquet` | 8.7 MB | 327 rows (one per unit × dose × plate) | unit, dose_uM, plate, n_lines (40–48, median 47), n_cells_trt, then genes; values are the mean over lines |

- **Units per dose:** 29 at 0.05 µM, 30 at 0.5 µM, 93 at 5 µM.
- **Plates per unit × dose:** 130 have 2 plates, 21 have 3 and 1 has 4. All plates are kept.
- **Salt variants:** variants of a unit on the same plate and line are averaged into one row (600 rows).
- **Check against the per-line table:** for a single-drug unit, averaging these rows over plates reproduces `tahoe_perline.parquet` exactly. For multi-salt units, the per-line table averages per drug first, so the weighting differs slightly.
