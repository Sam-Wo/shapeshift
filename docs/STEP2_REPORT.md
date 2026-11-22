# ShapeShift – Step 2 report: do morphology and the transcriptome agree on compound similarity?

*2026-11-22. Three-tier design. Main script: `step2/compare_tiers.py`. Results: `step2/results/`; figures: `step2/figures/`. Step 3 was not started.*

## Verdict

**Morphology and the transcriptome agree only weakly on compound similarity.** On full compound sets the agreement is close to zero. On the compounds where Tahoe's transcriptome is measured twice, the agreement is real and consistent across tiers, but still modest. The main thing limiting it is noise in Tahoe's single-well signatures. Cell-line mismatch is not the limit.

1. **No pre-registered primary test passes.** The primary test is a within-MoA Mantel test on annotated compounds at the top matched dose. It gives p = 0.25 (T1), 0.31 (T2) and 0.035 (T3). None survives Bonferroni correction across the three primaries (α = 0.0167). In T1 and T2 the observed r roughly equals the within-MoA null mean, so whatever agreement exists is what shared MoA class alone predicts.
2. **Removing Tahoe's plate effect reveals a small agreement in every tier.** Every drug on a Tahoe plate has the same DMSO control subtracted, which creates a shared-control (plate) component. With it removed, the all-compound Mantel r is 0.020 (T1, p = 0.04), 0.048 (T2, p ≤ 0.0005) and 0.037 (T3, p = 0.003). T1 and T2 have opposite confounds, yet they point the same way.
3. **Agreement depends on how well Tahoe measured the compound.** In each tier, roughly 45–61 compounds were profiled on two Tahoe plates, and the rest on one.

   | Tier | Two-plate units: r [CI] | One-plate units: r |
   |---|---|---|
   | T1 | 0.066 [−0.01, 0.16] | 0.009 |
   | T2 | 0.098 [−0.02, 0.22] | 0.027 |
   | T3 | 0.136 [0.02, 0.26] | 0.019 |

   - **Not stronger drugs.** The two-plate units do *not* have stronger signatures (median centred magnitude 14–15 vs 20–22), so the gap is not a potency effect. The likeliest reason is that a two-well average carries less noise than a single well.
   - **Ceiling on the same units.** The noise ceiling for the two-plate units is about 0.39–0.46. Their observed r is therefore about 15–35% of what the data could support, and less than that once the benefit of averaging two plates is counted.
   - **Single-plate ceiling unknown.** It can't be measured, because there is no second well. The single-plate signatures also look noisier (larger magnitude).
4. **Each modality carries MoA information on its own.** Morphology retrieves MoA at mAP 0.26–0.37 and the transcriptome at 0.17–0.19, against a null of 0.11. The transcriptome result survives plate-centring. Neither space is uninformative; they just don't rank compounds the same way beyond class membership.
5. **Matching the cell line does not help.**
   - T1 (A549 against A549) agrees no better than T3 (U2OS against the Tahoe consensus), on the full sets or on the same 92 compounds.
   - With plate-centred signatures, A549 ranks 10th of 48 Tahoe lines for T1 and HepG2/C3A 27th for T2. Both are within about 0.005 of the median line, so there is no detectable matched-line advantage.
   - CPJUMP1 shows that changing the cell line keeps 39–52% of the morphology similarity structure. That is a real effect, but it is not what limits agreement here.

**For Step 3: proceed only in a narrowed form.** Predicting the whole transcriptome from morphology (or the reverse) across compound libraries has little shared structure to learn from data of this quality: r ≈ 0.02–0.05 on full sets. The signal is concentrated in compounds with replicated, high-dose transcriptomic measurements, where r ≈ 0.07–0.14. A Step 3 should therefore restrict itself to such compounds, or to a transcriptome source with replicate wells, and treat a positive outcome as uncertain.

## What was run

| Tier | Morphology | Cell line vs Tahoe | Exposure | Units analysed (annotated, classes ≥3) | Primary pairing |
|---|---|---|---|---|---|
| **T1** | LINCS cpg0004: 136 plates, 6 doses | A549 = A549 (exact) | 48 h vs 24 h | 187 (85, 14) | LINCS 3.33 µM vs Tahoe 5 µM, A549 signature |
| **T2** | EU-OPENSCREEN cpg0036: 4 sites × 28 plates | HepG2 vs HepG2/C3A subline | 24 h vs 24 h | 171 of 172 (82, 14)¹ | 10 µM vs Tahoe 5 µM, HepG2/C3A signature |
| **T3** | JUMP cpg0016: 416 plates | U2OS, not in Tahoe | 48 h vs 24 h | 175 (82, 15) | 10 µM vs Tahoe 5 µM, cross-line consensus |
| Calibration | CPJUMP1 cpg0000: 16 plates | A549 vs U2OS | 24 h and 48 h | 306 compounds | – |

¹ One T2 unit has no HepG2/C3A signature at 5 µM.

**Compound sets**
- A unit is Tahoe's standardized parent structure, so salt and hydrate forms collapse into one unit (`step2/build_units.py`, `step2/units.csv`).
- Matching reuses `overlap.py` (Step 1) and `cellcheck/cellline_match.py`; nothing was re-matched.
- The "measured at both JUMP doses" filter applies only to T3.
- **Exclusions:** auranofin (its JUMP match lacks the gold), and estrone sulfate → estrone.
- **Stereo/racemate pairs** are flagged: 7 in T1, 3 in T3. EU-OPENSCREEN has no compound names, so its pairs cannot be detected. Flagged pairs stay in the main analysis; the exact-achiral rerun is the sensitivity check.
- **Two T1 name-bridge prodrugs are kept:** ixazomib citrate and methylprednisolone succinate. Ixazomib and ixazomib citrate share the same LINCS wells, which gives one morphology pair at cosine 1.

**T2 gate: passed.** HepG2/C3A has a median of 937 treated cells per drug at 5 µM, and 99% of conditions have at least 100 cells. The threshold is 100 cells because that pools thousands of UMIs per well-expressed gene; the Tahoe-wide median is about 1,260. See `step2/gate_t2.py` and `gate_t2.json`.

**Gallery accessions** were resolved from the gallery's own S3 index (`step2/gallery_index.py`), not hardcoded.

**Transcriptome side (Tahoe-100M)**
- **Source:** the official pseudobulk DE tables on Hugging Face, streamed by range request (14.9 GB); no DE file was downloaded whole.
- **Unit of measurement:** one well holds one drug × dose × plate, with all 50 lines pooled. Each row is that well's log2FC for one line against the same plate's DMSO in the same line.
- **Aggregation:** wells are averaged per drug × dose × line, then a cross-line consensus is taken.
- **Cell threshold:** conditions with fewer than 100 treated cells were dropped. This removed NCI-H596 and NCI-H2122 entirely, leaving 48 lines.
- **Genes:** 5,000 well-expressed genes.
- **Replicate scheme unconfirmed:** what DESeq used as replicates ("pseudoreplicates") could not be confirmed from public documentation, so only fold changes are used, never Tahoe p-values.
- **Cells sharing a well are not independent.** All 50 lines of a drug × dose share one well. After centring, the same drug in different lines of the *same* well correlates at r ≈ 0.31. In a *different* well it correlates at r ≈ 0.20. Any within-well cross-line statistic is therefore inflated.

**Morphology side.** Each dataset was processed strictly within itself, and features were never merged across tiers. Details are in `step2/{lincs,euos,jump,cpjump1}_notes.md` and `step2/MORPH_SPEC.md`.
1. Per-plate robust z-scores against the negative controls.
2. Feature selection: blocklist, outliers, low variance, then |r| > 0.9 pruning.
3. PCA-whitening fitted on the control wells.
4. Consensus profile: median of replicate wells, then mean across sites or sources.

## Headline: agreement per tier

**How the numbers are computed**
- **Statistic:** Spearman Mantel r between the two cosine-similarity matrices.
- **Within-MoA null:** labels are permuted only within each MoA class, with singletons fixed; 9,999 permutations, one-sided.
- **Plain null:** labels are permuted freely; 1,999 permutations.
- **CI:** 95% bootstrap over units (1,000 resamples), with duplicate self-pairs dropped.
- **Ceiling:** √(morphology split-half reliability, Spearman–Brown corrected × Tahoe split-plate reliability).

| | T1 LINCS A549 | T2 EU-OS HepG2 | T3 JUMP U2OS |
|---|---|---|---|
| **PRIMARY: annotated, within-MoA null** | **r = 0.054** [−0.021, 0.133]; null mean 0.043; **p = 0.25** | **r = 0.058** [−0.015, 0.152]; null mean 0.046; **p = 0.31** | **r = 0.092** [0.016, 0.173]; null mean 0.044; **p = 0.035** |
| Annotated, plain null | 0.054, p = 0.015 | 0.058, p = 0.014 | 0.092, p = 0.004 |
| All units, plain null | 0.015 [−0.015, 0.050], p = 0.11 | 0.007 [−0.031, 0.049], p = 0.34 | 0.044 [0.002, 0.088], p = 0.008 |
| *Plate-centred*, annotated, within-MoA null | 0.041 [−0.031, 0.117], p = 0.020 | 0.058 [−0.012, 0.136], p = 0.017 | 0.020, p = 0.14 |
| *Plate-centred*, all units, plain null | 0.020 [−0.011, 0.058], p = 0.039 | 0.048 [0.010, 0.091], p ≤ 0.0005 | 0.037 [0.002, 0.079], p = 0.003 |
| *Plate-centred*, pairs sharing no plate | 0.029, p = 0.013 | 0.042, p = 0.002 | 0.026, p = 0.026 |
| *Centred*, two-plate units only (n) | 0.066 (61), p = 0.012 | 0.098 (56), p = 0.002 | 0.136 (45), p = 0.003 |
| *Centred*, one-plate units only (n) | 0.009 (126), p = 0.25 | 0.027 (115), p = 0.05 | 0.019 (130), p = 0.10 |
| Ceiling on the two-plate units (same units, centred) | 0.45 | 0.46 | 0.39 |
| Ceiling, full set (assumes single-plate reliability transfers), raw / centred | 0.57 / 0.43 | 0.72 / 0.52 | 0.50 / 0.45 |
| Sensitivity: exact-achiral only, all units (annotated) | −0.025 (87 / 41); within-MoA p = 0.77 | −0.033 (83 / 41); p = 0.95 | 0.027 (88 / 42); p = 0.70 |

`step2/figures/agreement_per_tier.png` shows the primary r with its CI against each tier's within-MoA null distribution. All null arrays are in `step2/results/nulls.npz`.

**How to read the table**
- **Why plate-centring is needed.** In the uncentred Tahoe matrix, same-plate pairs have mean cosine 0.26–0.39, against 0.07–0.16 across plates, and 26–29% of all pairs share a plate. Same-MoA pairs share a plate 40% of the time, against 30% for other pairs. So the uncentred matrix is partly a plate matrix, and that plate structure partly lines up with MoA.
- **What centring does.** It subtracts each plate × line × dose's mean signature over all ~80 drugs on that plate. It was adopted after the primary test was declared, so it is labelled secondary and both versions are reported.
- **The two error estimates disagree in places.** Several bootstrap CIs include 0 while the permutation p is below 0.05 (e.g. T1 centred within-MoA, and T1 and T2 on two-plate units). The bootstrap is the more conservative of the two; give both weight.
- **The exact-achiral subsets show nothing,** but with only 41–42 annotated units they have little power.
- **Within-MoA tests could not run on the two-plate subsets,** which have only 1–2 MoA classes with ≥3 members. Those figures are all-unit, plain-null results.

### MoA retrieval in each space separately

Retrieval is leave-one-compound-out nearest-neighbour on annotated compounds. Queries come from classes with ≥3 members, and the null permutes labels.

| Tier | Morphology: top-1 / mAP (null) | Transcriptome raw: top-1 / mAP (null) | Transcriptome plate-centred: top-1 / mAP |
|---|---|---|---|
| T1 | 0.43 / 0.34 (0.07 / 0.12) | 0.21 / 0.19 (0.07 / 0.12) | 0.24 / 0.19 |
| T2 | 0.47 / 0.37 (0.07 / 0.12) | 0.22 / 0.19 (0.07 / 0.12) | 0.15 / 0.19 |
| T3 | 0.27 / 0.26 (0.06 / 0.11) | 0.17 / 0.19 (0.06 / 0.11) | 0.13 / 0.17 |

- **mAP:** p ≤ 0.002 throughout.
- **Top-1:** p ≤ 0.002, except raw T3 transcriptome (p = 0.007), centred T2 transcriptome (p = 0.06) and centred T3 transcriptome (p = 0.04).

### Reliability

| | T1 LINCS | T2 EU-OS | T3 JUMP |
|---|---|---|---|
| Morphology similarity-matrix split-half reliability (Spearman–Brown) | 0.82 at 3.33 µM, 0.96 at 10 µM, 0.52 at 0.04 µM | 0.95 | 0.77 at 10 µM, 0.64 at 0.625 µM |
| Position-matched replicate reproducibility (from the pipeline agents) | 64–82% active at 3.33–10 µM (`rep_pos`) | 73–86% per site (`vs_pos`) | 74% at 10 µM, 55% at 0.625 µM (`strict`) |

**The morphology split-half figures are upper bounds.** Replicate wells share their well position in 100% of LINCS and EU-OPENSCREEN pairs, because each layout repeats across plates. Position effects therefore survive into both halves. For the same reason, the simple DMSO-based "reproducible" fractions saturate at 85–99%, so the position-matched figures above are the ones to trust.

**Tahoe similarity-matrix reliability** was measured on the 93 units observed on two plates at 5 µM, each copy with independent wells and controls:

| Tahoe signature (5 µM) | Raw | Plate-centred |
|---|---|---|
| Consensus | 0.33 | 0.26 |
| A549 | 0.40 | 0.23 |
| HepG2/C3A | 0.53 | 0.28 |

- **At 0.05 and 0.5 µM, reliability is about 0** (about 30 units each).
- **Only 4–8% of conditions are "active"** (robust z > 3).
- **Median signature magnitude does not rise with dose:** 23 / 22 / 23.

## Dose-response (T1 only)

**Only three of the six LINCS doses are within 1.5× of a Tahoe dose:** 0.041→0.05, 0.37→0.5 and 3.33→5 µM. The other three are 2.2–2.4× away; the brief said all six.

**Agreement does not change with concentration.** `step2/figures/T1_dose_response.png` pairs each LINCS dose with the nearest Tahoe dose.
- The all-unit Mantel r is between −0.013 and +0.021 at every dose, and every bootstrap CI includes 0.
- The full 6 × 3 grid (`DOSE_grid` rows in `T1_LINCS_A549_results.csv`) is also flat, from −0.024 to +0.021.
- Two annotated within-MoA tests reached nominal p ≈ 0.01: at 0.37 µM (r = 0.055) and 1.11 µM (r = 0.117 [0.02, 0.21]). That is 2 of about 12 dose tests, consistent with chance.
- Tahoe has no measurable reliability at 0.05 or 0.5 µM, so the low-dose half of this curve cannot be interpreted.

**Does one modality respond at lower concentrations than the other? The data cannot answer this.**

| Activity calls | Active at matched levels (morphology / Tahoe) | Morphology first | RNA first | Same | Neither |
|---|---|---|---|---|---|
| As specified: morphology position-matched replicate test; Tahoe norm z > 3 | 47–64% / 4–9% | 134 | 2 | 6 | 45 |
| Strictness matched: Tahoe z residualized on cell count; morphology threshold set to Tahoe's top-dose rate | 5–8% / 3–6% | 12 | 11 | 3 | 161 |
| Same, block-reliable LINCS units only | | 7 | 8 | 2 | 87 |

- **The specified calls are not comparable.** The morphology call fires about 10× as often as the transcriptome call, so "morphology first" in 134 cases reflects detection sensitivity, not biology.
- **The Tahoe norm tracks treated-cell count** (Spearman −0.70 within A549), so a Tahoe "active" call partly means cytotoxicity or a poor well.
- **At matched strictness it is a coin flip** (12 vs 11), with only 11 compounds active in both modalities at any matched dose.
- **Further confounds:** a 48 h vs 24 h exposure difference, and LINCS position confounds. Only 104 of 187 units sit in plate blocks where activity calls are reliable (`lincs_notes.md`).

## CPJUMP1 calibration: how much does cell line alone move morphology?

The same 306 compounds were imaged at 5 µM in A549 and in U2OS under one protocol.

| | 48 h | 24 h |
|---|---|---|
| r between the A549 and U2OS similarity matrices | 0.29 [0.25, 0.32] | 0.34 |
| Attenuation factor (disattenuated, Spearman–Brown) | **0.39** [0.35, 0.43] | 0.49 |
| Attenuation, compounds active beyond well position in both lines | 0.44 [0.38, 0.52] | – |
| Attenuation, strongest third of compounds | 0.52 [0.42, 0.62] | – |

**Changing the cell line keeps about 40–50% of the reliable similarity structure.** The all-compound value is contaminated by well position, since every plate uses the same platemap. Read the result as a bracket, 0.39–0.52.

**This calibration cannot explain the low agreement in T3.** If the cell-line mismatch were what holds T3 back, the cell-line-matched T1 should agree about twice as strongly. It does not:

| Compound set | T1 r | T3 r |
|---|---|---|
| Full sets | 0.015 | 0.044 |
| Same 92 compounds (T1∩T3) | −0.008 | 0.041 |

So the attenuation factor is context, not a correction to apply to T3. The gap between modalities, together with differences in exposure and protocol, matters more than the cell line.

## Identical-composition comparisons (intersections)

Every tier uses 10 µM morphology against Tahoe 5 µM, with its own line's signature (consensus for T3). Values are all-unit r under the plain null unless noted. The within-MoA test requires at least 3 MoA classes with 3 or more members each.

| Intersection | Units (annotated, classes ≥3) | T1 r | T2 r | T3 r |
|---|---|---|---|---|
| T1∩T2∩T3 | 58 (25, 2 classes); within-MoA not testable | −0.001 [−0.09, 0.09] | −0.012 [−0.08, 0.08] | 0.074 [−0.02, 0.18], p = 0.014 |
| T1∩T3 | 92 (39, 5) | −0.008 | – | 0.041, p = 0.054 |
| T1∩T2 | 110 (51, 6) | −0.018 | 0.005 | – |
| T2∩T3 | 90 (42, 7) | – | 0.020 | 0.067, p = 0.006 (within-MoA 0.135, p = 0.026) |

On identical compounds, T3 consistently comes out highest and T1 lowest, which again says the cell line is not what limits agreement. These are secondary, uncorrected tests and are not findings on their own.

## Secondary question: does cross-line consistency predict agreement?

For each compound, row-wise agreement is the Spearman correlation between its similarity row in each space. This was correlated with Tahoe cross-line consistency at 5 µM.

| Tier | Raw consistency, ρ [95% CI] | Plate-centred consistency, ρ [95% CI] | Partial ρ controlling for magnitude (centred) |
|---|---|---|---|
| T1 | −0.011 [−0.17, 0.13] | 0.102 [−0.05, 0.26] | 0.079 [−0.07, 0.21] |
| T2 | −0.033 [−0.19, 0.12] | −0.029 [−0.17, 0.12] | −0.030 [−0.18, 0.11] |
| T3 | −0.063 [−0.21, 0.08] | 0.016 [−0.14, 0.17] | 0.018 [−0.14, 0.17] |

- **No relationship is detectable;** every CI spans 0 (`step2/figures/consistency_vs_agreement.png`).
- **The test is effectively untestable here:**
  - with overall agreement around 0.02–0.05, per-compound row agreement is mostly noise;
  - within-well consistency is inflated by shared-well effects.
- **The two-plate result answers this question better.** Agreement is higher where the transcriptome is measured more reliably (verdict point 3). It points the way the question expected.

## Other checks

- **Signature choice (3×3, `SIG_CHOICE` rows).** No tier gets better with its matched line's signature than with the consensus or the other line.
- **T3 low dose** (JUMP 0.625 µM vs Tahoe 0.5 µM): all-unit r = 0.006, within-MoA p = 0.65. This is expected, since Tahoe reliability is about 0 at 0.5 µM.

## Caveats, stated plainly

1. **Multiple testing.** The headline rests on 3 pre-registered tests. Everything else in `all_results.csv` is secondary and uncorrected: 262 rows carry a p-value, and 77 of them are below 0.05.
   - **Most hits are plain-permutation tests.** Those succeed by design whenever both spaces share MoA class structure.
   - **14 of 71 within-MoA tests are nominal hits,** against about 3.5 expected by chance. These tests are strongly correlated (overlapping compounds, neighbouring doses), so this excess is supportive rather than independent evidence of a small agreement beyond MoA class.
   - **Single hits are not findings,** for example LINCS 3.33 µM paired with Tahoe 0.5 µM (r = 0.125, p = 0.004), a dose pair that isn't even matched.
2. **Post hoc analyses.** Plate-centring and the two-plate vs one-plate split were added after the primary test was declared, in response to confounds found in the data. Both are labelled secondary.
3. **The two-plate vs one-plate split is not randomized.** Which drugs Tahoe placed on two plates may differ in ways beyond noise. Their centred magnitude is lower, and their activity rate is similar (T1, T2) or higher (T3).
4. **Shared wells.** Tahoe cells within a well, across all 50 lines, are not independent. The DESeq replicate scheme is unconfirmed, so no Tahoe p-values were used.
5. **Exposure time.** LINCS and JUMP are 48 h; Tahoe and EU-OPENSCREEN are 24 h. In CPJUMP1 the attenuation factor differs between 24 h and 48 h (0.49 vs 0.39).
6. **Cell lines.** HepG2 (EU-OPENSCREEN) is the parent of Tahoe's HepG2/C3A subline. U2OS is not a Tahoe line.
7. **Dose.**
   - Tahoe 5 µM is paired with 10 µM morphology in T2 and T3, and with 3.33 µM in T1.
   - Below 5 µM the Tahoe transcriptome has little reliable signal.
8. **Well-position confounds in morphology.**
   - Replicate wells share a position in all LINCS and EU-OPENSCREEN pairs, and in about 47% (high dose) and 100% (low dose) of JUMP compound wells.
   - Morphology reliabilities are therefore upper bounds.
   - Simple DMSO-based activity calls are saturated in every dataset.
9. **EU-OPENSCREEN sites disagree.** Replicate cosine is 0.37 within a site and 0.15 across sites.
10. **Stereochemistry.** All datasets are stereo-flat, so enantiomers and epimers collapse into one unit, e.g. doxorubicin/epirubicin and quinidine/quinine.
11. **MoA labels.** These are Tahoe `moa-fine`, about 50% of which are "unclear". Within-MoA tests therefore use about 82–85 annotated units per tier.

## Data, bytes and disk (as printed by `compare_tiers.py`)

| Pipeline | Bytes downloaded (final runs) | Peak local disk |
|---|---|---|
| Tahoe-100M DE (range reads) | 14.94 GB, of which footers 5.9 GB | 1.52 GB cache, plus 0.77 GB per-line table |
| JUMP U2OS (416 plates) | 6.10 GB | 2.0 GB plate cache, deleted afterwards |
| EU-OPENSCREEN (range reads into the Zenodo zip) | 1.08 GB | 0.89 GB, deleted afterwards |
| LINCS (136 plates) | 0.28 GB | 0.03 GB |
| CPJUMP1 (16 plates) | 0.11 GB | 0.25 GB, deleted afterwards |
| **Total** | **22.5 GB** | **7.49 GB upper bound** |

- **Development reruns** add roughly 1.5–2 GB of downloads on top of the 22.5 GB, which counts final runs only.
- **The peak-disk figure is an upper bound:** it sums every pipeline's peak plus the current project size, but the pipelines did not all peak at the same time.
- **The project directory is now 2.8 GB,** of which `step2/tahoe_cache/` is 1.5 GB. Keeping that cache makes Tahoe reruns cheap; delete it to reclaim the space.

## Reproduce

```bash
py -3.13 -m pip install -r requirements.txt               # pinned versions (Python 3.13)
py -3.13 -I -X utf8 overlap.py                            # Step 1: Tahoe x JUMP matching -> overlap.csv
py -3.13 -I -X utf8 cellcheck/obs_cell_lines.py           # which 50 Tahoe lines were screened
py -3.13 -I -X utf8 cellcheck/cellline_match.py           # LINCS / EU-OS / CPJUMP1 matching -> cellline_matches.csv
py -3.13 -I -X utf8 step2/build_units.py                  # T3 working set + compound units per tier
py -3.13 -I -X utf8 step2/tahoe_signatures.py             # Tahoe (stages pilot/fetch2/index/build/perplate; resumes from cache)
py -3.13 -I -X utf8 step2/gate_t2.py                      # T2 cell-count gate (reads the Tahoe index)
py -3.13 -I -X utf8 step2/jump_profiles.py                # T3 morphology
py -3.13 -I -X utf8 step2/morph_lincs.py                  # T1 morphology
py -3.13 -I -X utf8 step2/morph_euos.py                   # T2 morphology
py -3.13 -I -X utf8 step2/morph_cpjump1.py                # CPJUMP1 calibration
py -3.13 -I -X utf8 step2/compare_tiers.py                # all comparisons, figures, accounting (~25 min)
```

Every pipeline raises an error when a source is unreachable; no tier is skipped silently. Downloaded metadata is cached under `step2/*_cache/`.

**Outputs**
- Results tables: `step2/results/all_results.csv`, plus per-tier `T1_LINCS_A549_results.csv`, `T2_EUOS_HepG2_results.csv`, `T3_JUMP_U2OS_results.csv` and `TAHOE_results.csv`.
- Supporting tables: `CAL_cpjump1_calibration.csv`, `*_line_ranking.csv`, `T1_lowest_active_per_unit*.csv`, `consistency_vs_agreement.csv`, `nulls.npz` and `accounting.json`.
- Figures (`step2/figures/`): `agreement_per_tier.png`, `T1_dose_response.png`, `consistency_vs_agreement.png`.
