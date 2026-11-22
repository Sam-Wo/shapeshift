# JUMP morphology profiles: notes (written by step2/jump_profiles.py)

## Data access and format quirks
- Path pattern verified: `https://cellpainting-gallery.s3.amazonaws.com/cpg0016-jump/{source}/workspace/profiles/{batch}/{plate}/{plate}.parquet` (anonymous HTTPS GET).
- Each file holds only `Metadata_Source`, `Metadata_Plate`, `Metadata_Well` plus 4,762 float32 CellProfiler features (source_15: 3,805; it lacks 957 features, mostly texture/correlation/radial-distribution and Image_* ones). Well IDs are zero-padded (`A01`), matching well.csv; 1536-well plates (source_1, source_9) use `AA01`-style IDs for rows beyond Z.
- Profile files can contain wells absent from well.csv; we inner-join on source/plate/well.
- Some features are mostly NaN on some plates (e.g. `Cells_AreaShape_FormFactor`); handled by the 1% NaN rule.
- Plates downloaded: 416 of 416 needed; failures: 0.
- Bandwidth: 6.10 GB of plate files (+~0.17 GB of manual probes while developing). The filtered per-plate cache `step2/jump_cache/` (2.00 GB) is deleted at the end of the run unless `--keep-cache` is passed; a re-run then re-downloads. Per-plate download log: `step2/jump_download_manifest.csv`.
- Sum of per-plate download times 3002 s over 16 threads.
- Last full fresh download: 416 plates, 6.10 GB in 203 s wall (30 MB/s).
- Wall time of this run: 4.2 min (processing from the cache when it was present).

## Parameters
- NaN/inf rule: drop if > 1% of wells; constant features dropped; remaining NaN -> plate DMSO median.
- Robust z: (x - median_DMSO_plate) / max(1.4826*MAD_DMSO_plate, floor); floor = 0.01 x median-over-plates of 1.4826*MAD (min 1e-12); clip +-50.
- Selection: pycytominer blocklist (55), outlier features with pre-clip max |z| > 500, variance threshold (freq_cut 0.05, unique_cut 0.01), greedy |r| > 0.9 pruning (features visited in ascending total |r|, kept if |r| <= 0.9 with all kept ones).
- Sphering: PCA-whitening fit on DMSO wells only, centred on DMSO mean, components covering 99% of DMSO variance, scaled by 1/sqrt(eigenvalue + 3.0). eps was chosen from the grid (0.1, 1.0, 3.0, 10.0) by the replicate-vs-other-compound AUC over different-position pairs (position-robust), requiring it to match or beat the unsphered selected features at both doses (eps 0.1 was worse than no sphering). This is a mild use of the replicate signal for tuning.
- Consensus: HIGH = median of wells within each source, then mean across sources; LOW = median of wells.
- Activity: mean pairwise cosine among a compound's replicate wells (sphered) vs a plate-matched DMSO null (2000 draws: for every replicate well one DMSO well from the same plate, without replacement within plate); p = (1 + #null >= obs)/(1 + 2000); active = p < 0.05. BH q is also reported.
- The primary DMSO test is close to saturated: other-compound pairs already have mean cosine ~0.03-0.04 vs ~0 for DMSO pairs, so almost any set of compound wells beats DMSO. Stricter flags are therefore given alongside `active` (the caller chooses):
  - `p_diffpos` / `active_diffpos`: replicate score over pairs at different positions only (not same source+well), against the same plate-matched DMSO null with the same pair mask (controls positional similarity).
  - `p_cpdnull` / `active_vs_cpd`: all pairs vs a same-source other-compound null (each replicate well replaced by a well of a different working-set compound of the same source, at most one well per compound; whole-dataset pool for sources with too few compounds). Conservative: the null compounds are mostly bioactive drugs too.
  - `p_strict` / `active_strict`: different-position pairs vs the other-compound null (`strict_basis` = allpairs when a compound has no different-position pairs).

## HIGH (10.0 uM)
- Sources: source_1, source_10, source_11, source_15, source_2, source_3, source_5, source_6, source_8, source_9
- Plates 306; wells 14063 = 1477 compound + 12586 DMSO; compounds 175/175.
- Expected wells missing from profile files: 0 (compound wells 0) 
- Features: union 4762 -> intersection 3805 (missing by source: {'source_15': 957}) -> NaN>1% -0, constant -17 = 3788 -> blocklist -36, outlier -567, low variance -25, correlated -2274 = **886 selected** (185 of them Image_* features).
- Imputed values: 237; fraction of plate x feature MADs at floor: 0.0221.
- Sphering: **532 dimensions** (eigenvalues 278.02 ... 0.130); replicate AUC by eps: {0.1: {'all': 0.7132, 'diff_pos': 0.7058}, 1.0: {'all': 0.7308, 'diff_pos': 0.7241}, 3.0: {'all': 0.7391, 'diff_pos': 0.7328}, 10.0: {'all': 0.7437, 'diff_pos': 0.7375}}.
- Replicate check (selected features, unsphered): replicate-pair vs other-compound-pair cosine AUC 0.724 (cross-plate 0.724; different-position pairs 0.718); mean cosine replicate 0.283 (same position 0.481, different 0.272) vs other compound 0.075 (6684 replicate pairs, 346 same position, 26 same plate).
- Replicate check (sphered): replicate-pair vs other-compound-pair cosine AUC 0.739 (cross-plate 0.739; different-position pairs 0.733); mean cosine replicate 0.168 (same position 0.358, different 0.158) vs other compound 0.041 (6684 replicate pairs, 346 same position, 26 same plate).
- Active (p < 0.05): **168/175 = 96.0%**; BH q < 0.05: 96.0%; replicate cosine above null mean for 98.9%; median null percentile 100.0.
- Position effect: DMSO pairs on different plates of the same source have mean cosine 0.068 at the same well vs 0.002 at different wells.
- Stricter calls: active_diffpos 96.0% (of the 175 compounds with different-position replicate pairs), active_vs_cpd 72.6%, **active_strict 73.7%**; compounds without different-position pairs: 0; other-compound null used the whole-dataset pool for 1.1% of wells.
- Time: 175 s.

## LOW (0.625 uM)
- Sources: source_7
- Plates 110; wells 5865 = 1290 compound + 4575 DMSO; compounds 175/175.
- Expected wells missing from profile files: 0 (compound wells 0) 
- Features: union 4762 -> intersection 4762 (missing by source: none) -> NaN>1% -2, constant -30 = 4730 -> blocklist -55, outlier -74, low variance -8, correlated -3567 = **1026 selected** (241 of them Image_* features).
- Imputed values: 0; fraction of plate x feature MADs at floor: 0.0067.
- Sphering: **599 dimensions** (eigenvalues 309.66 ... 0.088); replicate AUC by eps: {0.1: {'all': 0.619, 'diff_pos': 0.5288}, 1.0: {'all': 0.637, 'diff_pos': 0.5355}, 3.0: {'all': 0.6466, 'diff_pos': 0.5379}, 10.0: {'all': 0.6528, 'diff_pos': 0.5359}}.
- Replicate check (selected features, unsphered): replicate-pair vs other-compound-pair cosine AUC 0.644 (cross-plate 0.645; different-position pairs 0.531); mean cosine replicate 0.219 (same position 0.348, different 0.114) vs other compound 0.081 (5455 replicate pairs, 2466 same position, 45 same plate).
- Replicate check (sphered): replicate-pair vs other-compound-pair cosine AUC 0.647 (cross-plate 0.648; different-position pairs 0.538); mean cosine replicate 0.111 (same position 0.179, different 0.055) vs other compound 0.026 (5455 replicate pairs, 2466 same position, 45 same plate).
- Active (p < 0.05): **166/175 = 94.9%**; BH q < 0.05: 94.9%; replicate cosine above null mean for 99.4%; median null percentile 100.0.
- Position effect: DMSO pairs on different plates of the same source have mean cosine 0.032 at the same well vs -0.001 at different wells.
- Stricter calls: active_diffpos 59.1% (of the 66 compounds with different-position replicate pairs), active_vs_cpd 60.0%, **active_strict 55.4%**; compounds without different-position pairs: 109; other-compound null used the whole-dataset pool for 0.0% of wells.
- Time: 78 s.

## Coverage
- Compounds missing at HIGH: none; at LOW: none.
- Caveat: 47% (HIGH) and 100% (LOW) of compound wells share their well position with another replicate of the same compound (same source), so part of the replicate similarity can be positional (see the DMSO position effect and the `active_diffpos` / `active_strict` flags above). `well` is kept in the wells files.

## Outputs
- `morph_profiles.parquet` (root): both doses, f0..fN float32, NaN in columns beyond a dose's dimensionality. Also `morph_profiles_high.parquet` and `morph_profiles_low.parquet` (root) with only that dose's columns. f-columns are NOT comparable between doses (different spaces).
- `step2/morph_activity.csv`, `step2/morph_wells_{high,low}.parquet`, `step2/jump_run.json` (full run info incl. selected feature names).

## Post-run checks (scratch diagnostics, not part of the pipeline)
- HIGH "same position" = same source + well is the right definition: DMSO pairs at the same well ID but different sources have mean cosine 0.005 (vs -0.001), against 0.070 within a source.
- Position leakage into consensus profiles, for compound pairs that share no plate: HIGH, pairs sharing a (source, well) position have cosine 0.156 (n=119) vs 0.148 (n=13,566), negligible; LOW 0.135 (n=86) vs 0.080 (n=13,631), modest and few pairs. Positions per compound can be taken from the wells files if the comparison wants to exclude these pairs.
- LOW `active_strict` by basis: 22.7% of the 66 compounds with different-position replicate pairs (position-controlled) vs 75.2% of the 109 whose replicates all sit at one position. For those 109, the "allpairs" basis cannot separate compound effect from position, so treat their LOW activity as unresolved.
- 1536-well IDs confirmed as `AA01`... in well.csv and profiles.
