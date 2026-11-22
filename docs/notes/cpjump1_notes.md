# CPJUMP1 calibration notes (written by step2/morph_cpjump1.py)

## Sources
- Gallery accession `cpg0000-jump-pilot` (resolved with `gallery_index.resolve('cpjump1')`), batch `2020_11_04_CPJUMP1`.
- Profiles: `cpg0000-jump-pilot/source_4/workspace/profiles/2020_11_04_CPJUMP1/<plate>/<plate>.csv.gz` (verified; raw well-level CellProfiler aggregates, only `Metadata_Plate`/`Metadata_Well` as metadata, 384 wells, 5,792 features incl. brightfield channels). Other files in the folder (`_augmented`, `_normalized*`, `_spherized`) are not used.
- Metadata (cached in `step2/cpjump1_cache/`, sha256 prefix): `experiment-metadata.tsv_sha256` 86d845dc152cbc35; `JUMP-Target-1_compound_metadata_additional_annotations.tsv_sha256` 23399e009100e268; `JUMP-Target-1_compound_metadata.tsv_sha256` 8a17b7c5136ca9c4; `JUMP-Target-1_compound_platemap.tsv_sha256` b22a82080fe75a76
  - experiment metadata and additional annotations: GitHub jump-cellpainting/2024_Chandrasekaran_NatureMethods (`benchmark/output/experiment-metadata.tsv`, `benchmark/input/JUMP-Target-1_compound_metadata_additional_annotations.tsv`, main);
  - compound metadata and platemap: GitHub jump-cellpainting/JUMP-Target (master).
- `target` (single designed target gene per compound; 146 genes x 2 compounds + 14 singletons) comes from the NatureMethods additional annotations: the JUMP-Target compound_metadata.tsv itself only has `target_list` (multi-gene). The compound key for everything is the full `broad_sample`; units map through `units.csv` `dataset_ids` (not InChIKey: the two metadata files disagree on some InChIKeys, e.g. amlodipine).

## Plates (one protocol)
- Filter = the NatureMethods benchmark's own: batch == 2020_11_04_CPJUMP1 (not `_DL`), Perturbation == compound, Density == 100, Antibiotics absent, Cas9 compound plates dropped. ORF/CRISPR plates excluded.
- A549 24 h: BR00116991, BR00116992, BR00116993, BR00116994 (anomaly flags: WGA, WGA, WGA, WGA; images: 27648, 27640, 27352, 27576)
- A549 48 h: BR00117015, BR00117016, BR00117017, BR00117019 (anomaly flags: none, none, none, none; images: 49144, 49152, 49144, 49152)
- U2OS 24 h: BR00116995, BR00117024, BR00117025, BR00117026 (anomaly flags: WGA, WGA, WGA, WGA; images: 27648, 27648, 27648, 27648)
- U2OS 48 h: BR00117010, BR00117011, BR00117012, BR00117013 (anomaly flags: Mitotracker, none, none, Mitotracker; images: 27648, 27648, 27648, 27648)
- Excluded compound plates of the batch (density/Cas9): BR00117008 (A549 48h, density 80, Parental), BR00117009 (A549 48h, density 80, Parental), BR00117050 (A549 48h, density 100, Cas9), BR00117051 (A549 48h, density 100, Cas9), BR00117052 (A549 48h, density 100, Cas9), BR00117053 (A549 48h, density 100, Cas9), BR00117054 (A549 48h, density 120, Parental), BR00117055 (A549 48h, density 120, Parental)
- Caveats: all 24 h plates carry the `WGA` anomaly flag (both lines, balanced); 2 of 4 U2OS 48 h plates carry `Mitotracker`; the A549 48 h plates were imaged with 16 instead of 9 sites per well (~49k images). Kept as published by the benchmark.
- Every plate uses the same JUMP-Target-1 platemap (306 compounds, 292 once and 14 twice per plate, 64 DMSO wells), so each compound sits at the same well in all replicates and in both cell lines. Positional effects can therefore inflate replicate mAP, the split-half ceiling and r_cross alike; see the diagnostics below.

## Processing (per cell line x time point, against its own DMSO)
- Raw columns: one intersection per cell line across both time points. Feature-selection scope: `line_time` (constant `SELECTION_SCOPE`; 'line' would select once per cell line). Sphering is always per line x time.
- Pre-filter (additions to the spec): brightfield-channel features dropped (Brightfield/HighZBF/LowZBF, to keep the standard 5-channel Cell Painting space of the other tiers) and image-coordinate features (AreaShape_Center_*, BoundingBoxMin/Max, Location_*).
- NaN/inf > 1% or constant dropped; NaN -> plate DMSO median; robust z vs plate DMSO with MAD floor (0.01 x median plate scale); clip +-50; pycytominer 55-feature blocklist; outlier features (pre-clip max|z| > 500); variance threshold (freq_cut 0.05, unique_cut 0.01) evaluated on the UNCLIPPED z (on clipped z, strong-responder features where many treated wells sit at +-50 were flagged as low-variance because the clip value became the modal value: 426 features in U2OS 48 h, 187 in U2OS 24 h); greedy |r| > 0.9 pruning.
- Sphering: PCA whitening fit on DMSO wells only, 99% DMSO variance, eigenvalue + 3.0 (fixed a priori, same as the JUMP tier; not tuned on any replicate or cross-line statistic).
- Consensus: median of replicate wells (single site, source_4). DMSO reference row = median of all DMSO wells.
- Replicate mAP: each replicate queries its other replicates + all DMSO wells on the same plates by cosine (sphered); null = 10000 random DMSO sets with the same per-plate composition (queried against the remaining DMSO wells); p = (1+#null>=obs)/(1+N); reproducible / active_rep if p < 0.05. active_z: robust z of the consensus L2 norm vs all treated compounds > 3.

| dataset | plates | wells (DMSO) | raw int. | after prefilter | NaN/const drop | after | blocklist | outlier | low var | corr | selected | sphere dim | reproducible | active_z | active_pos | reproducible (matched) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A549 24 h | 4 | 1535 (256) | 5792 | 3399 (-2277 bf, -116 pos) | -0/-3 | 3396 | -55 | -60 | -0 | -2693 | 588 | 84 | 0.928 | 0.232 | 0.621 | 0.957 |
| A549 48 h | 4 | 1536 (256) | 5792 | 3399 (-2277 bf, -116 pos) | -0/-3 | 3396 | -55 | -16 | -0 | -2747 | 578 | 117 | 0.967 | 0.255 | 0.680 | 1.000 |
| U2OS 24 h | 4 | 1453 (242) | 5792 | 3399 (-2277 bf, -116 pos) | -0/-3 | 3396 | -55 | -233 | -0 | -2384 | 724 | 76 | 0.886 | 0.222 | 0.569 | 0.957 |
| U2OS 48 h | 4 | 1536 (256) | 5792 | 3399 (-2277 bf, -116 pos) | -0/-3 | 3396 | -55 | -9 | -0 | -2739 | 593 | 108 | 0.850 | 0.265 | 0.549 | 1.000 |

- Imputed values / floored plate x feature MAD fraction: A549 24h 0 / 0.0006; A549 48h 0 / 0.0; U2OS 24h 0 / 0.0006; U2OS 48h 0 / 0.0006
- Compounds: 306 treated (all 5 uM), 4 wells per line x time (8 for the 14 twice-plated), except where wells are absent from the gallery profile files: BR00116993: 1 wells (N01); BR00116995: 83 wells (A01, A02, A03, A04, A05, A06, ...). BR00116995 (U2OS 24 h) is a partial plate (83 of 384 wells missing in odd rows A-O, cols 1-14); kept, so those compounds have 3 U2OS 24 h replicates. Matched subset: 21 CAL_CPJUMP1 units = 23 broad_samples (dexamethasone and quinidine/quinine map 2 IDs to one unit; within-unit pairs are dropped from the matched Mantel vectors and bootstrap).

## Calibration (cpjump1_calibration.csv; 48 h primary)
- Mantel: Spearman over the upper triangle of the compound x compound cosine matrices of the consensus profiles (DMSO excluded; compounds present in both lines). p from 9999 joint row/column permutations of the U2OS matrix (on ranks), and from 9999 permutations within JUMP-Target `target` classes (singletons fixed). CI: 1000 paired bootstrap draws over compounds (self-pairs from duplicated compounds dropped); the same draw gives r_cross, both reliabilities and the ratios.
- Reliability: 100 random split halves per line x time (each compound's wells split independently, 2+2 or 4+4), Mantel Spearman between the half matrices; mean over splits; Spearman-Brown 2r/(1+r). Within each bootstrap draw 10 rotating splits are used. Reliabilities are recomputed on each compound set.
- Active set: replicate-mAP p < 0.05 in both lines at that time point. CAUTION: this call is position-confounded (all replicates share one well; 36-70% of DMSO well positions pass the same test, see diagnostics), so `active` is close to `all`.
- `beyond_position` set: consensus norm > 95th percentile of the 64 DMSO-position consensus norms in BOTH lines (`active_pos` in the activity file); a position-robust activity call.

| time | set | n | r_cross [95% CI] (perm p) | p within-target | rel_SB A549 | rel_SB U2OS | attenuation_SB | attenuation (half rel) | row agreement median |
|---|---|---|---|---|---|---|---|---|---|
| 24 h | all | 306 | 0.336 [0.298, 0.377] (p<=1e-04) | <=1e-04 | 0.743 [0.711, 0.772] | 0.642 [0.597, 0.682] | 0.487 [0.432, 0.545] | 0.636 [0.558, 0.718] | 0.328 [0.306, 0.353] |
| 24 h | active | 258 | 0.353 [0.308, 0.400] (p<=1e-04) | <=1e-04 | 0.801 [0.771, 0.832] | 0.722 [0.675, 0.767] | 0.464 [0.408, 0.523] | 0.574 [0.503, 0.652] | 0.346 [0.319, 0.373] |
| 24 h | beyond_position | 162 | 0.421 [0.349, 0.491] (p<=1e-04) | <=1e-04 | 0.910 [0.885, 0.930] | 0.876 [0.838, 0.906] | 0.472 [0.392, 0.547] | 0.523 [0.435, 0.603] | 0.413 [0.378, 0.451] |
| 24 h | strong_tercile | 102 | 0.471 [0.378, 0.548] (p<=1e-04) | <=1e-04 | 0.962 [0.945, 0.974] | 0.950 [0.930, 0.965] | 0.493 [0.397, 0.573] | 0.514 [0.416, 0.596] | 0.474 [0.433, 0.528] |
| 24 h | weak_tercile | 102 | 0.401 [0.332, 0.462] (p<=1e-04) | <=1e-04 | 0.602 [0.547, 0.649] | 0.408 [0.342, 0.471] | 0.809 [0.691, 0.927] | 1.207 [1.030, 1.408] | 0.416 [0.375, 0.452] |
| 24 h | matched | 23 | 0.336 [0.096, 0.548] (p<=1e-04) | not testable (2 swappable classes, min p 0.25) | 0.820 [0.653, 0.909] | 0.701 [0.488, 0.839] | 0.443 [0.143, 0.724] | 0.549 [0.185, 0.922] | 0.320 [0.136, 0.426] |
| 48 h | all | 306 | 0.287 [0.253, 0.320] (p<=1e-04) | <=1e-04 | 0.835 [0.809, 0.857] | 0.650 [0.609, 0.688] | 0.390 [0.347, 0.432] | 0.489 [0.432, 0.544] | 0.281 [0.260, 0.300] |
| 48 h | active | 257 | 0.314 [0.271, 0.355] (p<=1e-04) | <=1e-04 | 0.874 [0.851, 0.895] | 0.731 [0.693, 0.768] | 0.393 [0.342, 0.445] | 0.470 [0.408, 0.533] | 0.297 [0.281, 0.327] |
| 48 h | beyond_position | 162 | 0.406 [0.340, 0.476] (p<=1e-04) | <=1e-04 | 0.947 [0.932, 0.961] | 0.887 [0.860, 0.910] | 0.443 [0.375, 0.516] | 0.480 [0.407, 0.556] | 0.368 [0.333, 0.383] |
| 48 h | strong_tercile | 103 | 0.505 [0.408, 0.599] (p<=1e-04) | 0.0002 | 0.981 [0.970, 0.988] | 0.963 [0.948, 0.974] | 0.520 [0.419, 0.617] | 0.534 [0.431, 0.633] | 0.478 [0.431, 0.534] |
| 48 h | weak_tercile | 102 | 0.340 [0.275, 0.407] (p<=1e-04) | <=1e-04 | 0.710 [0.662, 0.756] | 0.434 [0.381, 0.490] | 0.612 [0.495, 0.731] | 0.870 [0.695, 1.052] | 0.356 [0.327, 0.376] |
| 48 h | matched | 23 | 0.517 [0.298, 0.669] (p<=1e-04) | not testable (2 swappable classes, min p 0.25) | 0.858 [0.738, 0.921] | 0.708 [0.533, 0.827] | 0.664 [0.425, 0.844] | 0.806 [0.546, 1.038] | 0.488 [0.411, 0.566] |

Diagnostics (positional confound):
- 24 h: DMSO-position cross-line Mantel 0.491 (p<=1e-04); A549 twice-plated replicate cosine same well 0.803 vs different well 0.758; U2OS twice-plated replicate cosine same well 0.750 vs different well 0.699; A549 DMSO positions passing the replicate test 0.453, DMSO-position/treated norm ratio 0.627; U2OS DMSO positions passing the replicate test 0.359, DMSO-position/treated norm ratio 0.644
- 48 h: DMSO-position cross-line Mantel 0.486 (p<=1e-04); A549 twice-plated replicate cosine same well 0.871 vs different well 0.831; U2OS twice-plated replicate cosine same well 0.755 vs different well 0.712; A549 DMSO positions passing the replicate test 0.703, DMSO-position/treated norm ratio 0.470; U2OS DMSO positions passing the replicate test 0.453, DMSO-position/treated norm ratio 0.690
- 24 h `active_pos` (beyond position): both lines 162, A549 only 28, U2OS only 12, neither 104. Sets requiring activity in both lines (beyond_position, strong_tercile) exclude line-specific activity, itself a cell-line effect, so they bracket the attenuation from above; the A549 48 h plates were imaged at 16 sites/well, so part of the asymmetry may be sensitivity.
- 48 h `active_pos` (beyond position): both lines 162, A549 only 46, U2OS only 6, neither 92. Sets requiring activity in both lines (beyond_position, strong_tercile) exclude line-specific activity, itself a cell-line effect, so they bracket the attenuation from above; the A549 48 h plates were imaged at 16 sites/well, so part of the asymmetry may be sensitivity.
- Strength strata (`weak_tercile`, `strong_tercile`): compounds ranked by the mean over lines of their consensus-norm percentile, fixed before any calibration statistic was looked at in those strata. If r_cross in the weak tercile drifts toward the DMSO-position Mantel, the all-compound number is position-contaminated.

## Outputs
- `step2/morph_cpjump1.parquet`: consensus per broad_sample x cell line x time (all 306 compounds + a DMSO reference row), keyed by `broad_id` (broad_sample prefix) and `unit` where matched; f0..fN float32, NaN beyond a dataset's dimensionality. f-columns are NOT comparable across the four datasets.
- `step2/morph_cpjump1_activity.csv`, `step2/morph_cpjump1_wells.parquet` (sphered wells), `step2/cpjump1_calibration.csv`, `step2/cpjump1_row_agreement.csv` (per-compound row Spearman), `step2/cpjump1_run.json` (selected feature names, per-plate download log, full log).

## Bytes, disk, time
- Downloaded this run: 111.6 MB (plates 111.4 MB newly downloaded, plus S3 listings and GitHub metadata). Development runs before this one: ~128 MB (a first run that downloaded all 16 plate files, ~113 MB, before stopping at a missing-well assertion; two probe downloads of one plate file, ~14 MB; listings and GitHub files). The final run started from an empty `step2/cpjump1_cache/`.
- Peak disk (cache + outputs): 250.1 MB. Plate cache `step2/cpjump1_cache/plates/` deleted at the end; metadata cache kept.
- Wall time: 1327 s.
- Packages: python 3.13.0, pandas 2.3.3, numpy 2.2.6, scipy 1.16.3, scikit-learn 1.8.0, pyarrow 23.0.1, fsspec 2026.2.0, requests 2.32.5.

## Deviations from MORPH_SPEC
- Outputs cover ALL 306 compounds, not only matched units (task-mandated: the calibration uses all compounds).
- Extra pre-filter: brightfield features and image-coordinate features dropped (counts above).
- Feature space: normalization/selection/sphering per line x time; the 'shared within a cell line' requirement is met at the raw-column level (one intersection per line). No statistic crosses time points.
- Variance threshold on unclipped z (see above); the JUMP tier ran it on clipped z.
- Plates with missing wells (BR00116995 partial, BR00116993 one well) kept; DMSO-count check relaxed.
- Within-unit pairs (dexamethasone x2, quinidine/quinine) are dropped from every Mantel vector, not only the matched set (2 of ~46k pairs in the all set).
- Single site/source, so the consensus is a plain median.
