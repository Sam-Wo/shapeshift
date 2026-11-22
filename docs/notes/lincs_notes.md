# T1 LINCS Cell Painting (A549, 48 h, 6 doses): morphology profiles

Script: `step2/morph_lincs.py` (`py -3.13 -I -X utf8 step2/morph_lincs.py`). Metadata cache: `step2/lincs_cache/` (barcode_platemap.csv, 28 platemaps, run_info.json).

## Sources
- Accession `cpg0004-lincs` (resolved via `gallery_index.resolve('lincs')`), batch `2016_04_01_a549_48hr_batch1` (the Batch2 pilot `2017_12_05_Batch2` is not used).
- Profiles: `cpg0004-lincs/broad/workspace/profiles/2016_04_01_a549_48hr_batch1/<plate>/<plate>_augmented.csv.gz` (well-level, non-normalized aggregated CellProfiler features + platemap metadata; CSV, gzip). Streamed into memory, never written to disk.
- Platemaps: `cpg0004-lincs/broad/workspace/metadata/platemaps/2016_04_01_a549_48hr_batch1/` (barcode_platemap.csv + platemap/*.txt); cross-checked well by well against the augmented `Metadata_broad_sample` and `Metadata_mmoles_per_liter` (all agree).
- Plates: 136 profiled of 140 in barcode_platemap; missing (no profiles in the gallery): {'SQ00015225': 'C-7161-01-LM6-012', 'SQ00015226': 'C-7161-01-LM6-013', 'SQ00015227': 'C-7161-01-LM6-013', 'SQ00015228': 'C-7161-01-LM6-013'}. Plates per platemap: {2: 1, 4: 1, 5: 26} (count of plates -> number of platemaps).
- Compound key: `units.csv` tier `T1_LINCS_A549` `dataset_ids` are full 22-char broad_sample IDs (e.g. BRD-K70401845-003-09-6); all match `Metadata_broad_sample` exactly, so wells are mapped by full ID (each prefix has a single batch ID in this batch). The 13-char prefix keys compounds for the `norm_z` reference and the other-compound null. IDs shared by two units (wells emitted under both): {'BRD-K78659596-001-03-9': ['MBOMYENWWXQSNW-UHFFFAOYSA-N', 'MXAYKZJJDUDWDS-UHFFFAOYSA-N']}. Unit IDs absent from the batch: none.

## Doses
- `Metadata_mmoles_per_liter` is in fact micromolar (10 = 10 uM top dose). Each well is mapped to the nearest nominal level in log space, [0.0412, 0.123, 0.37, 1.11, 3.33, 10.0] uM (the `units.csv` labels), if within 1.5x; otherwise the well is off-ladder (kept for normalization/feature selection/sphering, excluded from dose conditions and from the norm_z / other-compound references).
- Actual values per nominal dose: {0.0412: [0.028, 0.0334, 0.037, 0.0394, 0.0412, 0.0414, 0.0415, 0.0441, 0.0461, 0.0475, 0.0515, 0.0558, 0.0593], 0.123: [0.0839, 0.1001, 0.1111, 0.1183, 0.1235, 0.1241, 0.1244, 0.1323, 0.1383, 0.1426, 0.1545, 0.1673, 0.1778], 0.37: [0.2516, 0.3002, 0.3333, 0.3549, 0.3704, 0.3722, 0.373, 0.3969, 0.4148, 0.4277, 0.4633, 0.502, 0.5334], 1.11: [0.7549, 0.9007, 1.0, 1.0649, 1.1111, 1.1166, 1.1191, 1.1907, 1.2444, 1.283, 1.39, 1.5059, 1.6001], 3.33: [2.2646, 2.702, 3.1946, 3.3333, 3.3499, 3.3574, 3.5721, 3.7333, 3.8489, 4.1701, 4.5176, 4.8003], 10.0: [6.7937, 8.106, 9.5837, 10.0, 10.05, 10.072, 10.716, 11.2, 11.547]}.
- Off-ladder treated wells: 3383 (mostly the two 20 uM positive controls on every plate); bin collisions (same compound, two actual doses on one nominal level on one plate): 0. All matched units sit exactly on the standard ladder.

## Processing (one feature space for the whole batch)
- Feature columns: union 1781; 10 plates (['SQ00015116', 'SQ00015117', 'SQ00015118', 'SQ00015119', 'SQ00015120', 'SQ00015121', 'SQ00015122', 'SQ00015123', 'SQ00015125', 'SQ00015126']) lack 36 features (all DNA RadialDistribution), so the plate intersection is used: 1745 columns.
- Raw feature columns used: 1745. Dropped non-morphological bookkeeping columns (`*_Location_*`, `*_Number_Object_Number`, `*_Parent_*`, `*ObjectNumber*`): -76; NaN/inf > 1% of wells: -0; constant: -20 -> 1649. Imputed values (plate DMSO median): 6.
- Robust z per plate vs its 24 DMSO wells, 1.4826 MAD with floor 0.01 x median-over-plates scale; plate x feature scales at floor: 0.0552; clip +-50.
- Selection: pycytominer blocklist -55, outlier (pre-clip max|z| > 500) -139, variance threshold (freq_cut 0.05, unique_cut 0.01) -177, greedy |r| > 0.9 -568 -> **710 features**.
- Sphering: PCA whitening on 3264 DMSO wells, 99% variance -> **418 dims**, scale 1/sqrt(eigenvalue + 3.0) (eps fixed a priori to the JUMP T3 value, not tuned); eigenvalues first/last kept [235.87, 0.087], 387 kept eigenvalues below eps.
- Consensus: median of replicate wells per unit x nominal dose (single site `broad`; units with two IDs pool the wells of both).

## Read this first: well position dominates the DMSO-based calls
- All replicate plates of a platemap share one layout and every dose of a compound sits at a fixed well (100.0% of 11046 cross-plate replicate pairs of matched units share a well; the rest are the two-ID units). Each compound fills one row x 6-column block (10 -> 0.04 uM left to right); DMSO occupies only A01-06, B01-06, F07-12, J13-18.
- DMSO cosine (sphered), different plates: same position 0.0346 vs different position -0.0019 (same plate, different position 0.0479). Same position, by map and batch: same map, same batch 0.105 (n=3816); same map, other batch 0.027 (n=2592); other map, same batch 0.071 (n=39096); other map, other batch 0.025 (n=174816). The position effect is batch- and map-specific. Matched replicate pairs 0.2838; different compounds at the same dose 0.0229.
- Same-position DMSO sets (672, each DMSO position across the plates of one map) through the spec tests: `active_rep` fires in 56.1% (edge 64.8%, interior 52.5%), `active_norm_dmso` in 23.7% (edge 46.4%, interior 14.3%). **The spec calls `active_rep` and `active_norm_dmso` (and `active_vs_cpd`) are therefore not usable as 'departs from control' in this dataset**; they are kept because the spec mandates them.
- Position-matched DMSO null (added): `rep_pos_p` / `norm_pos_p` compare the condition with whole same-position DMSO sets of other platemaps (same edge/interior class, trimmed to the condition's batch composition: batch-1 plate + later-batch plates); null sizes [175, 459] sets (fewer than 1000; deviation). Leave-one-map-out on the DMSO sets gives 4.5% / 5.5% (edge 4.1% / 4.1%, interior 4.6% / 6.1%); near 5% by construction, so not the false-positive rate at compound wells.
- Operational bound: the library (every on-ladder compound x dose, 9383 conditions) at the lowest dose. Rate of `active_rep_pos` / `active_norm_pos` at 0.0412 uM: all wells 37.0% / 38.2%; reference region (rows C-N, columns 1-18, chosen from the per-well heatmap) 23.1% / 25.2%; reliable blocks 22.5% / 24.6%. Lowest-dose rate by region (rep_pos / norm_pos): cols 19-24 59.9% / 58.7%; outer rows (A,B,O,P) cols 1-18 48.0% / 50.0%; rows C-N cols 1-18 23.1% / 25.2%. No DMSO well samples columns 19-24 or rows M-P, so no DMSO null can calibrate them, and a treated-well null there would be dose-relative.
- Block reliability: a block (row x 6 columns) is unreliable if its lowest-dose well's library rate (over all 28 maps) of `active_rep_pos` or `active_norm_pos` exceeds the reference-region rate (one-sided binomial p < 0.01; 35 of 56 pass) AND lies in columns 1-18 (columns 19-24 have no DMSO well; rule read from the heatmap, where the binomial test lacks power with 28 compounds per well). 32 of 56 blocks reliable; unreliable: ['A_2', 'A_3', 'A_4', 'B_2', 'B_4', 'C_3', 'D_4', 'E_4', 'F_4', 'G_4', 'H_4', 'I_4', 'J_4', 'K_4', 'L_4', 'M_1', 'M_4', 'N_4', 'O_2', 'O_3', 'O_4', 'P_2', 'P_3', 'P_4']. **Matched units in reliable blocks: 104 of 187.** Wells flagged individually (Bonferroni, vs reference region at the same dose): rep_pos 31, norm_pos 36.
- Bound for the recommended column: fraction of library compounds in reliable blocks (all 6 doses tested) whose sustained active dose is 0.0412 uM: {'active_rep_pos': {'n_compounds': 891, 'frac_sustained_at_0.0412': 0.1493}, 'active_norm_pos': {'n_compounds': 891, 'frac_sustained_at_0.0412': 0.174}}. This includes real low-dose activity, so it is an upper bound on the positional false-positive rate, not an estimate.
- Persistence inside reliable blocks: {'P(act0.123|act0.0412)': 0.766, 'P(act0.123|not0.0412)': 0.105} (adjacent wells of one block, so not a clean artifact test).
- `active_vs_cpd_pos` (other maps' compounds at the same well and batch, distinct maps per slot) lacks the same-map component of the position effect (0.106 vs 0.071) and is anticonservative; 0 matched conditions have no feasible draw (NaN).
- Positive controls (bortezomib, MG-132 at 20 uM; one same-position set per map and well): [{'broad_sample': 'BRD-K50691590-001-02-2', 'pert_iname': 'bortezomib', 'n_sets': 336, 'median_rep_mAP': 1.0, 'frac_rep_p05': 1.0, 'frac_rep_pos_p05': 1.0, 'median_norm': 76.4769, 'frac_norm_dmso_p05': 1.0, 'frac_norm_pos_p05': 1.0}, {'broad_sample': 'BRD-K60230970-001-10-0', 'pert_iname': 'MG-132', 'n_sets': 336, 'median_rep_mAP': 1.0, 'frac_rep_p05': 1.0, 'frac_rep_pos_p05': 1.0, 'median_norm': 82.4503, 'frac_norm_dmso_p05': 1.0, 'frac_norm_pos_p05': 1.0}].

**Recommended for the lowest-dose comparison:** `sustained_active_dose_rep_pos` (or `_norm_pos`) restricted to `block_reliable` units, i.e. the lowest dose from which every higher tested dose is active under the position-matched DMSO null; each dose is a different well, so isolated positional false positives drop out. Treat `lowest_active_dose_rep_pos` as an optimistic bound.

## Calls (columns)
- `rep_mAP`: each replicate queries its other replicates + all DMSO wells on the condition's plates (cosine). `rep_p`: 2000 plate-matched random DMSO sets (one DMSO well per replicate plate; the drawn wells are removed from the negatives). `active_rep` = p < 0.05 (spec).
- `norm_z`: robust z of the consensus L2 norm vs all on-ladder compounds at that dose (n per dose {0.0412: 1566, 0.123: 1566, 0.37: 1566, 1.11: 1566, 3.33: 1561, 10.0: 1558}); `active_z` = z > 3 (spec).
- `norm_dmso_p` / `active_norm_dmso` (addendum a): consensus norm vs 2000 DMSO pseudo-consensus norms (median of one random DMSO well per replicate plate).
- `vs_cpd_p` / `active_vs_cpd` (addendum b): replicate mAP vs 2000 sets of other compounds' wells on the same plates at the same dose (one distinct compound per slot).
- `rep_pos_p`, `norm_pos_p`, `vs_cpd_pos_p` and their `active_*_pos` calls: position-matched variants (above). `pos_flag_*`, `block`, `block_reliable`, `block_lowdose_rate_*`: position quality.
- Calls are NaN where the test could not be run (n_wells < 2 or no feasible null).

### Per dose (matched units)
```
dose_uM                                0.0412    0.1230    0.3700    1.1100    3.3300    10.0000
n_units                                187.000   187.000   187.000   187.000   187.000   187.000
median_n_wells                           5.000     5.000     5.000     5.000     5.000     5.000
min_n_wells                              2.000     2.000     2.000     2.000     2.000     2.000
n_wells                                938.000   938.000   938.000   938.000   938.000   938.000
frac_reproducible                        0.824     0.775     0.824     0.850     0.866     0.957
frac_active_z                            0.182     0.198     0.230     0.278     0.289     0.305
frac_active_norm_dmso                    0.765     0.818     0.807     0.888     0.877     0.957
frac_active_vs_cpd                       0.604     0.610     0.663     0.684     0.684     0.754
frac_active_rep_pos                      0.465     0.513     0.561     0.636     0.642     0.824
frac_active_norm_pos                     0.497     0.556     0.594     0.647     0.684     0.861
frac_active_vs_cpd_pos                   0.449     0.471     0.540     0.578     0.668     0.754
median_rep_mAP                           0.216     0.233     0.274     0.371     0.441     0.899
n_units_reliable_block                 104.000   104.000   104.000   104.000   104.000   104.000
frac_active_rep_pos_reliable_blocks      0.356     0.356     0.471     0.635     0.635     0.837
frac_active_norm_pos_reliable_blocks     0.346     0.404     0.510     0.596     0.663     0.837
wells_all_cpds                        7606.000  7606.000  7606.000  7606.000  7584.000  7568.000
plates                                 136.000   136.000   136.000   136.000   136.000   136.000
unique_wells                           933.000   933.000   933.000   933.000   933.000   933.000
plates_matched                         136.000   136.000   136.000   136.000   136.000   136.000
```

### Library active fraction by dose (all on-ladder compounds; compound nulls not run for the library)
```
         active_rep  active_z  active_norm_dmso  active_vs_cpd  active_rep_pos  active_norm_pos  active_vs_cpd_pos
dose_uM                                                                                                           
0.0412        0.746     0.138             0.688            NaN           0.370            0.382                NaN
0.1230        0.746     0.141             0.699            NaN           0.404            0.448                NaN
0.3700        0.780     0.178             0.725            NaN           0.444            0.494                NaN
1.1100        0.813     0.204             0.767            NaN           0.496            0.527                NaN
3.3300        0.844     0.199             0.802            NaN           0.564            0.596                NaN
10.0000       0.931     0.223             0.906            NaN           0.739            0.759                NaN
```
Reliable blocks only:
```
         active_rep  active_z  active_norm_dmso  active_vs_cpd  active_rep_pos  active_norm_pos  active_vs_cpd_pos
dose_uM                                                                                                           
0.0412        0.603     0.068             0.507            NaN           0.225            0.246                NaN
0.1230        0.625     0.087             0.551            NaN           0.254            0.287                NaN
0.3700        0.702     0.120             0.606            NaN           0.332            0.369                NaN
1.1100        0.758     0.145             0.679            NaN           0.436            0.439                NaN
3.3300        0.809     0.174             0.748            NaN           0.533            0.547                NaN
10.0000       0.924     0.201             0.890            NaN           0.727            0.733                NaN
```
Units covered: 187 of 187; with all 6 doses tested: 187. Units with a non-monotonic call (active at a dose, inactive at a higher tested dose): {'active_rep': 48, 'active_z': 33, 'active_norm_dmso': 40, 'active_vs_cpd': 56, 'active_rep_pos': 43, 'active_norm_pos': 53, 'active_vs_cpd_pos': 55}.

## Outputs
- `step2/morph_lincs.parquet`: unit, dose_uM, cell_line, time_h, n_wells, n_plates, plates, sites, dose_actual_uM, f0..fN (float32); plus a `DMSO` reference row (dose 0, all DMSO wells).
- `step2/morph_lincs_library_activity.csv`: DMSO-based calls for every on-ladder compound x dose (13-char prefix), used for the norm_z reference and position diagnostics.
- `step2/morph_lincs_activity.csv`, `step2/morph_lincs_lowest_active.csv`, `step2/morph_lincs_wells.parquet` (8892 rows, 21.6 MB; matched-unit wells + all DMSO wells; sphered).

## Resources
- Bytes downloaded: 280,564,937 (profiles 279,328,291); peak disk 26,780,346 bytes (metadata cache + outputs; profiles are held in memory only); download 42.8 s, total 368.8 s.
- Python 3.13.0; pandas 2.3.3, numpy 2.2.6, scipy 1.16.3, scikit-learn 1.8.0, pyarrow 23.0.1, requests 2.32.5, fsspec 2026.2.0

## Deviations from MORPH_SPEC
- Non-morphological bookkeeping columns (Location, object numbers, parent IDs) dropped before step 1; they encode image coordinates, not morphology.
- Features absent on some plates (36 DNA RadialDistribution columns missing on 10 plates) are dropped for the whole batch (plate intersection) to keep one feature space.
- pycytominer blocklist applied in feature selection, as in the JUMP T3 run.
- Off-ladder wells (beyond 1.5x of a nominal level, chiefly the 20 uM positive controls) are excluded from dose conditions instead of being forced to the nearest level (20 uM would otherwise land on 10 uM).
- dose_uM uses the `units.csv` labels (0.0412, 0.123, 0.37, 1.11, 3.33, 10) rather than rounded values.
- In the DMSO null of `rep_mAP` the drawn DMSO wells are removed from the negatives, so null pools are smaller by n than the observed pools (slightly conservative).
- Added `active_norm_dmso` and `active_vs_cpd` per the coordinator addendum, plus position-matched variants, block reliability and sustained lowest doses (position effect, see top). The position-matched DMSO null has 175-459 sets, fewer than the 1000 draws the spec asks for (only 24 DMSO positions x 27 other maps exist).
- Wells absent from the profile files: ['SQ00015207 P01'] (52,223 of 52,224 wells).

## Log
```
accession: cpg0004-lincs; batch 2016_04_01_a549_48hr_batch1
plates: 136 profiled of 140 in barcode_platemap; missing {'SQ00015225': 'C-7161-01-LM6-012', 'SQ00015226': 'C-7161-01-LM6-013', 'SQ00015227': 'C-7161-01-LM6-013', 'SQ00015228': 'C-7161-01-LM6-013'}
feature columns: union 1781, intersection 1745; 10 plates lack [36] features
wells: 52223, raw feature columns 1745, DMSO wells 3264
off-ladder treated wells: 3383; bin collisions 0
features: raw 1745, bookkeeping -76, nan>1% -0, const -20 -> 1649; blocklist -55, outlier -139, lowvar -177, corr -568 -> 710
sphering: 418 dims, eig first/last [235.87, 0.087], eps 3.0
matched conditions: 1122 (187 units)
matched-condition tests: 267s
library tests (9383 compound x dose): 10s
library active fraction by dose:
         active_rep  active_z  active_norm_dmso  active_vs_cpd  active_rep_pos  active_norm_pos  active_vs_cpd_pos
dose_uM                                                                                                           
0.0412        0.746     0.138             0.688            NaN           0.370            0.382                NaN
0.1230        0.746     0.141             0.699            NaN           0.404            0.448                NaN
0.3700        0.780     0.178             0.725            NaN           0.444            0.494                NaN
1.1100        0.813     0.204             0.767            NaN           0.496            0.527                NaN
3.3300        0.844     0.199             0.802            NaN           0.564            0.596                NaN
10.0000       0.931     0.223             0.906            NaN           0.739            0.759                NaN
flagged wells (vs reference region, Bonferroni): rep_pos 31, norm_pos 36
blocks: 56, reliable 32; unreliable ['A_2', 'A_3', 'A_4', 'B_2', 'B_4', 'C_3', 'D_4', 'E_4', 'F_4', 'G_4', 'H_4', 'I_4', 'J_4', 'K_4', 'L_4', 'M_1', 'M_4', 'N_4', 'O_2', 'O_3', 'O_4', 'P_2', 'P_3', 'P_4']
library active fraction by dose, reliable blocks only:
         active_rep  active_z  active_norm_dmso  active_vs_cpd  active_rep_pos  active_norm_pos  active_vs_cpd_pos
dose_uM                                                                                                           
0.0412        0.603     0.068             0.507            NaN           0.225            0.246                NaN
0.1230        0.625     0.087             0.551            NaN           0.254            0.287                NaN
0.3700        0.702     0.120             0.606            NaN           0.332            0.369                NaN
1.1100        0.758     0.145             0.679            NaN           0.436            0.439                NaN
3.3300        0.809     0.174             0.747            NaN           0.533            0.547                NaN
10.0000       0.924     0.201             0.890            NaN           0.727            0.733                NaN
library (reliable blocks) sustained-from-0.0412 fraction: {'active_rep_pos': {'n_compounds': 891, 'frac_sustained_at_0.0412': 0.1493}, 'active_norm_pos': {'n_compounds': 891, 'frac_sustained_at_0.0412': 0.174}}
positive controls: [{'broad_sample': 'BRD-K50691590-001-02-2', 'pert_iname': 'bortezomib', 'n_sets': 336, 'median_rep_mAP': 1.0, 'frac_rep_p05': 1.0, 'frac_rep_pos_p05': 1.0, 'median_norm': 76.4769, 'frac_norm_dmso_p05': 1.0, 'frac_norm_pos_p05': 1.0}, {'broad_sample': 'BRD-K60230970-001-10-0', 'pert_iname': 'MG-132', 'n_sets': 336, 'median_rep_mAP': 1.0, 'frac_rep_p05': 1.0, 'frac_rep_pos_p05': 1.0, 'median_norm': 82.4503, 'frac_norm_dmso_p05': 1.0, 'frac_norm_pos_p05': 1.0}]
same-position DMSO calibration: {'n_sets': 672, 'n_edge_sets': 196, 'pos_null_sizes': [175, 459], 'rep_p<0.05': 0.561, 'rep_p<0.05_edge': 0.648, 'rep_p<0.05_interior': 0.5252, 'norm_dmso_p<0.05': 0.2366, 'norm_dmso_p<0.05_edge': 0.4643, 'norm_dmso_p<0.05_interior': 0.1429, 'rep_pos_p<0.05': 0.0446, 'rep_pos_p<0.05_edge': 0.0408, 'rep_pos_p<0.05_interior': 0.0462, 'norm_pos_p<0.05': 0.0551, 'norm_pos_p<0.05_edge': 0.0408, 'norm_pos_p<0.05_interior': 0.0609}
DMSO cos same-pos/diff-plate 0.0346, diff-pos/diff-plate -0.0019; replicate cos 0.2838 (same well 1.000); other-cpd same-dose cos 0.0229
         n_units  median_n_wells  min_n_wells  n_wells  frac_reproducible  frac_active_z  frac_active_norm_dmso  frac_active_vs_cpd  frac_active_rep_pos  frac_active_norm_pos  frac_active_vs_cpd_pos  median_rep_mAP  n_units_reliable_block  frac_active_rep_pos_reliable_blocks  frac_active_norm_pos_reliable_blocks  wells_all_cpds  plates  unique_wells  plates_matched
dose_uM                                                                                                                                                                                                                                                                                                                                                                        
0.0412       187             5.0            2      938              0.824          0.182                  0.765               0.604                0.465                 0.497                   0.449           0.216                     104                                0.356                                 0.346            7606     136           933             136
0.1230       187             5.0            2      938              0.775          0.198                  0.818               0.610                0.513                 0.556                   0.471           0.233                     104                                0.356                                 0.404            7606     136           933             136
0.3700       187             5.0            2      938              0.824          0.230                  0.807               0.663                0.561                 0.594                   0.540           0.274                     104                                0.471                                 0.510            7606     136           933             136
1.1100       187             5.0            2      938              0.850          0.278                  0.888               0.684                0.636                 0.647                   0.578           0.371                     104                                0.635                                 0.596            7606     136           933             136
3.3300       187             5.0            2      938              0.866          0.289                  0.877               0.684                0.642                 0.684                   0.668           0.441                     104                                0.635                                 0.663            7584     136           933             136
10.0000      187             5.0            2      938              0.957          0.305                  0.957               0.754                0.824                 0.861                   0.754           0.899                     104                                0.837                                 0.837            7568     136           933             136
```
