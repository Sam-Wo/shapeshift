# Shared morphology-processing spec (T1 LINCS, T2 EU-OS, CPJUMP1 calibration; T3 JUMP already follows this)

Environment: Windows, `py -3.13 -I -X utf8`. Use only installed packages: pandas, numpy, scipy, scikit-learn, pyarrow, fsspec, requests, rdkit and matplotlib. Write pinned versions to your notes. Do not pip install.

Storage: stream per-plate profile files, keep only the wells you need, and delete raw downloads when you're done. Keep at most 1.5 GB on disk at any time. Count every downloaded byte; `step2/gallery_index.py` provides `Bytes`, or keep your own counter. Record peak disk use. Fail loudly with a raised error if a source is unreachable; never skip silently.

Gallery accessions: resolve them with `step2/gallery_index.py` (`resolve('lincs'|'cpjump1'|'euos'|'jump')`, `s3_list`, `s3_url`). Never hardcode `cpgNNNN`.

Compound key: `step2/units.csv` (columns `tier`, `unit`, `tahoe_drugs`, `dataset_ids`, `dataset_names`, `moa`, `stereo_flag`, `exact_achiral`, ...). `dataset_ids` holds that tier's dataset compound IDs, semicolon-joined; one unit may map to several IDs (salts, batches). The per-drug provenance is in `cellcheck/cellline_matches.csv` (`dataset_cpd_ids`, doses, wells). Map dataset wells to units through `dataset_ids`. Process ALL treated compounds of the dataset, not only matched ones: sphering and feature selection need the full plate context. Output only matched units plus the DMSO/negative-control reference.

Processing, strictly within your dataset (CellProfiler versions differ between datasets, so never merge features across datasets):
1. Use precomputed well-level CellProfiler profiles, never images. Drop `Metadata_*` columns. Drop features that are NaN or inf in more than 1% of wells, or constant. Impute the remaining NaNs with the plate's negative-control median.
2. Normalize per plate: robust z against the plate's negative-control wells (DMSO or the dataset's vehicle), `(x − median_ctrl) / (1.4826 · MAD_ctrl)`, with a MAD floor. Clip to ±50.
3. Select features per dataset: drop near-zero-variance features, drop extreme-outlier features, and drop one of each pair with |r| > 0.9 (greedy). Report the count after each step.
4. Sphere per dataset: PCA-whitening fitted on negative-control wells only, regularized, keeping 99% of the variance. Report the dimensionality.
5. Build the consensus per unit × dose (× cell line × time point, where relevant): median across replicate wells. If there are several sites or sources, take the median within each site, then the mean across sites. Record `n_wells`, the plates and the sites.
6. Replicate reproducibility, the metric that makes a weak tier visible: for each unit × dose, compute replicate mAP. Each replicate well queries a pool of its other replicates plus all negative-control wells on the same plates, ranked by cosine similarity. Build a null by repeating this with random sets of negative-control wells of the same size (≥ 1000 draws). Call a condition reproducible if p < 0.05. Report the fraction reproducible per dose.
7. Activity, i.e. "departs from control", per unit × dose. Use the replicate-mAP p < 0.05 from step 6 as `active_rep`. As a second call, `active_z`, compute a robust z of the consensus L2 norm against all treated compounds at that dose: z = (norm − median) / (1.4826 · MAD), active if z > 3. Report both.

Outputs, where `<tier>` is `lincs`, `euos` or `cpjump1`:
- `step2/morph_<tier>.parquet`: unit, dose_uM, [cell_line, time_h], n_wells, n_plates, plates, sites, then feature columns f0..fN (float32).
- `step2/morph_<tier>_activity.csv`: unit, dose_uM, [cell_line, time_h], n_wells, rep_mAP, rep_p, active_rep, norm_z, active_z.
- `step2/morph_<tier>_wells.parquet`: sphered well-level profiles plus metadata, under ~300 MB.
- `step2/<tier>_notes.md`: paths resolved, file formats, feature counts, sphering dimensions, counts per dose, the reproducibility fraction per dose, bytes downloaded, peak disk, timing, pinned package versions, and any deviation from this spec with the reason.
- `step2/<tier>_download_summary.json`: `{"bytes": ..., "peak_disk_bytes": ..., "seconds": ...}`.

Do not modify other agents' files: `step2/tahoe_*`, `step2/jump_*`, `morph_profiles*.parquet`, `step2/compare*`, `units.csv`, `overlap.*`, `REPORT*.md`, `CELL_LINE_CHECK.md`.
