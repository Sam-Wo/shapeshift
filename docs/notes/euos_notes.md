# EU-OS T2 morphology profiles (HepG2, 24 h, 10 uM): notes (written by step2/morph_euos.py)

## Sources
- Gallery accession (resolved): `cpg0036-EU-OS-bioactives`. Per-site workspace prefixes: FMP: load_data_csv, load_data_csv_orig, metadata; IMTM: load_data_csv, load_data_csv_orig, metadata; MEDINA: load_data_csv, load_data_csv_orig, metadata; USC: load_data_csv, load_data_csv_orig, metadata. There is **no `workspace/profiles` or `backend` prefix** (only images, load_data_csv, metadata; plus `all/workspace/metadata_harmonized/` with one 0.88 GB draft metadata CSV, not used).
- Profiles: Zenodo concept 13309565 -> latest record **19347244 (v3.0.0, 2026-03-31)**, DOI 10.5281/zenodo.19347244. File `Aggregated_Profiles.zip` (1,348,543,344 bytes): https://zenodo.org/api/records/19347244/files/Aggregated_Profiles.zip/content
- Versions of the record (REST API): 14776021 vNone (2024-01-31): Aggregated_Profiles.zip 1,339,488,920 B; 13309566 v1.0.0 (2024-08-13): Aggregated_Profiles.zip 1,339,488,920 B; 19347244 v3.0.0 (2026-03-31): Aggregated_Profiles.zip 1,348,543,344 B. Only v3.0.0 carries the IMTM plate-swap fix.
- The zip is never downloaded whole: its central directory is read with two HTTP range requests and each HepG2 plate CSV member is fetched with one range request, inflated in memory and stored as float32 parquet in `step2/euos_cache/plates/` (deleted at the end of the run).
- Zip members used: {'FMP': 28, 'IMTM': 29, 'MEDINA': 28, 'USC': 28} HepG2 CSVs (the 28 `FMP_U2OS/` members are skipped). Platemaps: gallery `<acc>/<site>/workspace/metadata/platemaps/*HepG2/platemap/*.txt` (cached in `step2/euos_cache/metadata/platemaps/`).

## File format and quirks
- One CSV per plate, one row per well (CellProfiler features aggregated per well by the authors): `Metadata_Batch, Metadata_Plate, Metadata_Well, Metadata_Object_Count` + 2977 features (prefixes `Nuc_`, `Cells_`, `Cyto_`). The feature set is **identical at all four sites** (same 2981-column header), so no intersection was needed.
- FMP files store the plate as `Metadata_Plate='B1001'` with the replicate in `Metadata_Batch='R1'`; the other sites store `B1001_R1`. The script builds the canonical `B1001_R1` and checks it against the file name.
- IMTM plate swap (Zenodo README + `2026-02-19_SwappedID_IMTM.pdf`): IMTM has 29 files and `B1001_R4` appears twice. Resolved by object-count pattern (per-well log count vs. the leave-one-out mean of the other replicates of each layout): {"IMTM:B1001_R4": [{"name": "Aggregated_Profiles/aggregated_data/IMTM_HepG2/2023-07-17_HepG2_10uM_B1001_R4_CP_Profiles_Aggregated.csv", "r_own": 0.927, "best_other": "B1007", "r_other": 0.023}, {"name": "Aggregated_Profiles/aggregated_data/IMTM_HepG2/2023-07-07_HepG2_10uM_B1001_R4_CP_Profiles_Aggregated.csv", "r_own": 0.036, "best_other": "B1007", "r_other": 0.983}]}. Dropped: ['Aggregated_Profiles/aggregated_data/IMTM_HepG2/2023-07-07_HepG2_10uM_B1001_R4_CP_Profiles_Aggregated.csv'] (the stale pre-fix file named in the README). Label check on all 112 kept plates: per-well log-count correlation with own layout median 0.893 vs best other layout 0.082; smallest margin 0.285; no suspects.
- Wells: 43008 in platemaps, 42847 in files, 161 platemap wells absent from the files (cause not recorded; plausibly wells without usable cells or images): FMP 55, IMTM 84, MEDINA 17, USC 5; of these, matched-unit wells: 17. Units with fewer than 16 wells: GUGOEEXESWIERI-UHFFFAOYSA-N (EOS101413, 13 wells); KKGQTZUTZRNORY-UHFFFAOYSA-N (EOS101578, 13 wells); PHXJVRSECIGDHY-UHFFFAOYSA-N (EOS100606, 14 wells); XDXDZDZNSLXDNA-UHFFFAOYSA-N (EOS101383, 7 wells).
- Layout: DMSO = 28 wells per plate, all in the edge columns 23/24; compounds in columns 1-22; nocodazole/tetrandrine 5 uM in columns 23/24 (positions differ by site). Within a site the layout of B100x_R1..R4 is identical; across sites only the control positions differ. Negative control = `Metadata_EOS == 'DMSO'` (never concentration 0).
- The 'concentration 0' control wells (14 nocodazole + 14 tetrandrine at each of MEDINA and USC) behave like the 5 uM controls (mean cosine 0.71-0.75 to the 5 uM consensus of the same control and site, similar norms; table below), i.e. a metadata quirk. They are reported, never used.

## Processing (MORPH_SPEC steps 1-7)
- Features: raw 2977 -> bookkeeping -24 (BoundingBox min/max, Center_X/Y, closest-object numbers; deviation, see below) -> NaN/inf >1% -0, constant -70 = 2883 -> blocklist -36 (pycytominer list translated to `Nuc_`), outliers (pre-clip max|z| > 500) -204, low variance (freq_cut 0.05, unique_cut 0.01) -342, correlated (greedy |r| > 0.9, r estimated on 20000 random wells) -1585 = **716 selected**.
- Imputed values: 0; fraction of plate x feature MADs at the floor: 0.0298.
- Normalization: robust z per plate (site x plate) against its 28 DMSO wells, MAD floor = 0.01 x median plate scale, clip +-50.
- Sphering: PCA-whitening fit on all 3136 DMSO wells pooled over the four sites, 99% variance, eps 3.0 (same as T3): **530 dimensions** (eigenvalues 153.13 ... 0.0999).
- **Site handling:** one feature selection and one sphering for all four sites (the feature sets are identical, so nothing is lost to an intersection). Consensus = median of wells within each site, then mean across sites, in that single sphered space. Why: a cross-site mean is only meaningful if the sites share one space; per-plate robust z against each plate's own DMSO reduces (does not remove) site and plate offsets before the joint selection and sphering: replicate-pair cosine is 0.3676 within a site vs 0.1449 across sites. Per-site reproducibility is reported too (within-site mAP in the shared space). Separate per-site sphering would make the cross-site mean meaningless, so it was not used.
- Per site: FMP: 28 plates, 10697 wells (784 DMSO), 2977 raw features, sphered DMSO mean variance 0.21; IMTM: 28 plates, 10668 wells (784 DMSO), 2977 raw features, sphered DMSO mean variance 0.175; MEDINA: 28 plates, 10735 wells (784 DMSO), 2977 raw features, sphered DMSO mean variance 0.185; USC: 28 plates, 10747 wells (784 DMSO), 2977 raw features, sphered DMSO mean variance 0.177
- Replicate mAP: each replicate well queries its other replicates + all DMSO wells of the replicate plates (cosine, sphered). Null: the replicates are replaced by DMSO wells drawn from the same plates (one per replicate well), 1000 draws, cached per plate set; reproducible if p < 0.05. Overall = all wells of the unit (all sites); per site = that site's 4 wells.
- active_z: robust z of the consensus L2 norm vs all 2464 treated compounds (median 5.9156, 1.4826 MAD 2.851), active if z > 3.
- Addendum calls: `active_norm_dmso` (consensus norm vs DMSO pseudo-consensus: one DMSO well per replicate well from the same plate, median within site, mean across sites, 1000 draws); `active_vs_cpd` (mean pairwise cosine of all replicate wells vs the same score for other compounds' wells drawn from the same plates, distinct compounds within a draw); `active_vs_pos` (extra, position-matched: mean within-site replicate-pair cosine vs wells at the SAME well position, site and replicate index on plates of other layouts, i.e. different compounds at the same position, distinct layouts within a site); `active_cpd_and_pos` = active_vs_cpd AND active_vs_pos. What each null matches: `active_vs_cpd` matches plates (hence site and batch) but not well position; `active_vs_pos` matches site, replicate index and well position but not batch where a layout's replicates were imaged together (see batch sharing below).

## Counts (dose 10 uM, HepG2, 24 h)
- Units: **172 / 172** covered (EOS ids 173; one unit has two EOS ids). Wells per unit: {7: 1, 13: 2, 14: 1, 16: 167, 32: 1} (16 = 4 sites x 4 replicates; fewer = wells absent from the files; 32 = two EOS ids).
- All wells: {"other_compound": 36512, "negative_control": 3136, "matched_compound": 2751, "positive_control": 448}; all 2464 treated compounds are used for normalization context, selection and the active_z reference.

## Reproducibility and activity (172 units)
- **Headline: site effects dominate the replicate structure.** Replicate-pair cosine 0.3676 within a site vs 0.1449 across sites; median replicate mAP 0.4374 over all 16 wells vs FMP 0.9583, IMTM 1.0, MEDINA 0.7978, USC 0.6348 within site. Same unit, consensus at two different sites: cosine 0.2307 (by site pair {'FMP-IMTM': 0.2233, 'FMP-MEDINA': 0.2799, 'FMP-USC': 0.2257, 'IMTM-MEDINA': 0.1914, 'IMTM-USC': 0.1977, 'MEDINA-USC': 0.2661}) vs 0.0402 for different units at different sites, so every site pair carries shared compound signal.
- Replicate mAP reproducible (p<0.05), overall: 100.0% (DMSO-null mean mAP 0.0426); saturated, see caveat.
- Per site (within-site mAP): FMP 95.3% (median mAP 0.9583), IMTM 94.2% (median mAP 1.0), MEDINA 94.2% (median mAP 0.7978), USC 98.3% (median mAP 0.6348). Number of sites reproducible per unit: {1: 2, 2: 5, 3: 15, 4: 150}.
- active_z 27.3%; active_norm_dmso 100.0%; active_vs_cpd 81.4%; active_vs_pos 95.3% (per site, within-site score: FMP 76.2%, IMTM 86.1%, MEDINA 75.6%, USC 73.3%); active_cpd_and_pos 80.2%.
- Batch sharing of within-site replicate pairs (platemap batch folders): {'FMP': 0.0, 'IMTM': 0.855, 'MEDINA': 0.238, 'USC': 0.0}. At FMP batch = replicate index, so `active_vs_pos` is batch-matched there; at IMTM a layout's 4 replicates were imaged in one batch, so replicate pairs share a batch and null pairs do not. Compare the FMP and IMTM per-site rates above to judge the batch contribution.
- **Saturation caveat:** `active_rep` and `active_norm_dmso` are saturated and should not be used as activity calls here. DMSO sits only in edge columns 23/24 while compounds sit in columns 1-22, and all replicates of a compound share one well position at every site (replicate pairs at the same well: 100%), so any compound set differs from DMSO by an edge/position offset. DMSO wells at the same position on different plates: cosine 0.0355 vs 0.0003 at different positions; different compounds at the same position: 0.124 vs 0.0865 at different positions. Neither addendum null controls everything (position vs batch, see above); `active_z` is the most conservative call.

## Positive controls (5 uM; excluded from the compound set)
```
    control  dose_uM   site  n_wells  n_plates  rep_mAP    rep_p  null_mAP_mean   norm  norm_z_vs_cpds  n_wells_conc0  conc0_cos_to_5uM_consensus  conc0_norm
 Nocodazole      5.0    FMP       56        28   1.0000 0.000999         0.0776 30.882            8.76              0                         NaN         NaN
 Nocodazole      5.0   IMTM       56        28   1.0000 0.000999         0.0770 35.913           10.52              0                         NaN         NaN
 Nocodazole      5.0 MEDINA       42        21   0.9658 0.000999         0.0785 31.044            8.81             14                       0.719      31.504
 Nocodazole      5.0    USC       42        21   1.0000 0.000999         0.0776 23.613            6.21             14                       0.725      21.174
 Nocodazole      5.0    all      196        98   0.9867      NaN            NaN 25.730            6.95             28                       0.636      24.135
Tetrandrine      5.0    FMP       56        28   1.0000 0.000999         0.0774 43.193           13.08              0                         NaN         NaN
Tetrandrine      5.0   IMTM       56        28   0.9996 0.000999         0.0769 33.097            9.53              0                         NaN         NaN
Tetrandrine      5.0 MEDINA       42        21   0.9992 0.000999         0.0783 25.790            6.97             14                       0.711      39.758
Tetrandrine      5.0    USC       42        21   0.9983 0.000999         0.0778 20.933            5.27             14                       0.744      28.790
Tetrandrine      5.0    all      196        98   0.9983      NaN            NaN 26.354            7.17             28                       0.631      32.136
```

## Outputs
- `step2/morph_euos.parquet`: 172 units + a `DMSO` reference row (dose_uM 0), f0..f529 float32.
- `step2/morph_euos_activity.csv`: spec columns + norm, norm_dmso_p, active_norm_dmso, rep_cos, vs_cpd_p, active_vs_cpd, rep_cos_within_site, vs_pos_p, active_vs_pos, active_cpd_and_pos, vs_pos_p_<site>, per-site rep_mAP/rep_p, n_sites_reproducible, eos_ids, wells.
- `step2/morph_euos_controls.csv` (positive controls); `step2/morph_euos_wells.parquet` (sphered wells of matched units and DMSO, with site/batch/plate/layout/rep/well/eos/unit/role/conc/object_count); `step2/euos_cache/euos_run.json` (full run info incl. selected feature names); `step2/euos_cache/metadata/` (Zenodo record JSON, zip directory, platemaps).

## Bytes, disk, time, versions
- This run: 1,083,485,997 bytes over HTTP (plate members 1,081,580,277); peak disk 885,926,580 bytes (filtered float32 plate cache 884,190,063 bytes, deleted: True); raw zip/CSV bytes never touch the disk. Wall time 273.5 s (download 100.6 s).
- Cumulative over development: two full plate fetches (about 2 x 1.08 GB) plus about 13 MB of manual probes (zip directory listings, one IMTM CSV, the swap PDF, header range reads): about 2.2 GB in total.
- Python 3.13; packages: {'pandas': '2.3.3', 'numpy': '2.2.6', 'pyarrow': '23.0.1', 'requests': '2.32.5'}.

## Deviations from MORPH_SPEC
- Profiles come from Zenodo (not the gallery), fetched by range requests from one zip, because the gallery holds no profiles.
- 24 bookkeeping columns (object bounding-box coordinates, centers, closest-object numbers) are dropped before the NaN/constant step: in these raw per-well aggregates they are not morphology.
- Feature correlation for the |r| > 0.9 pruning is estimated on a 20,000-well random subsample (memory); all other steps use all wells.
- Blocklist names translated from `Nuclei_` to `Nuc_` (this dataset's compartment name).
- Extra activity calls `active_vs_pos` (position-matched) and `active_cpd_and_pos`, because every replicate shares a well position and DMSO is confined to the edge columns, which saturates the DMSO-based calls.
- The two IMTM `B1001_R4` files are resolved from the data (object-count pattern); the result agrees with the Zenodo README.
