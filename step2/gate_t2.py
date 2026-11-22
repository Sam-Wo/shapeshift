"""T2 gate: are there enough HepG2/C3A (CVCL_1098) cells per drug x dose for a stable pseudobulk?
Reads the cached DE combination index (step2/tahoe_cache/de_index_combos.parquet; built from parquet footers only).
Threshold: >= 100 treated cells per condition, and a median >= 100 across T2 drugs at 5 uM.
Rationale: with ~100 cells the pseudobulk of a gene with baseMean ~40 (our gene-subset floor) pools
thousands of UMIs, so the lfc of well-expressed genes is stable; Tahoe-wide medians are ~1,200 cells,
and the four weakest lines (13-80 cells) are the cases this threshold is meant to catch."""
import glob, json, os
import pandas as pd

THRESH = 100
# final combination index written by step2/tahoe_signatures.py (--stage index); one row per drug x dose x line x plate
combo = pd.read_parquet("step2/tahoe_cache/de_index_combos.parquet")
combo = combo.assign(drug=combo.drug_stripped, cellosaur_min=combo.cellosaur, conc_min=combo.dose_uM.astype(float),
                     n_trt_min=combo.n_trt, n_ctrl_min=combo.n_ctrl)
u = pd.read_csv("step2/units.csv")
t2 = {d for s in u[u.tier == "T2_EUOS_HepG2"].tahoe_drugs for d in s.split(";")}
h = combo[(combo.cellosaur_min == "CVCL_1098") & combo.drug.isin(t2)]
out = {}
for conc, g in h.groupby("conc_min"):
    out[f"{conc:g}uM"] = dict(n_conditions=len(g), n_drugs=g.drug.nunique(), median_cells=float(g.n_trt_min.median()),
                              p10_cells=float(g.n_trt_min.quantile(.1)), min_cells=int(g.n_trt_min.min()),
                              frac_ge_threshold=float((g.n_trt_min >= THRESH).mean()),
                              median_ctrl_cells=float(g.n_ctrl_min.median()))
lines = combo.groupby("cellosaur_min").n_trt_min.median().sort_values()
out["all_lines_median_cells_lowest5"] = lines.head(5).to_dict()
out["all_lines_median_of_medians"] = float(lines.median())
out["threshold"] = THRESH
out["T2_drugs_total"] = len(t2)
out["T2_drugs_missing_at_5uM"] = sorted(t2 - set(h[h.conc_min.round(3) == 5].drug))
out["verdict"] = "PASS" if out.get("5uM", {}).get("median_cells", 0) >= THRESH else "UNDERPOWERED"
json.dump(out, open("step2/gate_t2.json", "w"), indent=1)
print(json.dumps(out, indent=1))
