"""ShapeShift cell-line check: Tahoe-100M drugs vs. Cell Painting datasets in Tahoe cell lines
(A549 = LINCS cpg0004 + CPJUMP1 cpg0000; HepG2 = EU-OPENSCREEN cpg0036).

Usage (from the project root):  py -3.13 -I -X utf8 cellcheck/cellline_match.py
Inputs: small metadata files in cellcheck/data/ (fetched by hand, see CELL_LINE_CHECK.md), overlap.csv,
step2/compounds.csv, data/tahoe/drug_metadata.parquet.  Nothing in data/ or step2/ is modified.
Outputs: cellcheck/cellline_matches.csv, cellcheck/data/dataset_compounds_*.csv, cellcheck/data/summary.json.

Matching uses overlap.py's own std_mol / std_key / graph_hash / norm_name / METALS so the handling is identical:
  tier 1  : std InChIKey equal (Tahoe SMILES or PubChem isomeric SMILES vs dataset SMILES and dataset's own key);
            1a if PubChem SMILES has no stereo, 1b otherwise
  tier 2  : InChIKey block 1 equal
  tier 2b : same molecular formula + same heavy-atom graph (tautomer / bond-shift drawings)
  tier 3  : metal-containing drug (demoted), or name bridge (norm_name) for drugs without a structure match
"""
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from overlap import METALS, graph_hash, norm_name, std_key, std_mol  # noqa: E402

D = ROOT / "cellcheck" / "data"
TAHOE_DOSES = (0.05, 0.5, 5.0)


def fmt_doses(ds):
    return ";".join(f"{d:g}" for d in sorted(ds))


# ---------------------------------------------------------------- dataset compound tables
def lincs_structures():
    """broad_sample -> (smiles, inchikey, name, how) using all Repurposing Hub sample releases."""
    rel = {}
    for d in ["20200324", "20180907", "20180516", "20170327"]:  # newest first
        r = pd.read_csv(D / f"lincs/repurposing_samples_{d}.txt", sep="\t", comment="!", encoding="latin-1")
        r["release"] = d
        rel[d] = r
    allr = pd.concat(rel.values(), ignore_index=True)
    exact = allr.drop_duplicates("broad_id").set_index("broad_id")
    dep = {}
    for r in allr.dropna(subset=["deprecated_broad_id"]).itertuples():
        for x in str(r.deprecated_broad_id).split("|"):
            dep.setdefault(x.strip(), r)
    allr["brd13"] = allr["broad_id"].str[:13]
    by13 = allr.groupby("brd13")

    def resolve(bs):
        if bs in exact.index:
            r = exact.loc[bs]
            return r.smiles, r.InChIKey, r.pert_iname, f"hub_samples_{r.release}"
        if bs in dep:
            r = dep[bs]
            return r.smiles, r.InChIKey, r.pert_iname, f"hub_deprecated_id_{r.release}"
        if bs[:13] in by13.groups:
            g = by13.get_group(bs[:13])
            keys = {std_key(s) for s in g.smiles.dropna()} - {None}
            if len(keys) == 1:
                r = g.iloc[0]
                return r.smiles, r.InChIKey, r.pert_iname, "hub_brd13_prefix"
            return None, None, g.pert_iname.iloc[0], f"hub_brd13_prefix_ambiguous({len(keys)} keys)"
        return None, None, None, None
    return resolve


def load_lincs():
    pm = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(str(D / "lincs/platemap/*.txt"))])
    pm = pm.dropna(subset=["broad_sample"])
    ext = pd.read_csv(D / "lincs/external_metadata.tsv", sep="\t").set_index("broad_sample")
    resolve = lincs_structures()
    rows = []
    for bs, g in pm.groupby("broad_sample"):
        smi, ik, nm, how = resolve(bs)
        e = ext.loc[bs] if bs in ext.index else None
        rows.append(dict(
            dataset="LINCS_cpg0004_A549_48h", cpd_id=bs, smiles=smi, inchikey=ik,
            inchikey14=(e["InChIKey14"] if e is not None and isinstance(e["InChIKey14"], str) else None),
            name=nm if isinstance(nm, str) else (e["pert_iname"] if e is not None else None),
            structure_source=how or ("external_metadata_InChIKey14_only" if e is not None and isinstance(e["InChIKey14"], str) else None),
            doses_uM=fmt_doses({float(f"{x:.3g}") for x in g.mmoles_per_liter.dropna()}),
            n_wells=len(g), cell_lines="A549", time_h="48"))
    return pd.DataFrame(rows)


def load_lincs_batch2():
    """2017_12_05_Batch2 pilot (A549/MCF7/U2OS, 6/24/48 h). A549 platemaps only."""
    pm = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(str(D / "lincs/batch2_platemap/*_A549_*H.txt"))])
    pm = pm.dropna(subset=["broad_sample"])
    pm = pm[pm.broad_sample != "NA"]
    ext = pd.read_csv(D / "lincs/external_metadata.tsv", sep="\t")
    ext13 = ext.assign(b13=ext.broad_sample.str[:13]).drop_duplicates("b13").set_index("b13")
    resolve = lincs_structures()
    rows = []
    for bs, g in pm.groupby("broad_sample"):
        smi, ik, nm, how = resolve(bs)
        e = ext13.loc[bs[:13]] if bs[:13] in ext13.index else None
        rows.append(dict(
            dataset="LINCS_cpg0004_Batch2pilot_A549", cpd_id=bs, smiles=smi, inchikey=ik,
            inchikey14=(e["InChIKey14"] if e is not None and isinstance(e["InChIKey14"], str) else None),
            name=nm if isinstance(nm, str) else (e["pert_iname"] if e is not None else None),
            structure_source=how, doses_uM=fmt_doses({float(f"{x:.3g}") for x in g.mmoles_per_liter.dropna()}),
            n_wells=len(g), cell_lines="A549 (also MCF7, U2OS)",
            time_h=";".join(sorted({t.rstrip("H") for t in g.time_point}, key=float))))
    return pd.DataFrame(rows)


def load_cpjump1():
    meta = pd.read_csv(D / "cpjump1/JUMP-Target-1_compound_metadata.tsv", sep="\t")
    pm = pd.read_csv(D / "cpjump1/JUMP-Target-1_compound_platemap.tsv", sep="\t").dropna(subset=["broad_sample"])
    wells = pm.broad_sample.value_counts()
    return pd.DataFrame(dict(
        dataset="CPJUMP1_cpg0000_A549+U2OS", cpd_id=meta.broad_sample, smiles=meta.smiles, inchikey=meta.InChIKey,
        inchikey14=None, name=meta.pert_iname, structure_source="JUMP-Target-1_compound_metadata",
        doses_uM="5", n_wells=meta.broad_sample.map(wells).fillna(0).astype(int).values,
        cell_lines="A549;U2OS", time_h="24;48", pert_type=meta.pert_type, control_type=meta.control_type))


def load_euos():
    ext = pd.read_csv(D / "euos/FMP/workspace/metadata/external_metadata/external_metadata.tsv", sep="\t", index_col=0)
    ext = ext.set_index("Metadata_EOS")
    fs = glob.glob(str(D / "euos/*/workspace/metadata/platemaps/*/platemap/*.txt"))
    pm = pd.concat([pd.read_csv(f, sep="\t").assign(
        cell="U2OS" if "U2OS" in f else "HepG2") for f in fs])
    rows = []
    for eos, g in pm.groupby("Metadata_EOS"):
        if eos == "DMSO":
            continue
        e = ext.loc[eos] if eos in ext.index else None
        sites = sorted({f"{p}:{c}" for p, c in g[["Metadata_Partner", "cell"]].drop_duplicates().itertuples(index=False)})
        rows.append(dict(
            dataset="EU-OS_cpg0036_HepG2+U2OS", cpd_id=eos,
            smiles=e["smiles"] if e is not None else None, inchikey=e["inchikey"] if e is not None else None,
            inchikey14=None, name=None if e is not None else eos,  # only the 2 positive controls have names
            structure_source="external_metadata" if e is not None else "control_name_only",
            doses_uM=fmt_doses(set(g.Metadata_Concentration.astype(float))),
            n_wells=len(g), n_wells_hepg2=int((g.cell == "HepG2").sum()), n_wells_u2os=int((g.cell == "U2OS").sum()),
            cell_lines="HepG2;U2OS", time_h="24", sites=";".join(sites), pdid=e["pdid"] if e is not None else None))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- matching
def prep(ds):
    ds = ds.copy()
    # strip CXSMILES extensions (" |...|", e.g. anisodamine in the Hub file) that RDKit's plain parser rejects
    mols = ds.smiles.map(lambda x: std_mol(x.split(" |")[0]) if isinstance(x, str) else None)
    ds["std_key"] = [Chem.MolToInchiKey(m) or None if m is not None else None for m in mols]
    ds["formula"] = [Chem.MolToInchi(m).split("/")[1] if m is not None else None for m in mols]
    ds["ghash"] = [graph_hash(m) if m is not None else None for m in mols]
    ds["nname"] = ds.name.map(norm_name)
    return ds


def match(tahoe, ds):
    full, b1, by_f, names = {}, {}, {}, {}
    for r in ds.itertuples():
        for k in {r.std_key, r.inchikey} - {None}:
            if isinstance(k, str) and len(k) == 27:
                full.setdefault(k, set()).add(r.cpd_id)
                b1.setdefault(k[:14], set()).add(r.cpd_id)
        if isinstance(r.inchikey14, str):
            b1.setdefault(r.inchikey14, set()).add(r.cpd_id)
        if r.formula:
            by_f.setdefault((r.formula, r.ghash), set()).add(r.cpd_id)
        if r.nname:
            names.setdefault(r.nname, set()).add(r.cpd_id)
    out = []
    for t in tahoe.itertuples():
        keys = [k for k in (t.tahoe_key, t.pubchem_key) if isinstance(k, str)]
        hit, tier = set().union(*[full.get(k, set()) for k in keys]) if keys else set(), None
        if hit:
            tier = "1a_exact_achiral" if t.has_stereo is False else "1b_inchikey_stereo_unresolved"
        if not hit and keys:
            hit = set().union(*[b1.get(k[:14], set()) for k in keys])
            tier = "2_inchikey_block1" if hit else None
        if not hit:
            q = std_mol(t.canonical_smiles)
            if q is not None:
                hit = by_f.get((Chem.MolToInchi(q).split("/")[1], graph_hash(q)), set())
                tier = "2b_same_graph_tautomer_or_bond_variant" if hit else None
        if hit:
            m = Chem.MolFromSmiles(t.canonical_smiles)
            if m is not None and any(a.GetSymbol() in METALS for a in m.GetAtoms()):
                tier = "3_metal_stripped"
        if not hit:
            hit = names.get(norm_name(t.drug), set())
            tier = "3_name_bridge" if hit else None
        if hit:
            out.append((t.drug, tier, sorted(hit)))
    return out


def main():
    tahoe = pd.read_csv(ROOT / "overlap.csv")
    tahoe["drug"] = tahoe["drug"].str.strip()
    dm = pd.read_parquet(ROOT / "data/tahoe/drug_metadata.parquet")
    dm["drug"] = dm["drug"].str.strip()
    tahoe = tahoe.merge(dm[["drug", "canonical_smiles"]], on="drug", how="left")
    tahoe["has_stereo"] = tahoe["has_stereo"].map({True: True, False: False, "True": True, "False": False})
    step2 = pd.read_csv(ROOT / "step2/compounds.csv")
    step2_drugs = {d.strip() for s in step2.tahoe_drugs.dropna() for d in s.split(";")}
    step1_conf = set(tahoe.loc[tahoe.tier.fillna("").str[0].isin(["1", "2"]), "drug"])

    datasets = {"LINCS": load_lincs(), "LINCS_Batch2": load_lincs_batch2(), "CPJUMP1": load_cpjump1(),
                "EU-OS": load_euos()}
    rows, summary = [], {}
    for key, ds in datasets.items():
        ds = prep(ds)
        ds.to_csv(D / f"dataset_compounds_{key}.csv", index=False)
        info = ds.set_index("cpd_id")
        res = match(tahoe, ds)
        for drug, tier, ids in res:
            sub = info.loc[ids]
            t = tahoe.set_index("drug").loc[drug]
            doses = sorted({float(x) for s in sub.doses_uM for x in str(s).split(";")})
            # nearest dataset dose to each Tahoe dose, with fold difference (dataset / Tahoe)
            near = {td: min(doses, key=lambda d: abs(np.log(d / td))) for td in TAHOE_DOSES} if doses else {}
            rows.append(dict(
                drug=drug, dataset=ds.dataset.iloc[0], tier=tier, dataset_cpd_ids=";".join(ids),
                dataset_names=";".join(sorted({str(n) for n in sub.name.dropna()})),
                dataset_std_keys=";".join(sorted({k for k in sub.std_key.dropna()})),
                n_dataset_std_keys=sub.std_key.dropna().nunique(),
                dataset_doses_uM=fmt_doses(doses), tahoe_doses_uM="0.05;0.5;5",
                tahoe_dose_nearest=";".join(f"{td:g}->{d:g} ({d / td:.2g}x)" for td, d in near.items()),
                tahoe_doses_within_1p5x=";".join(f"{td:g}" for td, d in near.items() if 1 / 1.52 <= d / td <= 1.52) or "none",
                tahoe_doses_within_2x=";".join(f"{td:g}" for td, d in near.items() if 0.5 <= d / td <= 2) or "none",
                dataset_wells=int(sub.n_wells.sum()), cell_lines=sub.cell_lines.iloc[0], time_h=sub.time_h.iloc[0],
                is_control=bool(("control_type" in sub and sub.control_type.notna().any())
                                or sub.structure_source.eq("control_name_only").any()),
                step1_jump_tier=t["tier"], step1_jump_confident=drug in step1_conf, in_step2=drug in step2_drugs,
                moa_fine=t["moa-fine"], name_agrees=(norm_name(drug) in {norm_name(n) for n in sub.name.dropna()})
                if sub.name.notna().any() else None))
        m = pd.DataFrame([r for r in rows if r["dataset"] == ds.dataset.iloc[0]])
        conf = m[m.tier.str[0].isin(["1", "2"])].copy()
        # unique matched compounds: collapse stereoisomers / salts / dataset duplicates by std_key (else by id)
        conf["ucpd"] = conf.dataset_std_keys.where(conf.dataset_std_keys != "", conf.dataset_cpd_ids)
        moa = (conf[conf.moa_fine != "unclear"].groupby("moa_fine")["ucpd"].nunique().sort_values(ascending=False))
        summary[key] = {
            "dataset": ds.dataset.iloc[0],
            "dataset_compound_ids": int(ds.cpd_id.nunique()),
            "dataset_unique_structures_std_key": int(ds.std_key.nunique()),
            "dataset_ids_without_structure": int(ds.std_key.isna().sum()),
            "structure_source_counts": ds.structure_source.fillna("unresolved").value_counts().to_dict(),
            "dose_levels_uM": ds.doses_uM.value_counts().to_dict(),
            "matched_tahoe_drugs_by_tier": m.tier.value_counts().to_dict() if len(m) else {},
            "matched_tahoe_drugs_t12": len(conf),
            "matched_unique_cpds_t12": int(conf.ucpd.nunique()),
            "matched_t12_also_step1_jump_t12": int(conf.step1_jump_confident.sum()),
            "matched_t12_in_step2": int(conf.in_step2.sum()),
            "matched_t12_controls": conf.loc[conf.is_control, "drug"].tolist(),
            "matched_t12_not_in_step1": conf.loc[~conf.step1_jump_confident, ["drug", "tier", "step1_jump_tier"]].values.tolist(),
            "tier3": m.loc[m.tier.str.startswith("3"), ["drug", "tier", "dataset_cpd_ids", "dataset_names"]].values.tolist(),
            "name_disagrees_t12": conf.loc[conf.name_agrees == False, ["drug", "dataset_names"]].values.tolist(),  # noqa: E712
            "multi_std_key_hits": conf.loc[conf.n_dataset_std_keys > 1, ["drug", "dataset_cpd_ids"]].values.tolist(),
            "moa_classes_ge3_unique": moa[moa >= 3].to_dict(),
            "tahoe_doses_within_1p5x_t12": conf.tahoe_doses_within_1p5x.value_counts().to_dict(),
            "tahoe_doses_within_2x_t12": conf.tahoe_doses_within_2x.value_counts().to_dict(),
            "tahoe_dose_nearest_t12": conf.tahoe_dose_nearest.value_counts().head(5).to_dict(),
        }
    out = pd.DataFrame(rows).sort_values(["dataset", "tier", "drug"])
    out.to_csv(ROOT / "cellcheck/cellline_matches.csv", index=False)
    (D / "summary.json").write_text(json.dumps(summary, indent=1, default=str))
    print(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()
