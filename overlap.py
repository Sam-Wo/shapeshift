"""ShapeShift step 1: Tahoe-100M x JUMP-CP compound overlap (metadata only).

Usage:  py -3.13 overlap.py
Downloads small metadata tables into ./data (never expression data, images or profiles),
standardizes structures with RDKit, matches in three tiers and writes overlap.csv + stats
used by REPORT.md (printed to stdout and saved to data/summary.json).
"""
import ast
import json
import time
from pathlib import Path

import pandas as pd
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolHash
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).parent
DATA = ROOT / "data"

HF = "https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/metadata/"
JUMP = "https://raw.githubusercontent.com/jump-cellpainting/datasets/main/metadata/"
FILES = {
    "tahoe/drug_metadata.parquet": HF + "drug_metadata.parquet",
    "tahoe/sample_metadata.parquet": HF + "sample_metadata.parquet",  # 65 kB; NOT obs_metadata (2.3 GB)
    "tahoe/cell_line_metadata.parquet": HF + "cell_line_metadata.parquet",  # 19 kB, cell-line check
    "jump/compound.csv.gz": JUMP + "compound.csv.gz",
    "jump/well.csv.gz": JUMP + "well.csv.gz",
    "jump/plate.csv.gz": JUMP + "plate.csv.gz",
    "jump/perturbation_control.csv": JUMP + "perturbation_control.csv",
    # JUMP team's JCP2022-keyed join of the Broad Drug Repurposing Hub (names, MoA, targets; 287 kB)
    "jump/repurposing_hub_annotations_processed.tsv":
        "https://cellpainting-gallery.s3.amazonaws.com/cpg0016-jump/source_all/workspace/publication_data/"
        "2025_Chandrasekaran/jump_production_datastore/interim/repurposing_hub_annotations_processed.tsv",
}

# Nominal JUMP cpg0016 concentrations; no per-well dose is published. Source:
# https://github.com/broadinstitute/jump_hub/blob/main/explanations/quirks_details.md
def jump_conc_uM(source, plate_type, is_poscon):
    if is_poscon or plate_type in ("TARGET1", "TARGET2", "POSCON8"):
        return 5.0
    if plate_type == "COMPOUND":
        return 0.625 if source == "source_7" else 10.0
    return None


def norm_name(s):
    """Lower-case drug name without salt / hydrate words and punctuation, for name bridging."""
    if not isinstance(s, str):
        return None
    s = s.lower().replace("(", " ").replace(")", " ").replace(",", " ")
    return " ".join(w for w in s.split() if w not in SALT_WORDS).replace("-", "").replace(" ", "") or None
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PC_CACHE = DATA / "pubchem_cache.json"


# ---------------------------------------------------------------- download
def fetch_all():
    for rel, url in FILES.items():
        p = DATA / rel
        if not p.exists():
            p.parent.mkdir(parents=True, exist_ok=True)
            print("downloading", url)
            r = requests.get(url, timeout=120)
            r.raise_for_status()
            p.write_bytes(r.content)


# ---------------------------------------------------------------- standardization
_cleanup = rdMolStandardize.CleanupParameters()
_chooser = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
_uncharger = rdMolStandardize.Uncharger()


def std_mol(smiles):
    """SMILES -> neutralized, stereo-free parent (largest organic fragment), or None."""
    if not isinstance(smiles, str) or not smiles:
        return None
    try:
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return None
        m = rdMolStandardize.Cleanup(m, _cleanup)
        m = _chooser.choose(m)
        m = _uncharger.uncharge(m)
        Chem.RemoveStereochemistry(m)
        return m
    except Exception:
        return None


def std_key(smiles):
    """Stereo-free InChIKey of the standardized parent.
    Stereo is removed because JUMP's InChIKeys are stereo-flattened (~100% '-UHFFFAOYSA-')."""
    m = std_mol(smiles)
    return (Chem.MolToInchiKey(m) or None) if m is not None else None


def graph_hash(m):
    """Element + connectivity hash that ignores bond orders, charges and H placement."""
    return rdMolHash.MolHash(Chem.Mol(m), rdMolHash.HashFunction.ElementGraph)


# ---------------------------------------------------------------- PubChem (cached)
def _pc_load():
    return json.loads(PC_CACHE.read_text()) if PC_CACHE.exists() else {}


def _pc_get(cache, key, url, data=None):
    if key in cache:
        return cache[key]
    out = None
    for attempt in range(4):
        try:
            r = requests.post(url, data=data, timeout=60) if data else requests.get(url, timeout=60)
            if r.status_code == 200:
                out = r.json()
                break
            if r.status_code == 404:  # genuine "not found": cache it as empty
                cache[key] = {}
                return {}
        except Exception:
            pass
        time.sleep(2 * (attempt + 1))  # transient error / throttling (503): back off, don't cache
    if out is not None:
        cache[key] = out
    time.sleep(0.25)  # PubChem asks for <= 5 req/s
    return out


def pubchem_props(cids, cache):
    """CID -> {InChIKey, IsomericSMILES, Title}."""
    out = {}
    cids = [str(c) for c in cids]
    for i in range(0, len(cids), 100):
        chunk = cids[i:i + 100]
        j = _pc_get(cache, "props:" + ",".join(chunk),
                    f"{PUBCHEM}/compound/cid/property/InChIKey,IsomericSMILES,SMILES,Title/JSON",
                    data={"cid": ",".join(chunk)})
        for p in (j or {}).get("PropertyTable", {}).get("Properties", []):
            out[str(p["CID"])] = p
    return out


def pubchem_parent(cid, cache):
    j = _pc_get(cache, f"parent:{cid}", f"{PUBCHEM}/compound/cid/{cid}/cids/JSON?cids_type=parent")
    return [str(c) for c in (j or {}).get("IdentifierList", {}).get("CID", [])]


def pubchem_name(name, cache):
    j = _pc_get(cache, f"name:{name.lower()}",
                f"{PUBCHEM}/compound/name/{requests.utils.quote(name)}/cids/JSON")
    return [str(c) for c in (j or {}).get("IdentifierList", {}).get("CID", [])][:3]


# ---------------------------------------------------------------- main
def main():
    fetch_all()
    tahoe = pd.read_parquet(DATA / "tahoe/drug_metadata.parquet")
    tahoe["drug"] = tahoe["drug"].str.strip()
    samples = pd.read_parquet(DATA / "tahoe/sample_metadata.parquet")
    comp = pd.read_csv(DATA / "jump/compound.csv.gz")
    well = pd.read_csv(DATA / "jump/well.csv.gz")
    plate = pd.read_csv(DATA / "jump/plate.csv.gz")
    controls = pd.read_csv(DATA / "jump/perturbation_control.csv")

    # --- Tahoe doses (from the 65 kB sample table)
    doses = {}
    for v in samples["drugname_drugconc"].dropna().unique():
        for name, conc, unit in ast.literal_eval(v):
            doses.setdefault(name.strip(), set()).add((float(conc), unit))
    tahoe["tahoe_doses"] = tahoe["drug"].map(
        lambda d: ";".join(f"{c:g} {u}" for c, u in sorted(doses.get(d, ()))) or None)

    # --- JUMP: standardize every SMILES (also keep JUMP's own key)
    std_cache = DATA / "jump_std_keys.csv.gz"
    if std_cache.exists():
        comp = comp.merge(pd.read_csv(std_cache), on="Metadata_JCP2022", how="left")
    else:
        print(f"standardizing {len(comp):,} JUMP SMILES (one-off, cached) ...")
        comp["std_key"] = comp["Metadata_SMILES"].map(std_key)
        comp[["Metadata_JCP2022", "std_key"]].to_csv(std_cache, index=False)
    jump_keys = {}  # full key -> JCP ids ; block1 -> JCP ids
    jump_b1 = {}
    for jcp, k_own, k_std in comp[["Metadata_JCP2022", "Metadata_InChIKey", "std_key"]].itertuples(index=False):
        for k in {k_own, k_std} - {None, float("nan")}:
            if isinstance(k, str):
                jump_keys.setdefault(k, set()).add(jcp)
                jump_b1.setdefault(k[:14], set()).add(jcp)

    # --- Tahoe structures: Tahoe SMILES + PubChem (isomeric) structure of the CID
    cache = _pc_load()
    tahoe["cid"] = tahoe["pubchem_cid"].astype("Int64")
    props = pubchem_props(tahoe["cid"].dropna().unique(), cache)
    if len(props) < 0.9 * tahoe["cid"].notna().sum():
        raise RuntimeError(f"PubChem returned only {len(props)} CIDs (transient outage?) - rerun later")
    tahoe["tahoe_key"] = tahoe["canonical_smiles"].map(std_key)
    tahoe["pubchem_smiles"] = tahoe["cid"].map(
        lambda c: (props.get(str(int(c)), {}).get("IsomericSMILES") or props.get(str(int(c)), {}).get("SMILES"))
        if pd.notna(c) else None)
    tahoe["pubchem_key"] = tahoe["pubchem_smiles"].map(std_key)
    # stereo status from PubChem's isomeric SMILES (Tahoe's own SMILES are flat)
    tahoe["has_stereo"] = tahoe["pubchem_smiles"].map(
        lambda s: any(ch in s for ch in "@/\\") if isinstance(s, str) else None)
    tahoe["pubchem_title"] = tahoe["cid"].map(lambda c: props.get(str(int(c)), {}).get("Title") if pd.notna(c) else None)

    tier, jcps, via = [], [], []
    for r in tahoe.itertuples():
        keys = [k for k in (r.tahoe_key, r.pubchem_key) if k]
        hit = set().union(*[jump_keys.get(k, set()) for k in keys]) if keys else set()
        if hit:
            # both sides stereo-flat: only achiral drugs are a truly exact match
            tier.append("1a_exact_achiral" if r.has_stereo is False else "1b_inchikey_stereo_unresolved")
            jcps.append(hit); via.append("structure"); continue
        hit = set().union(*[jump_b1.get(k[:14], set()) for k in keys]) if keys else set()
        if hit:
            tier.append("2_inchikey_block1"); jcps.append(hit); via.append("structure"); continue
        tier.append(None); jcps.append(set()); via.append(None)
    tahoe["tier"], tahoe["jcp_set"], tahoe["match_via"] = tier, jcps, via

    # --- Tier 2b: same connectivity drawn differently. JUMP stores many drugs as unusual tautomers or
    # with shifted double bonds (e.g. sunitinib as a hydroxy-indole, dexamethasone as a delta-5 diene,
    # lactic acid as CC(=O)C(O)O), which changes even InChIKey block 1. Match on identical molecular
    # formula + identical heavy-atom graph (elements and bonds, ignoring bond orders / H placement).
    print("tier 2b: same formula + same heavy-atom graph ...")
    comp["formula"] = comp["Metadata_InChI"].str.split("/").str[1]
    by_formula = comp.groupby("formula")
    for i in tahoe.index[tahoe["tier"].isna()]:
        q = std_mol(tahoe.at[i, "canonical_smiles"])
        if q is None:
            continue
        f = Chem.MolToInchi(q).split("/")[1]
        if f not in by_formula.groups:
            continue
        hq = graph_hash(q)
        hit = {jcp for jcp, s in by_formula.get_group(f)[["Metadata_JCP2022", "Metadata_SMILES"]]
               .itertuples(index=False) if (m := std_mol(s)) is not None and graph_hash(m) == hq}
        if hit:
            tahoe.at[i, "tier"] = "2b_same_graph_tautomer_or_bond_variant"
            tahoe.at[i, "jcp_set"] = hit
            tahoe.at[i, "match_via"] = f"formula {f} + heavy-atom graph - REVIEW (structure-only)"

    # --- Demote matches where standardization discarded a metal (e.g. auranofin -> its thiosugar
    # ligand): the JUMP compound may be the bare ligand, not the drug.
    for i in tahoe.index[tahoe["tier"].notna()]:
        full = Chem.MolFromSmiles(tahoe.at[i, "canonical_smiles"])
        if full is not None and any(a.GetSymbol() in METALS for a in full.GetAtoms()):
            tahoe.at[i, "tier"] = "3_metal_stripped"
            tahoe.at[i, "match_via"] = f"{tahoe.at[i, 'match_via']}; metal removed by standardization - REVIEW"

    # --- Tier 3a: name bridge via JUMP-side names (Repurposing Hub join + JUMP control names)
    hub = pd.read_csv(DATA / "jump/repurposing_hub_annotations_processed.tsv", sep="\t")
    jump_names = {}
    for jcp, nm in pd.concat([
            hub[["Metadata_JCP2022", "Metadata_repurposing_name"]].set_axis(["jcp", "name"], axis=1),
            controls[["Metadata_JCP2022", "Metadata_Name"]].set_axis(["jcp", "name"], axis=1)]).itertuples(index=False):
        if norm_name(nm):
            jump_names.setdefault(norm_name(nm), set()).add(jcp)
    for i in tahoe.index[tahoe["tier"].isna()]:
        r = tahoe.loc[i]
        for n in (r["drug"], r["pubchem_title"]):
            hit = jump_names.get(norm_name(n))
            if hit:
                tahoe.at[i, "tier"] = "3_name_bridge"
                tahoe.at[i, "jcp_set"] = set(hit)
                tahoe.at[i, "match_via"] = f"name '{n}' = JUMP/Repurposing Hub name - REVIEW"
                break

    # --- Tier 3b: PubChem structure bridge for unmatched (parent CID of salts; name/synonym search)
    print("tier 3b: PubChem parent / name lookups for unmatched drugs ...")
    for i in tahoe.index[tahoe["tier"].isna()]:
        r = tahoe.loc[i]
        cand = []
        if pd.notna(r["cid"]):
            cand += [(c, "pubchem_parent_cid") for c in pubchem_parent(r["cid"], cache) if c != str(r["cid"])]
        names = {r["drug"]}
        if isinstance(r["pubchem_title"], str):
            names.add(r["pubchem_title"])
        # strip common salt words from the Tahoe name, e.g. "Ixazomib citrate" -> "Ixazomib"
        base = " ".join(w for w in r["drug"].split() if w.lower() not in SALT_WORDS)
        names.add(base)
        for n in names:
            cand += [(c, f"pubchem_name:{n}") for c in pubchem_name(n, cache)]
        if not cand:
            continue
        cprops = pubchem_props(sorted({c for c, _ in cand}), cache)
        for c, how in cand:
            k = std_key(cprops.get(c, {}).get("IsomericSMILES") or cprops.get(c, {}).get("SMILES"))
            if k and (jump_b1.get(k[:14])):
                tahoe.at[i, "tier"] = "3_pubchem_structure"
                tahoe.at[i, "jcp_set"] = set(jump_b1[k[:14]])
                tahoe.at[i, "match_via"] = f"{how} (CID {c}) - REVIEW"
                break
    PC_CACHE.write_text(json.dumps(cache))


    # --- JUMP well / plate / source counts per JCP, split by plate type
    w = well.merge(plate[["Metadata_Source", "Metadata_Plate", "Metadata_PlateType"]],
                   on=["Metadata_Source", "Metadata_Plate"], how="left")
    ctrl_ids = set(controls.loc[controls["Metadata_pert_type"] == "poscon", "Metadata_JCP2022"])  # JUMP positive controls
    agg = {}
    for jcp_ids in tahoe["jcp_set"]:
        for j in jcp_ids:
            agg[j] = None
    wsub = w[w["Metadata_JCP2022"].isin(agg)]

    def summarize(ids):
        s = wsub[wsub["Metadata_JCP2022"].isin(ids)]
        comp_s = s[s["Metadata_PlateType"] == "COMPOUND"]
        other = s[s["Metadata_PlateType"] != "COMPOUND"]
        return pd.Series({
            "jump_wells_compound": len(comp_s),
            "jump_plates_compound": comp_s["Metadata_Plate"].nunique(),
            "jump_sources_compound": comp_s["Metadata_Source"].nunique(),
            "jump_source_list": ";".join(sorted(comp_s["Metadata_Source"].unique())),
            "jump_wells_10uM": int((comp_s["Metadata_Source"] != "source_7").sum()),
            "jump_wells_0.625uM_source7": int((comp_s["Metadata_Source"] == "source_7").sum()),
            "jump_wells_other_platetypes": len(other),
            "jump_other_platetypes": ";".join(sorted(other["Metadata_PlateType"].dropna().unique())),
            "jump_is_poscon": bool(set(ids) & ctrl_ids),
            "jump_concentration_uM_nominal": ";".join(f"{c:g}" for c in sorted({
                c for c in (jump_conc_uM(src, pt, bool(set(ids) & ctrl_ids))
                            for src, pt in s[["Metadata_Source", "Metadata_PlateType"]].drop_duplicates()
                            .itertuples(index=False)) if c is not None})) or None,
            "jump_hub_name": ";".join(sorted(set(hub.loc[hub.Metadata_JCP2022.isin(ids), "Metadata_repurposing_name"]))) or None,
            "jump_hub_moa": ";".join(sorted(set(hub.loc[hub.Metadata_JCP2022.isin(ids), "Metadata_repurposing_moa"].dropna()))) or None,
        })

    stats = tahoe["jcp_set"].apply(lambda ids: summarize(ids) if ids else pd.Series(dtype=object))
    tahoe = pd.concat([tahoe, stats], axis=1)
    tahoe["jump_jcp_ids"] = tahoe["jcp_set"].map(lambda s: ";".join(sorted(s)) or None)
    # independent sanity check: does the JUMP-side Repurposing Hub name agree with the Tahoe name?
    tahoe["hub_name_agrees"] = [
        (norm_name(d) in {norm_name(x) for x in hn.split(";")} or
         (isinstance(pt, str) and norm_name(pt) in {norm_name(x) for x in hn.split(";")}))
        if isinstance(hn, str) else None
        for d, pt, hn in tahoe[["drug", "pubchem_title", "jump_hub_name"]].itertuples(index=False)]

    # --- outputs
    cols = ["drug", "tier", "match_via", "jump_jcp_ids", "moa-broad", "moa-fine", "targets", "human-approved",
            "cid", "has_stereo", "tahoe_key", "pubchem_key", "jump_wells_compound", "jump_plates_compound",
            "jump_sources_compound", "jump_source_list", "jump_wells_10uM", "jump_wells_0.625uM_source7",
            "jump_wells_other_platetypes", "jump_other_platetypes",
            "jump_is_poscon", "jump_concentration_uM_nominal", "tahoe_doses", "jump_hub_name", "hub_name_agrees",
            "jump_hub_moa"]
    out = tahoe[cols].sort_values(["tier", "drug"], na_position="last")
    out.to_csv(ROOT / "overlap.csv", index=False)

    m = tahoe[tahoe["tier"].notna()]
    conf = m[m["tier"].str[0].isin(["1", "2"])]
    # unique compounds: Tahoe salt/hydrate duplicates and stereoisomers (JUMP is stereo-flat) collapse
    conf = conf.assign(cpd=conf["jump_jcp_ids"])
    unique_parents_conf = conf["cpd"].nunique()
    t1 = conf[conf["tier"].str.startswith("1")]
    moa_tab = (conf[conf["moa-fine"] != "unclear"].groupby("moa-fine")
               .agg(matched_drugs=("drug", "size"), matched_unique_jump_cpds=("cpd", "nunique"))
               .rename_axis("moa_fine").reset_index())
    moa_tab = moa_tab.merge(tahoe["moa-fine"].value_counts().rename_axis("moa_fine").reset_index(name="tahoe_total"),
                            on="moa_fine", how="right").fillna({"matched_drugs": 0, "matched_unique_jump_cpds": 0})
    moa_tab = moa_tab[moa_tab.moa_fine != "unclear"].sort_values(
        ["matched_unique_jump_cpds", "tahoe_total"], ascending=False)
    summary = {
        "tahoe_drugs": len(tahoe),
        "tahoe_with_doses": int(tahoe["tahoe_doses"].notna().sum()),
        "tahoe_dose_count_dist": tahoe["tahoe_doses"].fillna("").str.count(";").add(1).where(tahoe["tahoe_doses"].notna()).value_counts().to_dict(),
        "tahoe_parse_fail": tahoe.loc[tahoe["tahoe_key"].isna(), "drug"].tolist(),
        "tahoe_vs_pubchem_key_block1_disagree": tahoe.loc[
            tahoe.tahoe_key.notna() & tahoe.pubchem_key.notna() & (tahoe.tahoe_key.str[:14] != tahoe.pubchem_key.str[:14]), "drug"].tolist(),
        "jump_compounds": len(comp),
        "jump_std_parse_fail": int(comp["std_key"].isna().sum()),
        "jump_own_vs_std_key_differ": int((comp["std_key"].notna() & (comp["std_key"] != comp["Metadata_InChIKey"])).sum()),
        "tier_counts": tahoe["tier"].fillna("unmatched").value_counts().to_dict(),
        "confident_t12_drugs": len(conf),
        "confident_t12_unique_jump_compounds": int(unique_parents_conf),
        "confident_t12_shared_jump_cpd": conf[conf.duplicated("cpd", keep=False)]
            .groupby("cpd")["drug"].apply(list).tolist(),
        "confident_hub_name_agrees": conf["hub_name_agrees"].value_counts(dropna=False).to_dict(),
        "confident_hub_name_disagrees": conf.loc[conf["hub_name_agrees"] == False, ["drug", "jump_hub_name"]].values.tolist(),
        "confident_on_compound_plates": int((conf["jump_wells_compound"] > 0).sum()),
        "confident_control_or_target2_only": conf.loc[conf["jump_wells_compound"] == 0, "drug"].tolist(),
        "confident_with_unclear_moa": int((conf["moa-fine"] == "unclear").sum()),
        "moa_classes_ge3_matched_unique": int((moa_tab["matched_unique_jump_cpds"] >= 3).sum()),
        "wells_compound_median": float(conf["jump_wells_compound"].median()),
        "sources_compound_dist": conf["jump_sources_compound"].value_counts().sort_index().to_dict(),
        "tier3": m.loc[m["tier"].str.startswith("3"), ["drug", "match_via", "jump_jcp_ids", "jump_hub_name"]].values.tolist(),
        "confident_unclear_moa_but_hub_moa": int(((conf["moa-fine"] == "unclear") & conf["jump_hub_moa"].notna()).sum()),
        "jump_conc_dist_confident": conf["jump_concentration_uM_nominal"].value_counts().to_dict(),
        **{f"dose_subset_{name}": {
            "drugs": int(sel.sum()),
            "unique_jump_cpds": int(conf.loc[sel, "cpd"].nunique()),
            "median_wells": float(conf.loc[sel, col].median()),
            "moa_classes_ge3_unique": int((conf[sel & (conf["moa-fine"] != "unclear")]
                                           .groupby("moa-fine")["cpd"].nunique() >= 3).sum()),
        } for name, col, sel in [  # positive controls excluded: their wells are 5 uM control wells
            ("10uM", "jump_wells_10uM", (conf["jump_wells_10uM"] > 0) & ~conf["jump_is_poscon"].fillna(False).astype(bool)),
            ("0.625uM_source7", "jump_wells_0.625uM_source7", (conf["jump_wells_0.625uM_source7"] > 0) & ~conf["jump_is_poscon"].fillna(False).astype(bool)),
            ("both", "jump_wells_10uM", (conf["jump_wells_10uM"] > 0) & (conf["jump_wells_0.625uM_source7"] > 0) & ~conf["jump_is_poscon"].fillna(False).astype(bool)),
        ]},
        "tier1_only": {
            "drugs": int(t1.shape[0]),
            "unique_jump_cpds": int(t1["cpd"].nunique()),
            "moa_classes_ge3_unique": int((t1[t1["moa-fine"] != "unclear"].groupby("moa-fine")["cpd"].nunique() >= 3).sum()),
            "achiral_exact": int((t1["tier"] == "1a_exact_achiral").sum()),
        },
        "unmatched_approved_with_moa": tahoe.loc[tahoe["tier"].isna() & (tahoe["human-approved"] == "yes"),
                                                 ["drug", "moa-fine"]].values.tolist(),
    }
    moa_tab.to_csv(DATA / "moa_coverage.csv", index=False)
    (DATA / "summary.json").write_text(json.dumps(summary, indent=1, default=str))
    print(json.dumps(summary, indent=1, default=str))
    print(moa_tab.to_string(index=False))


METALS = {"Pt", "Pd", "Au", "Ag", "As", "Sb", "Bi", "Hg", "Fe", "Cu", "Zn", "Co", "Ni", "Ru", "Rh", "Ir", "Os",
          "Ga", "Gd", "Ti", "V", "Cr", "Mn", "Sn", "Se", "Te", "Tl", "Pb", "Cd", "Mo", "W", "Re", "Tc"}

SALT_WORDS = {"hydrochloride", "hcl", "dihydrochloride", "sodium", "potassium", "citrate", "mesylate", "maleate",
              "tosylate", "sulfate", "phosphate", "acetate", "fumarate", "tartrate", "succinate", "besylate",
              "hydrobromide", "bromide", "chloride", "malate", "lactate", "dimesylate", "monohydrate", "hydrate",
              "(s)-", "(r)-", "calcium", "magnesium", "trihydrate", "disodium", "dihydrate", "free", "base"}

if __name__ == "__main__":
    main()
