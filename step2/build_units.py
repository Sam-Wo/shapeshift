"""Build the shared compound-unit table for the three-tier Step 2.

Unit = Tahoe-side standardized parent structure (overlap.py `std_key` of the Tahoe SMILES, stereo-flat),
so salt / hydrate variants of one Tahoe drug collapse, consistently across tiers. Matching itself comes
from Step 1 (overlap.csv, JUMP) and cellcheck/cellline_matches.csv (LINCS, EU-OS, CPJUMP1), both produced
with overlap.py's matching functions - nothing is re-matched here.

Output: step2/units.csv (one row per unit x tier) and step2/units_summary.json
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
TIER_DATASET = {"T1_LINCS_A549": "LINCS_cpg0004_A549_48h",
                "T2_EUOS_HepG2": "EU-OS_cpg0036_HepG2+U2OS",
                "CAL_CPJUMP1": "CPJUMP1_cpg0000_A549+U2OS"}
DROP_DRUGS = {"Auranofin"}                       # JUMP/any match lacks the gold centre
DROP_PAIRS = {("Estrone sulfate (potassium)", "LINCS")}  # estrone sulfate -> estrone: different molecule


def first_moa(s):
    s = [x for x in s if isinstance(x, str) and x != "unclear"]
    return s[0] if s else "unclear"


def make_t3_compounds():
    """T3 working set (step2/compounds.csv): Step-1 confident JUMP matches (tier 1/2) measured at BOTH JUMP doses
    (10 uM non-source_7 and 0.625 uM source_7), excluding positive controls and auranofin; one row per JCP id."""
    t = pd.read_csv(ROOT / "overlap.csv")
    conf = t[t.tier.fillna("").str[0].isin(["1", "2"])]
    pos = conf.jump_is_poscon.astype(str).str.lower().eq("true")
    b = conf[(conf.jump_wells_10uM > 0) & (conf["jump_wells_0.625uM_source7"] > 0) & ~pos & (conf.drug != "Auranofin")]
    g = (b.groupby("jump_jcp_ids")
         .agg(tahoe_drugs=("drug", lambda s: ";".join(s)), n_tahoe_drugs=("drug", "size"),
              moa=("moa-fine", first_moa),
              moa_conflict=("moa-fine", lambda s: len(set(x for x in s if x != "unclear")) > 1),
              tiers=("tier", lambda s: ";".join(sorted(set(s)))), hub_moa=("jump_hub_moa", "first"),
              wells_10uM=("jump_wells_10uM", "first"), wells_0625=("jump_wells_0.625uM_source7", "first"))
         .reset_index().rename(columns={"jump_jcp_ids": "compound_id"}))
    g["exact_achiral"] = g.tiers.eq("1a_exact_achiral")
    g.to_csv(ROOT / "step2" / "compounds.csv", index=False)
    return g


def main():
    make_t3_compounds()
    ov = pd.read_csv(ROOT / "overlap.csv")
    ov["unit"] = ov["tahoe_key"].fillna(ov["pubchem_key"])
    drug2unit = dict(zip(ov.drug, ov.unit))
    drug_moa = dict(zip(ov.drug, ov["moa-fine"]))
    rows = []

    # --- T1 / T2 / calibration from the cell-line check
    cm = pd.read_csv(ROOT / "cellcheck" / "cellline_matches.csv")
    for tier, ds in TIER_DATASET.items():
        d = cm[cm.dataset == ds]
        for r in d.itertuples():
            if r.drug in DROP_DRUGS or any(r.drug == a and b in ds for a, b in DROP_PAIRS):
                continue
            if not (str(r.tier).startswith("1") or str(r.tier).startswith("2") or str(r.tier).startswith("3_name")):
                continue
            rows.append(dict(tier=tier, unit=drug2unit[r.drug], tahoe_drug=r.drug, match_tier=r.tier,
                             dataset_ids=r.dataset_cpd_ids, dataset_names=r.dataset_names,
                             dataset_doses_uM=r.dataset_doses_uM, is_control=r.is_control,
                             name_agrees=r.name_agrees, moa=drug_moa.get(r.drug, "unclear")))

    # --- T3: JUMP U2OS, restricted to the 175 both-dose compounds (filter applies to T3 only)
    c3 = pd.read_csv(ROOT / "step2" / "compounds.csv")
    for r in c3.itertuples():
        for drug in r.tahoe_drugs.split(";"):
            if drug in DROP_DRUGS:
                continue
            o = ov[ov.drug == drug].iloc[0]
            rows.append(dict(tier="T3_JUMP_U2OS", unit=drug2unit[drug], tahoe_drug=drug, match_tier=o.tier,
                             dataset_ids=r.compound_id, dataset_names=o.jump_hub_name,
                             dataset_doses_uM="0.625;10", is_control=False,
                             name_agrees=o.hub_name_agrees, moa=drug_moa.get(drug, "unclear")))
    long = pd.DataFrame(rows)

    # stereo / racemate flag: a chiral (1b) match where the names differ only by a stereo prefix
    # (levobupivacaine/bupivacaine, dexmedetomidine/medetomidine, (S)-crizotinib/crizotinib), or the
    # dataset lists two stereoisomers under one structure (quinidine;quinine). Plain synonyms are not flagged.
    import sys
    sys.path.insert(0, str(ROOT))
    from overlap import norm_name
    STEREO_PREFIX = {"levo", "dex", "es", "s", "r", "l", "d", "rac", "dl"}

    def stereo_pair(r):
        if not str(r.match_tier).startswith("1b"):
            return False
        names = [n for n in str(r.dataset_names).replace(";", "|").split("|") if n and n != "nan"]
        normed = {norm_name(n) for n in names} - {None}
        if len(normed) > 1:
            return True
        t = norm_name(r.tahoe_drug)
        for n in normed:
            if n and t and n != t:
                long_, short = (t, n) if len(t) > len(n) else (n, t)
                if short in long_ and long_.replace(short, "", 1) in STEREO_PREFIX:
                    return True
        return False

    long["stereo_flag"] = [stereo_pair(r) for r in long.itertuples()]
    units = (long.groupby(["tier", "unit"])
             .agg(tahoe_drugs=("tahoe_drug", lambda s: ";".join(sorted(set(s)))),
                  n_tahoe_drugs=("tahoe_drug", "nunique"),
                  dataset_ids=("dataset_ids", lambda s: ";".join(sorted(set(";".join(map(str, s)).split(";"))))),
                  dataset_names=("dataset_names", lambda s: ";".join(sorted({str(x) for x in s if pd.notna(x)}))),
                  match_tiers=("match_tier", lambda s: ";".join(sorted(set(s)))),
                  moa=("moa", first_moa),
                  stereo_flag=("stereo_flag", "any"),
                  any_tier3=("match_tier", lambda s: any(str(x).startswith("3") for x in s)),
                  dataset_doses_uM=("dataset_doses_uM", "first"))
             .reset_index())
    units["exact_achiral"] = units.match_tiers.eq("1a_exact_achiral")
    # three-way intersection (T1 & T2 & T3); T2 membership kept even if T2 is later gated out
    sets = {t: set(units.loc[units.tier == t, "unit"]) for t in units.tier.unique()}
    inter3 = sets["T1_LINCS_A549"] & sets["T2_EUOS_HepG2"] & sets["T3_JUMP_U2OS"]
    inter13 = sets["T1_LINCS_A549"] & sets["T3_JUMP_U2OS"]
    units["in_T1T2T3"] = units.unit.isin(inter3)
    units["in_T1T3"] = units.unit.isin(inter13)
    units.to_csv(ROOT / "step2" / "units.csv", index=False)

    summ = {}
    for t, g in units.groupby("tier"):
        vc = g[g.moa != "unclear"].moa.value_counts()
        summ[t] = dict(units=len(g), tahoe_drugs=int(g.n_tahoe_drugs.sum()), exact_achiral=int(g.exact_achiral.sum()),
                       stereo_flag=int(g.stereo_flag.sum()), tier3=int(g.any_tier3.sum()),
                       annotated=int((g.moa != "unclear").sum()), moa_classes_ge3=int((vc >= 3).sum()))
    gi = units[units.in_T1T2T3].drop_duplicates("unit")
    vc = gi[gi.moa != "unclear"].moa.value_counts()
    summ["intersection_T1T2T3"] = dict(units=len(inter3), annotated=int((gi.moa != "unclear").sum()),
                                       moa_classes_ge3=int((vc >= 3).sum()),
                                       exact_achiral=int(units[units.in_T1T2T3].groupby("unit").exact_achiral.all().sum()))
    summ["intersection_T1T3"] = dict(units=len(inter13))
    all_drugs = sorted(set(long.tahoe_drug))
    summ["tahoe_drugs_needed_union"] = len(all_drugs)
    (ROOT / "step2" / "tahoe_drugs_needed.txt").write_text("\n".join(all_drugs))
    (ROOT / "step2" / "units_summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
