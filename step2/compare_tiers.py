"""ShapeShift Step 2 (three-tier design): does morphology agree with the transcriptome on compound similarity?

PRE-REGISTERED PRIMARY TESTS (fixed before any result was seen; everything else is secondary):
  For each tier: annotated units (Tahoe moa-fine != 'unclear'), Spearman Mantel r between the
  cosine-similarity matrices of the morphology consensus and the Tahoe signature, with a null that
  permutes unit labels WITHIN MoA class (singleton classes fixed); 9999 permutations; one-sided p;
  95% CI from a unit bootstrap (duplicate self-pairs dropped).
    T1  LINCS A549      morphology 3.33 uM  vs Tahoe 5 uM  A549 (CVCL_0023) per-line signature
    T2  EU-OS HepG2     morphology 10 uM    vs Tahoe 5 uM  HepG2/C3A (CVCL_1098) per-line signature
    T3  JUMP U2OS       morphology 10 uM    vs Tahoe 5 uM  cross-line consensus signature
  T1 and T2 have orthogonal confounds (time vs exact line) and are reported side by side, never pooled.

Secondary: plain permutation null, all units, every morphology x Tahoe-signature pairing (3x3), rank of
the matched line among all 50, plate-centered Tahoe signatures, pairs sharing no plate, MoA retrieval
(leave-one-out), split-half reliability, T1 dose-response grid + lowest active dose, intersections,
exact-achiral sensitivity, consistency vs agreement, CPJUMP1 attenuation.

Inputs (produced by the pipeline agents; schemas checked, script raises on missing inputs unless a
gate file step2/gate_<tier>.json says the tier is UNDERPOWERED / UNAVAILABLE):
  step2/units.csv, step2/tahoe_perline.parquet, tahoe_signatures.parquet, step2/tahoe_consistency.csv,
  step2/tahoe_activity.csv, step2/morph_lincs*.{parquet,csv}, step2/morph_euos*.{parquet,csv},
  morph_profiles_{high,low}.parquet + step2/morph_wells_{high,low}.parquet + step2/morph_activity.csv (T3),
  step2/cpjump1_calibration.csv
Outputs: step2/results/*.csv, step2/results/nulls.npz, step2/figures/*.png, printed summary.
Usage: py -3.13 -I -X utf8 step2/compare_tiers.py [--perms 9999] [--perms2 1999] [--boot 1000]
"""
import argparse
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parent.parent
S2 = ROOT / "step2"
OUT = S2 / "results"
FIG = S2 / "figures"
RNG = np.random.default_rng(20261122)
A549, HEPG2 = "CVCL_0023", "CVCL_1098"
TAHOE_DOSES = (0.05, 0.5, 5.0)
LINCS_DOSES = (0.0412, 0.123, 0.37, 1.11, 3.33, 10.0)        # dose_uM labels used by morph_lincs.py
MATCHED_T1 = [(0.0412, 0.05), (0.37, 0.5), (3.33, 5.0)]       # the three LINCS doses within 1.5x of Tahoe
TIERS = {
    "T1_LINCS_A549": dict(primary=(3.33, 5.0, "A549"), sig_line="A549"),
    "T2_EUOS_HepG2": dict(primary=(10.0, 5.0, "HepG2C3A"), sig_line="HepG2C3A"),
    "T3_JUMP_U2OS": dict(primary=(10.0, 5.0, "consensus"), sig_line="consensus"),
}
META = {"unit", "compound_id", "dose", "dose_uM", "n_wells", "n_plates", "sources", "plates", "sites", "cell_line",
        "time_h", "n_lines", "n_tahoe_drugs", "tahoe_plates", "magnitude", "tahoe_drugs", "plate", "n_cells_trt",
        "n_cells_ctrl", "cell_name", "cellosaurus", "n_cells", "jcp", "source", "batch", "well", "is_dmso", "site"}


# ================================================================ helpers
def die(msg):
    raise RuntimeError(msg)


def near(a, b, tol=0.02):
    return abs(np.log(a) - np.log(b)) < tol


def feat_cols(df):
    return [c for c in df.columns if c not in META and pd.api.types.is_numeric_dtype(df[c])
            and not c.startswith("Metadata")]


def as_matrix(df, units):
    d = df.set_index("unit").loc[units]
    X = d[feat_cols(df)].to_numpy(np.float64)
    X = X[:, ~np.all(np.isnan(X), axis=0)]
    return np.nan_to_num(X)


def cos(X):
    n = np.linalg.norm(X, axis=1, keepdims=True)
    n[n == 0] = 1
    Y = X / n
    return Y @ Y.T


def splitset(s):
    return set(str(s).split(";")) if isinstance(s, str) and s and s != "nan" else set()


def tri(n):
    return np.triu_indices(n, 1)


def rank_matrix(A):
    """Symmetric matrix of the ranks of A's upper-triangle entries (diag 0). Ranks are permutation-
    invariant, so permuting rows/cols of this matrix equals re-ranking the permuted A."""
    n = A.shape[0]
    iu = tri(n)
    R = np.zeros_like(A, dtype=np.float64)
    R[iu] = rankdata(A[iu])
    return R + R.T


def mantel(A, B, perms, strata=None, mask=None, keep_null=False):
    """Spearman Mantel r(A, B) over the upper triangle (optional fixed pair mask). Null permutes rows+cols
    of A freely or within strata (singletons fixed). One-sided p (null >= obs)."""
    n = A.shape[0]
    iu = tri(n)
    keep = np.ones(len(iu[0]), bool) if mask is None else mask[iu]
    groups = None
    if strata is not None:
        groups = [np.where(strata == s)[0] for s in np.unique(strata)]
        groups = [g for g in groups if len(g) > 1]
    if mask is None:
        RA, rb = rank_matrix(A), rankdata(B[iu])
        rb = (rb - rb.mean()) / rb.std()
        ra = RA[iu]
        r_obs = float(np.mean((ra - ra.mean()) / ra.std() * rb))
    else:
        rb = rankdata(B[iu][keep])
        r_obs = float(np.corrcoef(rankdata(A[iu][keep]), rb)[0, 1])
    null = np.empty(perms)
    for k in range(perms):
        p = np.arange(n)
        if groups is None:
            p = RNG.permutation(n)
        else:
            for g in groups:
                p[g] = RNG.permutation(g)
        if mask is None:
            ra = RA[np.ix_(p, p)][iu]
            null[k] = np.mean((ra - ra.mean()) / ra.std() * rb)
        else:
            null[k] = np.corrcoef(rankdata(A[np.ix_(p, p)][iu][keep]), rb)[0, 1]
    pval = (np.sum(null >= r_obs) + 1) / (perms + 1)
    n_perm_units = int(sum(len(g) for g in groups)) if groups is not None else n
    return dict(r=r_obs, p=pval, n_pairs=int(keep.sum()), n_permutable=n_perm_units,
                null_mean=float(null.mean()), null_sd=float(null.std()),
                null_q95=float(np.quantile(null, .95)), null=null if keep_null else None)


def boot_ci(A, B, n_boot, mask=None):
    """Unit bootstrap 95% CI for the Spearman Mantel r; pairs of a resampled unit with its own copy dropped."""
    n = A.shape[0]
    rs = []
    for _ in range(n_boot):
        idx = RNG.integers(0, n, n)
        i, j = np.triu_indices(n, 1)
        ok = idx[i] != idx[j]
        if mask is not None:
            ok &= mask[idx[i], idx[j]]
        a, b = A[idx[i][ok], idx[j][ok]], B[idx[i][ok], idx[j][ok]]
        if len(a) > 10:
            rs.append(spearmanr(a, b).statistic)
    return float(np.quantile(rs, .025)), float(np.quantile(rs, .975))


def nn_retrieval(S, labels, perms=1000, min_class=3):
    """Leave-one-compound-out MoA retrieval: each query (class size >= min_class) is removed from the database,
    the rest are ranked by similarity. Returns top-1 accuracy, mAP, n_queries, and label-permutation p-values."""
    labels = np.asarray(labels)

    def score(lab):
        vc = pd.Series(lab).value_counts()
        top1, aps = [], []
        for i in range(len(lab)):
            if vc[lab[i]] < min_class:
                continue
            o = np.argsort(-S[i], kind="stable")
            o = o[o != i]
            rel = (lab[o] == lab[i]).astype(float)
            top1.append(rel[0])
            aps.append((np.cumsum(rel) / np.arange(1, len(rel) + 1) * rel).sum() / rel.sum())
        return (np.mean(top1) if top1 else np.nan), (np.mean(aps) if aps else np.nan), len(aps)

    t1, m, nq = score(labels)
    null = np.array([score(RNG.permutation(labels))[:2] for _ in range(perms)])
    return dict(top1=t1, mAP=m, n_queries=nq, top1_null=float(null[:, 0].mean()), mAP_null=float(null[:, 1].mean()),
                p_top1=float((np.sum(null[:, 0] >= t1) + 1) / (perms + 1)),
                p_mAP=float((np.sum(null[:, 1] >= m) + 1) / (perms + 1)))


def testable(moa, min_classes=3, min_size=3):
    vc = pd.Series(moa[moa != "unclear"]).value_counts()
    return int((vc >= min_size).sum()) >= min_classes, int((vc >= min_size).sum())


# ================================================================ loading
def load_units():
    u = pd.read_csv(S2 / "units.csv")
    return u


def gate_status(tier):
    g = S2 / f"gate_{tier.split('_')[0].lower()}.json"
    if g.exists():
        v = json.loads(g.read_text()).get("verdict", "")
        if v in ("UNDERPOWERED", "UNAVAILABLE"):
            return v
    return "OK"


class Tahoe:
    def __init__(self):
        pl = S2 / "tahoe_perline.parquet"
        if not pl.exists():
            die("missing step2/tahoe_perline.parquet")
        self.pl = pd.read_parquet(pl)
        if "unit" not in self.pl:
            die("tahoe_perline.parquet has no 'unit' column (old schema)")
        # identify the Cellosaurus-id column
        cid = [c for c in self.pl.columns if self.pl[c].dtype == object and self.pl[c].astype(str).str.startswith("CVCL").mean() > .9]
        if not cid:
            die("no Cellosaurus column in tahoe_perline")
        self.cid = cid[0]
        self.genes = feat_cols(self.pl)
        self.pcol = "plates" if "plates" in self.pl else ("plate" if "plate" in self.pl else None)
        self.pl["dose_uM"] = self.pl["dose_uM"].astype(float)
        self.cons = pd.read_parquet(ROOT / "tahoe_signatures.parquet")
        if "unit" not in self.cons:
            die("tahoe_signatures.parquet has no 'unit' column (old schema)")
        self.cons["dose_uM"] = self.cons["dose_uM"].astype(float)
        self._cent = None

    def lines(self):
        return sorted(self.pl[self.cid].unique())

    def sig(self, which, dose, centered=False):
        """-> DataFrame with unit, plates, genes. which in {'A549','HepG2C3A','consensus', <CVCL id>}."""
        pl = self.centered() if centered else self.pl
        if which == "consensus":
            if not centered:
                d = self.cons[[near(x, dose) for x in self.cons.dose_uM]].copy()
                pc = "plates" if "plates" in d else "tahoe_plates"
                d["plates"] = d[pc].astype(str)
                return d
            d = pl[[near(x, dose) for x in pl.dose_uM]]
            g = d.groupby("unit")[self.genes].mean().reset_index()
            g["plates"] = d.groupby("unit")[self.pcol].apply(lambda s: ";".join(sorted({p for x in s for p in str(x).split(";")}))).values
            return g
        cid = {"A549": A549, "HepG2C3A": HEPG2}.get(which, which)
        d = pl[(pl[self.cid] == cid) & [near(x, dose) for x in pl.dose_uM]].copy()
        d["plates"] = d[self.pcol].astype(str)
        return d

    def centered(self):
        """Per (line, dose, plate) centering: subtract the mean signature of all units on that plate -
        removes the shared-DMSO-control component that every drug on a plate has in common."""
        # Uses the ALL-drug plate x dose x line mean (tahoe_plate_line_mean.parquet, ~80 drugs per plate),
        # averaged over the unit's own plates. (Grouping on the joined plate string would leave multi-plate
        # units in groups of 1-2 and zero them out.)
        if self._cent is None:
            f = S2 / "tahoe_plate_line_mean.parquet"
            if not f.exists():
                die("missing step2/tahoe_plate_line_mean.parquet (needed for plate centring)")
            pm = pd.read_parquet(f)
            pm["dose_uM"] = pm.dose_uM.astype(float).round(4)
            G = pm[self.genes].to_numpy(np.float32)
            idx = {(str(p), d, c): i for i, (p, d, c) in enumerate(zip(pm.plate.astype(str), pm.dose_uM, pm.cell_line))}
            d = self.pl.copy()
            X = d[self.genes].to_numpy(np.float32)
            miss = 0
            for r, (pl, ds, c) in enumerate(zip(d[self.pcol].astype(str), d.dose_uM.round(4), d[self.cid])):
                ii = [idx[(p, ds, c)] for p in pl.split(";") if (p, ds, c) in idx]
                if ii:
                    X[r] -= G[ii].mean(0)
                else:
                    miss += 1
            if miss:
                die(f"plate centring: {miss} rows without a plate x line mean")
            d[self.genes] = X
            self._cent = d
        return self._cent


def load_morph(tier, units):
    """-> (profiles DataFrame keyed by unit + dose_uM, activity DataFrame or None, wells DataFrame or None)."""
    if tier == "T3_JUMP_U2OS":
        u3 = units[units.tier == tier]
        m = {}
        for r in u3.itertuples():
            for j in str(r.dataset_ids).split(";"):
                if j in m and m[j] != r.unit:
                    die(f"T3 JCP {j} maps to two units")
                m[j] = r.unit
        if len(set(m.values())) != len(m):
            die("T3 JCP->unit mapping is not 1:1")
        prof, wells = [], []
        for dose, uM in (("high", 10.0), ("low", 0.625)):
            p = pd.read_parquet(ROOT / f"morph_profiles_{dose}.parquet")
            p["unit"], p["dose_uM"] = p.compound_id.map(m), uM
            if p.unit.isna().any():
                die("unmapped T3 compound")
            prof.append(p)
            w = pd.read_parquet(S2 / f"morph_wells_{dose}.parquet")
            w["unit"] = np.where(w.is_dmso, "DMSO", w.compound_id.map(m))
            w["dose_uM"] = uM
            wells.append(w)
        a = pd.read_csv(S2 / "morph_activity.csv")
        a["unit"] = a.compound_id.map(m)
        a["dose_uM"] = a.dose.map({"high": 10.0, "low": 0.625})
        a["active_main"] = a["active_strict"].astype(bool)
        return pd.concat(prof), a, pd.concat(wells)
    stem = {"T1_LINCS_A549": "lincs", "T2_EUOS_HepG2": "euos"}[tier]
    f = S2 / f"morph_{stem}.parquet"
    if not f.exists():
        die(f"missing {f}")
    p = pd.read_parquet(f)
    if "cell_line" in p and stem == "euos":
        p = p[p.cell_line.astype(str).str.lower().str.contains("hepg2")]
    a = pd.read_csv(S2 / f"morph_{stem}_activity.csv") if (S2 / f"morph_{stem}_activity.csv").exists() else None
    if a is not None:
        # LINCS: replicate positions are fixed across plates, so the DMSO-based calls are position-driven;
        # use the position-matched replicate call (morph_lincs.py, lincs_notes.md). EU-OS: DMSO calls saturate,
        # use the compound-null call. Activity only feeds the T1 lowest-active-dose analysis.
        col = next(c for c in ("active_rep_pos", "active_cpd_and_pos", "active_vs_cpd", "active_z") if c in a)
        a["active_main"] = a[col].fillna(False).astype(bool)
        a["active_main_na"] = a[col].isna()
        a.attrs["active_col"] = col
    wf = S2 / f"morph_{stem}_wells.parquet"
    w = pd.read_parquet(wf) if wf.exists() else None
    if w is not None and "dose_uM" not in w:           # EU-OS wells: single 10 uM compound dose
        w = w.copy()
        w["dose_uM"] = 10.0
    return p, a, w


# ================================================================ split-half reliability & replicate reproducibility
def well_unit_col(w):
    return "unit" if "unit" in w else "compound_id"


def split_half(w, units, dose, n_split=20):
    """Mantel r between similarity matrices built from two random halves of each unit's replicate wells,
    Spearman-Brown corrected. Only units with >= 2 wells."""
    if w is None:
        return dict(r_half=np.nan, r_sb=np.nan, n=0)
    uc = well_unit_col(w)
    d = w[(w[uc].isin(units)) & [near(x, dose) for x in w.dose_uM]]
    fc = [c for c in feat_cols(d) if d[c].notna().any()]
    groups = {u: g[fc].to_numpy(np.float64) for u, g in d.groupby(uc) if len(g) >= 2}
    us = sorted(groups)
    if len(us) < 10:
        return dict(r_half=np.nan, r_sb=np.nan, n=len(us))
    rs = []
    for _ in range(n_split):
        H1, H2 = [], []
        for u in us:
            X = groups[u][RNG.permutation(len(groups[u]))]
            k = len(X) // 2
            H1.append(np.median(X[:k], 0)), H2.append(np.median(X[k:], 0))
        A, B = cos(np.nan_to_num(np.array(H1))), cos(np.nan_to_num(np.array(H2)))
        iu = tri(len(us))
        rs.append(spearmanr(A[iu], B[iu]).statistic)
    r = float(np.mean(rs))
    return dict(r_half=r, r_sb=2 * r / (1 + r) if r > -1 else np.nan, n=len(us))


def replicate_repro(w, units, dose, n_null=500):
    """Uniform replicate reproducibility for every dataset: replicate mAP among the compound wells at this dose
    (DMSO excluded; pairs on the same plate excluded), null = random well sets of the same size.
    Returns fraction of units with p < 0.05 and the median mAP."""
    if w is None:
        return dict(frac_reproducible=np.nan, median_mAP=np.nan, n=0)
    uc = well_unit_col(w)
    d = w[(w[uc].isin(units)) & [near(x, dose) for x in w.dose_uM]].reset_index(drop=True)
    fc = [c for c in feat_cols(d) if d[c].notna().any()]
    X = np.nan_to_num(d[fc].to_numpy(np.float64))
    S = cos(X)
    lab = d[uc].to_numpy()
    plate = d["plate"].astype(str).to_numpy() if "plate" in d else np.arange(len(d)).astype(str)
    np.fill_diagonal(S, -np.inf)

    def ap_set(idx):
        aps = []
        idxs = set(idx)
        for i in idx:
            ok = plate != plate[i]
            ok[i] = False
            o = np.argsort(-S[i][ok])
            rel = np.isin(np.where(ok)[0][o], list(idxs - {i})).astype(float)
            if rel.sum() == 0:
                continue
            aps.append((np.cumsum(rel) / np.arange(1, len(rel) + 1) * rel).sum() / rel.sum())
        return np.mean(aps) if aps else np.nan

    res, nulls = [], {}
    for u, g in d.groupby(uc):
        idx = g.index.to_numpy()
        if len(idx) < 2:
            continue
        k = len(idx)
        if k not in nulls:
            nulls[k] = np.array([ap_set(RNG.choice(len(d), k, replace=False)) for _ in range(n_null)])
        m = ap_set(idx)
        nl = nulls[k][~np.isnan(nulls[k])]
        res.append((u, m, (np.sum(nl >= m) + 1) / (len(nl) + 1) if not np.isnan(m) else np.nan))
    r = pd.DataFrame(res, columns=["unit", "mAP", "p"])
    return dict(frac_reproducible=float((r.p < 0.05).mean()), median_mAP=float(r.mAP.median()), n=len(r))


# ================================================================ core analysis
def compare(tier, mdose, tdose, line, unit_list, morph, tahoe, units_tab, args, label, rows, nulls,
            secondary=True, centered=False, mask_plates=True, retrieval=True, full=False):
    m = morph[[near(x, mdose) for x in morph.dose_uM]]
    t = tahoe.sig(line, tdose, centered=centered)
    us = [u for u in unit_list if u in set(m.unit) and u in set(t.unit)]
    if len(us) < 10:
        rows.append(dict(tier=tier, label=label, morph_dose=mdose, tahoe_dose=tdose, tahoe_sig=line,
                         centered=centered, n_units=len(us), note="too few units"))
        return None
    A, B = cos(as_matrix(m, us)), cos(as_matrix(t, us))
    moa = units_tab.drop_duplicates("unit").set_index("unit").loc[us, "moa"].to_numpy()
    ann = np.where(moa != "unclear")[0]
    base = dict(tier=tier, label=label, morph_dose=mdose, tahoe_dose=tdose, tahoe_sig=line, centered=centered)
    ok, ncls = testable(moa)
    Aa, Ba = A[np.ix_(ann, ann)], B[np.ix_(ann, ann)]
    perms = args.perms if label == "PRIMARY" else args.perms2
    if ok:
        res = mantel(Aa, Ba, perms, strata=moa[ann], keep_null=True)
        lo, hi = boot_ci(Aa, Ba, args.boot)
        rows.append({**base, "analysis": "annotated_within_moa", "n_units": len(ann), "n_classes_ge3": ncls,
                     **{k: v for k, v in res.items() if k != "null"}, "ci_low": lo, "ci_high": hi})
        nulls[f"{tier}|{label}|{mdose}|{tdose}|{line}|{centered}|within"] = res["null"]
    else:
        rows.append({**base, "analysis": "annotated_within_moa", "n_units": len(ann), "n_classes_ge3": ncls,
                     "note": "not testable (<3 MoA classes with >=3 units)"})
    res = mantel(Aa, Ba, args.perms2)
    rows.append({**base, "analysis": "annotated_permutation", "n_units": len(ann), **{k: v for k, v in res.items() if k != "null"}})
    res = mantel(A, B, args.perms2, keep_null=True)
    lo, hi = boot_ci(A, B, args.boot)
    rows.append({**base, "analysis": "all_permutation", "n_units": len(us), **{k: v for k, v in res.items() if k != "null"},
                 "ci_low": lo, "ci_high": hi})
    nulls[f"{tier}|{label}|{mdose}|{tdose}|{line}|{centered}|all"] = res["null"]
    out = dict(us=us, A=A, B=B, moa=moa)
    if not secondary:
        return out
    if mask_plates:
        jp = [splitset(x) for x in m.set_index("unit").loc[us, "plates"]] if "plates" in m else [set()] * len(us)
        tp = [splitset(x) for x in t.set_index("unit").loc[us, "plates"]]
        n = len(us)
        same_t = np.array([[bool(tp[i] & tp[j]) for j in range(n)] for i in range(n)])
        same_m = np.array([[bool(jp[i] & jp[j]) for j in range(n)] for i in range(n)])
        iu = tri(n)
        rows.append({**base, "analysis": "tahoe_same_vs_diff_plate_similarity", "n_units": n,
                     "r": float(B[iu][same_t[iu]].mean() - B[iu][~same_t[iu]].mean()),
                     "note": f"mean cos same-plate {B[iu][same_t[iu]].mean():.3f} vs diff {B[iu][~same_t[iu]].mean():.3f}; "
                             f"frac pairs same tahoe plate {same_t[iu].mean():.3f}, same morph plate {same_m[iu].mean():.3f}"})
        res = mantel(A, B, args.perms2, mask=~(same_t | same_m))
        rows.append({**base, "analysis": "all_no_shared_plates", "n_units": n, **{k: v for k, v in res.items() if k != "null"}})
    if retrieval:
        for mod, S in (("morph", Aa), ("tahoe", Ba)):
            rr = nn_retrieval(S, moa[ann], perms=min(args.perms2, 1000))
            rows.append({**base, "analysis": f"moa_retrieval_{mod}", "n_units": len(ann), **rr})
    return out


def rowwise(A, B):
    n = A.shape[0]
    return np.array([spearmanr(A[i, np.arange(n) != i], B[i, np.arange(n) != i]).statistic for i in range(n)])


def consistency_test(tier, out, tdose, args, rows, scatter, ccol="consistency"):
    cons = pd.read_csv(S2 / "tahoe_consistency.csv")
    cons = cons[[near(x, tdose) for x in cons.dose_uM]].set_index("unit")
    us = out["us"]
    agree = rowwise(out["A"], out["B"])
    cc = cons.reindex(us)[ccol].to_numpy(float)
    mag = cons.reindex(us)["magnitude" if ccol == "consistency" else "magnitude_centred"].to_numpy(float)
    sfx = "" if ccol == "consistency" else "__" + ccol
    ok = ~np.isnan(cc) & ~np.isnan(agree)
    r = spearmanr(cc[ok], agree[ok]).statistic
    bs = []
    idx_all = np.where(ok)[0]
    for _ in range(args.boot):
        b = RNG.choice(idx_all, len(idx_all))
        bs.append(spearmanr(cc[b], agree[b]).statistic)
    # partial (rank-residualised on magnitude)
    rx, ry, rz = rankdata(cc[ok]), rankdata(agree[ok]), rankdata(mag[ok])
    Z = np.column_stack([np.ones_like(rz), rz])
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    pr = np.corrcoef(ex, ey)[0, 1]
    pbs = []
    for _ in range(args.boot):
        b = RNG.choice(len(ex), len(ex))
        pbs.append(np.corrcoef(ex[b], ey[b])[0, 1])
    rows.append(dict(tier=tier, label="SECONDARY_consistency", tahoe_dose=tdose, analysis="consistency_vs_rowwise_agreement" + sfx,
                     n_units=int(ok.sum()), r=r, ci_low=np.quantile(bs, .025), ci_high=np.quantile(bs, .975),
                     p=spearmanr(cc[ok], agree[ok]).pvalue))
    rows.append(dict(tier=tier, label="SECONDARY_consistency", tahoe_dose=tdose,
                     analysis="consistency_vs_agreement_partial_magnitude" + sfx, n_units=int(ok.sum()), r=pr,
                     ci_low=np.quantile(pbs, .025), ci_high=np.quantile(pbs, .975)))
    rows.append(dict(tier=tier, label="SECONDARY_consistency", tahoe_dose=tdose, analysis="consistency_vs_magnitude" + sfx,
                     n_units=int(ok.sum()), r=spearmanr(cc[ok], mag[ok]).statistic))
    if sfx:
        return
    scatter.append(pd.DataFrame(dict(tier=tier, unit=np.array(us)[ok], consistency=cc[ok], magnitude=mag[ok],
                                     rowwise_agreement=agree[ok], moa=out["moa"][ok])))


# ================================================================ T1 dose response & lowest active dose
def lowest_active(units_t1, morph_act, tahoe, rows, args, subset=""):
    ta = pd.read_csv(S2 / "tahoe_activity.csv")
    scol = next(c for c in ta.columns if c in ("signature", "which", "sig", "cell_line", "line"))
    ta = ta[ta[scol].astype(str).str.contains("A549|CVCL_0023", regex=True)]
    zcol = next(c for c in ta.columns if c in ("z_norm", "z", "norm_z", "robust_z"))
    mz = "norm_z" if "norm_z" in morph_act else None
    res = []
    for u in units_t1:
        rec = {"unit": u}
        for k, (ld, td) in enumerate(MATCHED_T1):
            ma = morph_act[(morph_act.unit == u) & [near(x, ld) for x in morph_act.dose_uM]]
            tt = ta[(ta.unit == u) & [near(x, td) for x in ta.dose_uM]]
            rec[f"m_act_{k}"] = bool(ma.active_main.any()) if len(ma) else np.nan
            rec[f"t_act_{k}"] = bool((tt[zcol] > 3).any()) if len(tt) else np.nan
            rec[f"m_z_{k}"] = float(ma[mz].max()) if (mz and len(ma)) else np.nan
            rec[f"t_z_{k}"] = float(tt[zcol].max()) if len(tt) else np.nan
        res.append(rec)
    d = pd.DataFrame(res)

    def first(row, pre):
        for k in range(3):
            if row[f"{pre}_{k}"] is True or row[f"{pre}_{k}"] == 1:
                return k
        return np.nan

    def sustained(row, pre):
        """lowest matched level from which every higher matched level is also active"""
        best = np.nan
        for k in (2, 1, 0):
            v = row[f"{pre}_{k}"]
            if v is True or v == 1:
                best = k
            else:
                break
        return best

    def classify(d, mcol, tcol, tag, fn=None):
        fn = fn or first
        tag = f"{tag}{subset}"
        lm = d.apply(lambda r: fn(r, mcol), axis=1)
        lt = d.apply(lambda r: fn(r, tcol), axis=1)
        cls = np.select([lm.isna() & lt.isna(), lt.isna() | (lm < lt), lm.isna() | (lt < lm)],
                        ["neither", "morph_first", "rna_first"], "same")
        vc = pd.Series(cls).value_counts()
        both = (~lm.isna()) & (~lt.isna())
        from scipy.stats import binomtest
        mf, rf = int(vc.get("morph_first", 0)), int(vc.get("rna_first", 0))
        bt = binomtest(mf, mf + rf) if mf + rf else None
        rows.append(dict(tier="T1_LINCS_A549", label="DOSE_lowest_active", analysis=f"lowest_active_dose_{tag}",
                         n_units=len(d), morph_first=mf, rna_first=rf, same=int(vc.get("same", 0)),
                         neither=int(vc.get("neither", 0)), p=bt.pvalue if bt else np.nan,
                         note=f"active fraction per matched level morph={[round(float(np.nanmean(d[f'{mcol}_{k}'].astype(float))),3) for k in range(3)]} "
                              f"tahoe={[round(float(np.nanmean(d[f'{tcol}_{k}'].astype(float))),3) for k in range(3)]}; "
                              f"lowest-dose shift among units active in both: median morph level {lm[both].median()} vs tahoe {lt[both].median()}"))
        return cls

    d["class_specified"] = classify(d, "m_act", "t_act", "specified_thresholds")
    d["class_sustained"] = classify(d, "m_act", "t_act", "sustained_specified_thresholds", fn=sustained)
    # call-rate-matched thresholds: Tahoe z threshold set so its active fraction at the top matched dose
    # equals morphology's there (or vice versa), then reclassify on the z scores.
    if not d["m_z_2"].isna().all():
        fm = np.nanmean(d["m_act_2"].astype(float))
        thr_t = np.nanquantile(d["t_z_2"], 1 - fm) if 0 < fm < 1 else 3
        for k in range(3):
            d[f"t_act2_{k}"] = d[f"t_z_{k}"] > thr_t
            d[f"m_act2_{k}"] = d[f"m_act_{k}"]
        d["class_matched"] = classify(d, "m_act2", "t_act2", f"callrate_matched_tahoe_z>{thr_t:.2f}")
        d["class_matched_sustained"] = classify(d, "m_act2", "t_act2", f"sustained_callrate_matched", fn=sustained)
    # strictness-matched: Tahoe z residualised on log(n_cells_trt) within dose (the norm tracks cell count),
    # active if residual robust z > 3; morphology threshold on the position-matched replicate p chosen so that
    # its active rate at the top matched dose equals Tahoe's there.
    a549 = tahoe.sig("A549", 5.0).iloc[:0]
    pl = tahoe.pl[tahoe.pl[tahoe.cid] == A549][["unit", "dose_uM", "n_cells_trt"]]
    tt = ta.copy() if "n_cells_trt" in ta else ta.merge(pl, on=["unit", "dose_uM"], how="left")
    tt = tt.reset_index(drop=True)
    tt["res"] = np.nan
    for dd, g in tt.groupby(tt.dose_uM.round(3)):
        ok = g.n_cells_trt.notna() & g[zcol].notna()
        X = np.column_stack([np.ones(ok.sum()), np.log(g.loc[ok, "n_cells_trt"].clip(lower=1))])
        b = np.linalg.lstsq(X, g.loc[ok, zcol].to_numpy(float), rcond=None)[0]
        e = g.loc[ok, zcol].to_numpy(float) - X @ b
        tt.loc[g.index[ok], "res"] = (e - np.median(e)) / (1.4826 * np.median(np.abs(e - np.median(e))))
    pcol = "rep_pos_p" if "rep_pos_p" in morph_act else None
    if pcol:
        for k, (ld, td) in enumerate(MATCHED_T1):
            t3 = tt[[near(x, td) for x in tt.dose_uM]].set_index("unit")["res"]
            d[f"t_act3_{k}"] = d.unit.map(t3) > 3
            ma = morph_act[[near(x, ld) for x in morph_act.dose_uM]].set_index("unit")
            # composite score: -log10(position-matched p) with norm_z as tie-breaker (p ties at its minimum)
            sc = -np.log10(ma[pcol].clip(lower=1e-6)) + 1e-3 * ma.get("norm_z", 0)
            d[f"m_p_{k}"] = d.unit.map(sc)
        rate_t = float(np.nanmean(d["t_act3_2"].astype(float)))
        thr_m = float(np.nanquantile(d["m_p_2"], 1 - rate_t)) if rate_t > 0 else np.inf
        for k in range(3):
            d[f"m_act3_{k}"] = d[f"m_p_{k}"] > thr_m
        d["class_strict_matched"] = classify(d, "m_act3", "t_act3", f"strict_matched_tahoe_resid_z>3_morph_score>{thr_m:.3f}")
        d["class_strict_matched_sustained"] = classify(d, "m_act3", "t_act3", "sustained_strict_matched", fn=sustained)
        both = sum(any(d.loc[i, f"m_act3_{k}"] is True or d.loc[i, f"m_act3_{k}"] == 1 for k in range(3)) and
                   any(d.loc[i, f"t_act3_{k}"] is True or d.loc[i, f"t_act3_{k}"] == 1 for k in range(3)) for i in d.index)
        rows.append(dict(tier="T1_LINCS_A549", label="DOSE_lowest_active", analysis=f"strict_matched_units_active_in_both{subset}",
                         n_units=len(d), r=both, note=f"tahoe top-dose rate {rate_t:.3f}; morph score threshold {thr_m:.4g}"))
    d.to_csv(OUT / f"T1_lowest_active_per_unit{subset}.csv", index=False)
    # n_cells check: does Tahoe norm track n_cells_trt? (fewer cells -> noisier, larger norms)
    a = tahoe.sig("A549", 5.0)
    if "n_cells_trt" in a:
        nrm = np.linalg.norm(np.nan_to_num(a[tahoe.genes].to_numpy(float)), axis=1)
        rows.append(dict(tier="T1_LINCS_A549", label="DOSE_lowest_active", analysis="tahoe_A549_norm_vs_ncells",
                         n_units=len(a), r=spearmanr(nrm, a.n_cells_trt).statistic, p=spearmanr(nrm, a.n_cells_trt).pvalue))


# ================================================================ Tahoe noise ceiling (split-plate reliability)
def tahoe_reliability(rows, units, morphs, tahoe, args):
    """Reliability of the Tahoe similarity matrix: units observed on >= 2 plates at a dose; build the matrix from
    plate copy 1 and from plate copy 2 (independent wells and DMSO controls) and Spearman-Mantel them. Done for the
    cross-line consensus and for A549 / HepG2C3A single lines, raw and plate-centred. Then the ceiling on each tier's
    primary r = sqrt(rel_morph_SB * rel_tahoe), and the primary rerun on units whose signature replicates."""
    fc = S2 / "tahoe_perplate_consensus.parquet"
    fl = S2 / "tahoe_perplate.parquet"
    if not (fc.exists() and fl.exists()):
        die("missing step2/tahoe_perplate*.parquet (needed for the Tahoe noise ceiling)")
    pc, pl = pd.read_parquet(fc), pd.read_parquet(fl)
    pm = pd.read_parquet(S2 / "tahoe_plate_line_mean.parquet")
    genes = tahoe.genes
    cid = [c for c in pl.columns if pl[c].dtype == object and pl[c].astype(str).str.startswith("CVCL").mean() > .9][0]
    # plate-centring references: per plate x dose x line, and the line-mean of those for the consensus
    pm["dose_uM"] = pm.dose_uM.astype(float).round(4)
    pmc = pm.groupby(["plate", "dose_uM"])[genes].mean()
    rel, rep_units = {}, set()
    for dose in TAHOE_DOSES:
        for which in ("consensus", "A549", "HepG2C3A"):
            if which == "consensus":
                d = pc[[near(x, dose) for x in pc.dose_uM]].copy()
            else:
                d = pl[(pl[cid] == {"A549": A549, "HepG2C3A": HEPG2}[which]) & np.array([near(x, dose) for x in pl.dose_uM])].copy()
            d["plate"] = d.plate.astype(str)
            for centred in (False, True):
                X = d[genes].to_numpy(np.float64)
                if centred:
                    if which == "consensus":
                        ref = pmc.reindex(list(zip(d.plate, d.dose_uM.astype(float).round(4))))
                    else:
                        ref = pm.set_index(["plate", "dose_uM", "cell_line"])[genes].reindex(
                            list(zip(d.plate, d.dose_uM.astype(float).round(4), d[cid])))
                    X = X - np.nan_to_num(ref.to_numpy(np.float64))
                d2 = d[["unit", "plate"]].copy().reset_index(drop=True)
                g = d2.groupby("unit").indices
                us = sorted(u for u, ix in g.items() if len(ix) >= 2)
                if len(us) < 10:
                    rows.append(dict(tier="TAHOE", label="TAHOE_RELIABILITY", analysis=f"split_plate_{which}",
                                     tahoe_dose=dose, centered=centred, n_units=len(us), note="too few multi-plate units"))
                    continue
                rs, unit_r = [], {u: [] for u in us}
                for _ in range(50):
                    i1, i2 = [], []
                    for u in us:
                        ix = RNG.permutation(g[u])[:2]
                        i1.append(ix[0]), i2.append(ix[1])
                    A, B = cos(X[i1]), cos(X[i2])
                    iu = tri(len(us))
                    rs.append(spearmanr(A[iu], B[iu]).statistic)
                    for k, u in enumerate(us):
                        unit_r[u].append(np.corrcoef(X[i1[k]], X[i2[k]])[0, 1])
                r = float(np.mean(rs))
                rel[(dose, which, centred)] = r
                rows.append(dict(tier="TAHOE", label="TAHOE_RELIABILITY", analysis=f"split_plate_{which}", tahoe_dose=dose,
                                 centered=centred, n_units=len(us), r=r, ci_low=float(np.quantile(rs, .025)),
                                 ci_high=float(np.quantile(rs, .975)),
                                 note=f"median per-unit cross-plate r {np.median([np.mean(v) for v in unit_r.values()]):.3f}"))
                if dose == 5.0 and which == "consensus" and centred:
                    rep_units = {u for u, v in unit_r.items() if np.mean(v) > 0.2}
    # ceilings and replicating-subset primaries
    res = pd.DataFrame(rows)
    for tier, cfg in TIERS.items():
        if tier not in morphs:
            continue
        mdose, tdose, line = cfg["primary"]
        msb = res[(res.tier == tier) & (res.analysis == "morph_split_half") & np.isclose(res.morph_dose, mdose)]["r_sb"]
        for centred in (False, True):
            tr = rel.get((tdose, line, centred))
            if tr is None or msb.empty:
                continue
            ceil = float(np.sqrt(max(msb.iloc[0], 0) * max(tr, 0)))
            rows.append(dict(tier=tier, label="CEILING", analysis="max_attainable_mantel_r", morph_dose=mdose, tahoe_dose=tdose,
                             tahoe_sig=line, centered=centred, r=ceil,
                             note=f"sqrt(morph split-half SB {msb.iloc[0]:.3f} x Tahoe single-plate reliability {tr:.3f}); "
                                  f"Tahoe reliability is for single-plate signatures (most units have one plate)"))
        ut = set(units.loc[units.tier == tier, "unit"])
        sub = sorted(ut & rep_units)
        rows.append(dict(tier=tier, label="REPLICATING", analysis="n_units_replicating_across_plates", n_units=len(sub),
                         note="centred consensus cross-plate r > 0.2 at 5 uM"))
        if len(sub) >= 10:
            compare(tier, mdose, tdose, line, sub, morphs[tier][0], tahoe, units, args, "SENS_replicating_units", rows, {},
                    secondary=False, centered=True)
        # Like-for-like ceiling: multi-plate units (where Tahoe reliability is measurable) vs the single-plate rest.
        multi_all = set(pc.loc[[near(x, tdose) for x in pc.dose_uM]].groupby("unit").plate.nunique().loc[lambda s: s >= 2].index)
        multi, single = sorted(ut & multi_all), sorted(ut - multi_all)
        rows.append(dict(tier=tier, label="MULTIPLATE", analysis="n_units_multi_vs_single_plate", n_units=len(multi),
                         note=f"multi-plate {len(multi)}, single-plate {len(single)}; replicating subset {len(sub)} "
                              f"({len(set(sub) & set(multi))} of them multi-plate)"))
        for lab, subset in (("SENS_multiplate_units", multi), ("SENS_singleplate_units", single)):
            if len(subset) >= 10:
                compare(tier, mdose, tdose, line, subset, morphs[tier][0], tahoe, units, args, lab, rows, {},
                        secondary=False, centered=True)
        tr_m = split_plate_rel(pc, pl, pm, pmc, genes, cid, multi, line, tdose, centred=True)
        sh_m = split_half(morphs[tier][2], multi, mdose)
        if tr_m is not None and not np.isnan(sh_m["r_sb"]):
            rows.append(dict(tier=tier, label="CEILING", analysis="max_attainable_mantel_r_multiplate_units", morph_dose=mdose,
                             tahoe_dose=tdose, tahoe_sig=line, centered=True, n_units=len(multi),
                             r=float(np.sqrt(max(sh_m["r_sb"], 0) * max(tr_m, 0))),
                             note=f"same units: morph split-half SB {sh_m['r_sb']:.3f} x Tahoe split-plate {tr_m:.3f}"))
        # do multi-plate units differ from single-plate ones? (if so the ceiling does not transfer)
        cons = pd.read_csv(S2 / "tahoe_consistency.csv")
        cons = cons[[near(x, tdose) for x in cons.dose_uM]].set_index("unit")
        act = pd.read_csv(S2 / "tahoe_activity.csv")
        ctx = {"A549": "A549", "HepG2C3A": "HepG2C3A", "consensus": "consensus"}[line]
        act = act[(act.context == ctx) & np.array([near(x, tdose) for x in act.dose_uM])].set_index("unit")
        for nm, subset in (("multi", multi), ("single", single)):
            rows.append(dict(tier=tier, label="MULTIPLATE", analysis=f"tahoe_profile_{nm}_plate_units", n_units=len(subset),
                             r=float(cons.reindex(subset)["magnitude_centred"].median()),
                             note=f"median centred magnitude; active fraction {act.reindex(subset)['active'].astype(float).mean():.3f}; "
                                  f"median centred consistency {cons.reindex(subset)['consistency_centred'].median():.3f}"))


def split_plate_rel(pc, pl, pm, pmc, genes, cid, subset, which, dose, centred=True, n_rep=50):
    """Tahoe split-plate similarity-matrix reliability restricted to `subset` units (see tahoe_reliability)."""
    if which == "consensus":
        d = pc[[near(x, dose) for x in pc.dose_uM]].copy()
    else:
        d = pl[(pl[cid] == {"A549": A549, "HepG2C3A": HEPG2}[which]) & np.array([near(x, dose) for x in pl.dose_uM])].copy()
    d = d[d.unit.isin(subset)].reset_index(drop=True)
    d["plate"] = d.plate.astype(str)
    X = d[genes].to_numpy(np.float64)
    if centred:
        if which == "consensus":
            ref = pmc.reindex(list(zip(d.plate, d.dose_uM.astype(float).round(4))))
        else:
            ref = pm.set_index(["plate", "dose_uM", "cell_line"])[genes].reindex(
                list(zip(d.plate, d.dose_uM.astype(float).round(4), d[cid])))
        X = X - np.nan_to_num(ref.to_numpy(np.float64))
    g = d.groupby("unit").indices
    us = sorted(u for u, ix in g.items() if len(ix) >= 2)
    if len(us) < 10:
        return None
    rs = []
    for _ in range(n_rep):
        i1, i2 = zip(*[tuple(RNG.permutation(g[u])[:2]) for u in us])
        A, B = cos(X[list(i1)]), cos(X[list(i2)])
        iu = tri(len(us))
        rs.append(spearmanr(A[iu], B[iu]).statistic)
    return float(np.mean(rs))

# ================================================================ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=9999)
    ap.add_argument("--perms2", type=int, default=1999)
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--tiers", default="T1_LINCS_A549,T2_EUOS_HepG2,T3_JUMP_U2OS")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)
    units = load_units()
    tahoe = Tahoe()
    rows, nulls, scatter, curves = [], {}, [], []
    morphs = {}
    active_tiers = []
    for tier in args.tiers.split(","):
        g = gate_status(tier)
        if g != "OK":
            rows.append(dict(tier=tier, label="GATE", analysis="gate", note=g))
            print(f"{tier}: gated out ({g})")
            continue
        morphs[tier] = load_morph(tier, units)
        active_tiers.append(tier)

    for tier in active_tiers:
        prof, act, wells = morphs[tier]
        cfg = TIERS[tier]
        mdose, tdose, line = cfg["primary"]
        ut = units[units.tier == tier]
        all_u = ut.unit.tolist()
        print(f"== {tier}: {len(all_u)} units", flush=True)
        # PRIMARY
        out = compare(tier, mdose, tdose, line, all_u, prof, tahoe, units, args, "PRIMARY", rows, nulls)
        # 3x3: every Tahoe signature choice; plate-centred variant
        for sig in ("A549", "HepG2C3A", "consensus"):
            if sig != line:
                compare(tier, mdose, tdose, sig, all_u, prof, tahoe, units, args, "SIG_CHOICE", rows, nulls,
                        secondary=False)
            if sig != line:
                compare(tier, mdose, tdose, sig, all_u, prof, tahoe, units, args, "PLATE_CENTERED", rows, nulls,
                        secondary=False, centered=True)
        compare(tier, mdose, tdose, line, all_u, prof, tahoe, units, args, "PRIMARY_CENTERED", rows, nulls,
                secondary=True, centered=True)
        # rank of the matched line among all 50 (plain all-unit Spearman Mantel r, no permutations)
        m = prof[[near(x, mdose) for x in prof.dose_uM]]
        lr = []
        for cid in tahoe.lines():
            t = tahoe.sig(cid, tdose, centered=True)   # plate-centred: raw per-line matrices are plate-dominated
            us = [u for u in all_u if u in set(m.unit) and u in set(t.unit)]
            if len(us) >= 10:
                A, B = cos(as_matrix(m, us)), cos(as_matrix(t, us))
                iu = tri(len(us))
                lr.append((cid, spearmanr(A[iu], B[iu]).statistic, len(us)))
        lr = pd.DataFrame(lr, columns=["line", "r", "n"]).sort_values("r", ascending=False).reset_index(drop=True)
        lr.to_csv(OUT / f"{tier}_line_ranking.csv", index=False)
        target = {"A549": A549, "HepG2C3A": HEPG2}.get(line)
        for name, cid in (("A549", A549), ("HepG2C3A", HEPG2)):
            if cid in set(lr.line):
                rows.append(dict(tier=tier, label="LINE_RANK", analysis=f"rank_of_{name}_among_lines", morph_dose=mdose,
                                 tahoe_dose=tdose, r=float(lr.loc[lr.line == cid, "r"].iloc[0]),
                                 note=f"rank {int(lr.index[lr.line == cid][0]) + 1} of {len(lr)}; median line r {lr.r.median():.3f}"))
        # split-half reliability and uniform replicate reproducibility at every morphology dose
        for d in sorted(x for x in prof.dose_uM.unique() if x > 0):
            sh = split_half(wells, all_u, d)
            rp = replicate_repro(wells, all_u, d)
            rows.append(dict(tier=tier, label="RELIABILITY", analysis="morph_split_half", morph_dose=d,
                             r=sh["r_half"], r_sb=sh["r_sb"], n_units=sh["n"]))
            rows.append(dict(tier=tier, label="RELIABILITY", analysis="morph_replicate_reproducibility", morph_dose=d,
                             r=rp["frac_reproducible"], note=f"median replicate mAP {rp['median_mAP']:.3f}", n_units=rp["n"]))
        # exact-achiral sensitivity (primary configuration)
        ach = ut.loc[ut.exact_achiral, "unit"].tolist()
        compare(tier, mdose, tdose, line, ach, prof, tahoe, units, args, "SENS_exact_achiral", rows, nulls,
                secondary=False)
        # secondary: consistency vs agreement (primary configuration)
        if out is not None:
            consistency_test(tier, out, tdose, args, rows, scatter)
            consistency_test(tier, out, tdose, args, rows, scatter, ccol="consistency_centred")
        # T1 dose response: every LINCS dose x nearest Tahoe dose (+ full 6x3 grid, A549 and consensus)
        if tier == "T1_LINCS_A549":
            for ld in LINCS_DOSES:
                for td in TAHOE_DOSES:
                    nearest = min(TAHOE_DOSES, key=lambda x: abs(np.log(x / ld)))
                    for sig in ("A549", "consensus"):
                        lab = "DOSE_matched" if td == nearest else "DOSE_grid"
                        o = compare(tier, ld, td, sig, all_u, prof, tahoe, units, args, lab, rows, nulls,
                                    secondary=False)
            if act is not None:
                lowest_active(all_u, act, tahoe, rows, args)
                if "block_reliable" in act:
                    rel = sorted(set(act.loc[act.block_reliable.astype(bool), "unit"]) & set(all_u))
                    lowest_active(rel, act, tahoe, rows, args, subset="__block_reliable")
        # T3 low dose pair (0.625 vs 0.5) as secondary
        if tier == "T3_JUMP_U2OS":
            compare(tier, 0.625, 0.5, "consensus", all_u, prof, tahoe, units, args, "SECONDARY_low_dose", rows, nulls)

    # intersections: identical composition, 10 uM morphology in every tier vs Tahoe 5 uM
    sets = {t: set(units.loc[units.tier == t, "unit"]) for t in active_tiers}
    combos = {"T1T2T3": ["T1_LINCS_A549", "T2_EUOS_HepG2", "T3_JUMP_U2OS"], "T1T3": ["T1_LINCS_A549", "T3_JUMP_U2OS"],
              "T1T2": ["T1_LINCS_A549", "T2_EUOS_HepG2"], "T2T3": ["T2_EUOS_HepG2", "T3_JUMP_U2OS"]}
    for name, ts in combos.items():
        if not all(t in sets for t in ts):
            continue
        inter = sorted(set.intersection(*[sets[t] for t in ts]))
        for t in ts:
            prof = morphs[t][0]
            for sig in dict.fromkeys((TIERS[t]["sig_line"], "consensus")):
                compare(t, 10.0, 5.0, sig, inter, prof, tahoe, units, args, f"INTERSECT_{name}", rows, nulls,
                        secondary=False)

    tahoe_reliability(rows, units, morphs, tahoe, args)
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "all_results.csv", index=False)
    for t in res.tier.dropna().unique():
        res[res.tier == t].to_csv(OUT / f"{t}_results.csv", index=False)
    np.savez_compressed(OUT / "nulls.npz", **{k: v for k, v in nulls.items() if v is not None})
    if scatter:
        pd.concat(scatter).to_csv(OUT / "consistency_vs_agreement.csv", index=False)
    calibration(res, active_tiers)
    figures(res, nulls, scatter)
    accounting()
    prim = res[(res.label == "PRIMARY") & res.analysis.isin(["annotated_within_moa", "all_permutation"])]
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(prim[["tier", "analysis", "n_units", "n_classes_ge3", "r", "p", "ci_low", "ci_high", "null_mean", "null_q95"]].to_string(index=False))


def calibration(res, tiers):
    f = S2 / "cpjump1_calibration.csv"
    if not f.exists():
        die("missing step2/cpjump1_calibration.csv")
    cal = pd.read_csv(f)
    cal.to_csv(OUT / "CAL_cpjump1_calibration.csv", index=False)


def figures(res, nulls, scatter):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # 1) agreement per tier: primary within-MoA r with CI and its null distribution
    p = res[(res.label == "PRIMARY") & (res.analysis == "annotated_within_moa") & res.r.notna()]
    if len(p):
        fig, ax = plt.subplots(figsize=(6.5, 3.8))
        for i, r in enumerate(p.itertuples()):
            key = next((k for k in nulls if k.startswith(f"{r.tier}|PRIMARY|") and k.endswith("|within")), None)
            if key is not None:
                v = ax.violinplot(nulls[key], positions=[i], widths=.7, showextrema=False)
                for b in v["bodies"]:
                    b.set_alpha(.25)
            ax.errorbar(i, r.r, yerr=[[r.r - r.ci_low], [r.ci_high - r.r]], fmt="o", c="#1a73e8", capsize=4)
            ax.text(i, r.ci_high + .01, f"p={r.p:.3g}", ha="center", fontsize=8)
        ax.axhline(0, c="k", lw=.6)
        ax.set_xticks(range(len(p)), [t.replace("_", "\n", 1) for t in p.tier], fontsize=8)
        ax.set_ylabel("Mantel r (annotated units)\nnull: permute within MoA")
        ax.set_title("Morphology vs transcriptome agreement, primary test per tier", fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG / "agreement_per_tier.png", dpi=150)
    # 2) T1 dose response
    d = res[(res.tier == "T1_LINCS_A549") & (res.label == "DOSE_matched") & (res.analysis == "all_permutation")]
    if len(d):
        fig, ax = plt.subplots(figsize=(6, 3.8))
        for sig, c in (("A549", "#1a73e8"), ("consensus", "#e8710a")):
            s = d[d.tahoe_sig == sig].sort_values("morph_dose")
            ax.errorbar(s.morph_dose, s.r, yerr=[s.r - s.ci_low, s.ci_high - s.r], fmt="o-", c=c, capsize=3,
                        label=f"Tahoe {sig} (nearest dose)")
        ax.set_xscale("log")
        ax.axhline(0, c="k", lw=.6)
        ax.set_xlabel("LINCS concentration (µM)")
        ax.set_ylabel("Mantel r, all units")
        ax.set_title("T1 LINCS A549: agreement vs concentration", fontsize=10)
        ax.legend(fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(FIG / "T1_dose_response.png", dpi=150)
    # 3) consistency vs agreement
    if scatter:
        sc = pd.concat(scatter)
        ts = sc.tier.unique()
        fig, axes = plt.subplots(1, len(ts), figsize=(4.2 * len(ts), 3.8), sharey=True, squeeze=False)
        for ax, t in zip(axes[0], ts):
            s = sc[sc.tier == t]
            a = s.moa != "unclear"
            ax.scatter(s.consistency[~a], s.rowwise_agreement[~a], s=12, c="#9aa0a6", label="MoA unclear")
            ax.scatter(s.consistency[a], s.rowwise_agreement[a], s=14, c="#1a73e8", label="annotated")
            rr = res[(res.tier == t) & (res.analysis == "consistency_vs_rowwise_agreement")].iloc[0]
            ax.set_title(f"{t}\nρ={rr.r:.2f} [{rr.ci_low:.2f}, {rr.ci_high:.2f}]", fontsize=9)
            ax.axhline(0, c="k", lw=.6)
            ax.set_xlabel("Tahoe cross-line consistency")
        axes[0][0].set_ylabel("per-compound agreement\n(row-wise Spearman)")
        axes[0][0].legend(fontsize=7, frameon=False)
        fig.tight_layout()
        fig.savefig(FIG / "consistency_vs_agreement.png", dpi=150)


def accounting():
    tot, peaks = 0, {}
    for f in S2.glob("*download_summary.json"):
        j = json.loads(f.read_text())
        tot += float(j.get("bytes", j.get("total_bytes_streamed", 0)) or 0)
        pk = [v for k, v in j.items() if "peak" in k.lower() and isinstance(v, (int, float))]
        peaks[f.stem] = float(max(pk)) if pk else 0.0
    du = sum(p.stat().st_size for p in ROOT.rglob("*") if p.is_file())
    s = dict(total_bytes_downloaded=tot, per_pipeline_peak_disk=peaks, peak_disk_upper_bound=sum(peaks.values()) + du,
             project_disk_now=du)
    (OUT / "accounting.json").write_text(json.dumps(s, indent=1))
    print(f"TOTAL BYTES DOWNLOADED: {tot / 1e9:.2f} GB; project size now {du / 1e9:.2f} GB; "
          f"peak disk upper bound {s['peak_disk_upper_bound'] / 1e9:.2f} GB")


if __name__ == "__main__":
    main()
