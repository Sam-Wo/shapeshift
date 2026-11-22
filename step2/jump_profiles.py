"""JUMP Cell Painting morphology profiles for the 175 working-set compounds.

Reproducible end to end:
    py -3.13 -I -X utf8 step2/jump_profiles.py

Pipeline (agreed design; see step2/jump_notes.md, which this script writes):
  0. Stream per-plate CellProfiler well-level profiles from the public Cell Painting
     Gallery bucket (cpg0016-jump), keep only DMSO + working-set wells, cache that
     subset per plate in step2/jump_cache/ (safe to delete; re-runs re-download).
  1. Two separate datasets: HIGH = 10 uM (all non-source_7 COMPOUND plates holding a
     working-set compound), LOW = 0.625 uM (source_7 plates).
  2. Feature intersection across plates of the dataset; drop Metadata_*; drop features
     NaN/inf in >1% of wells or constant; impute remaining NaN with plate DMSO median.
  3. Per-plate robust z vs DMSO: (x - med_DMSO) / max(1.4826*MAD_DMSO, floor); clip +-50.
  4. Feature selection on normalized wells (compounds + DMSO): pycytominer blocklist,
     outlier features (max |z| > 500 before clipping), variance threshold
     (pycytominer freq_cut=0.05 / unique_cut=0.01), greedy Pearson |r| > 0.9 pruning.
  5. PCA whitening fitted on DMSO wells only (99% variance kept, eigenvalue + eps).
  6. Consensus: HIGH = median within source, mean across sources; LOW = median.
  7. Activity: mean pairwise cosine of replicate wells vs a plate-matched DMSO null.
"""
from __future__ import annotations

import io
import json
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------- config
ROOT = Path(__file__).resolve().parent.parent
STEP2 = ROOT / "step2"
CACHE = STEP2 / "jump_cache"
MANIFEST = STEP2 / "jump_download_manifest.csv"
DL_SUMMARY = STEP2 / "jump_download_summary.json"
URL = ("https://cellpainting-gallery.s3.amazonaws.com/cpg0016-jump/"
       "{source}/workspace/profiles/{batch}/{plate}/{plate}.parquet")
DMSO = "JCP2022_033924"
LOW_SOURCE = "source_7"
DOSE_UM = {"high": 10.0, "low": 0.625}
N_THREADS = 16

NAN_FRAC_MAX = 0.01      # drop features NaN/inf in more than this fraction of wells
MAD_SCALE = 1.4826
MAD_FLOOR_REL = 0.01     # floor = this x (median over plates of the plate DMSO scale)
MAD_FLOOR_ABS = 1e-12
OUTLIER_Z = 500.0        # pre-clip max |z| above which a feature is dropped
CLIP_Z = 50.0
FREQ_CUT, UNIQUE_CUT = 0.05, 0.01
CORR_THRESH = 0.9
VAR_KEEP = 0.99          # PCA-whitening: keep components explaining this DMSO variance
SPHERE_EPS = 3.0         # added to eigenvalues (DMSO features ~unit variance); chosen from
EPS_GRID = (0.1, 1.0, 3.0, 10.0)  # this diagnostic grid: smallest eps whose replicate AUC beats
                                  # the unsphered features at both doses (grid logged in notes)
N_NULL = 2000            # draws for the plate-matched DMSO null
P_ACTIVE = 0.05
SEED = 0

CHANNELS = ["AGP", "DNA", "ER", "Mito", "RNA"]
# pycytominer/data/blocklist_features.txt (55 features), regenerated from its pattern
BLOCKLIST = (
    [f"Nuclei_Correlation_Manders_{a}_{b}" for a, b in permutations(CHANNELS, 2)]
    + [f"Nuclei_Correlation_RWC_{a}_{b}" for a, b in permutations(CHANNELS, 2)]
    + [f"Nuclei_Granularity_{g}_{c}" for g in (14, 15, 16) for c in CHANNELS]
)
assert len(BLOCKLIST) == 55

LOG_LINES: list[str] = []
KEEP_CACHE = False


def log(msg: str = "") -> None:
    print(msg, flush=True)
    LOG_LINES.append(msg)


# ----------------------------------------------------------------------------- metadata
def load_metadata() -> pd.DataFrame:
    """Wells to keep (DMSO + working set) on COMPOUND plates that hold a working-set cpd."""
    cpds = pd.read_csv(STEP2 / "compounds.csv")
    ws = set(cpds.compound_id)
    assert len(ws) == 175, len(ws)
    well = pd.read_csv(ROOT / "data/jump/well.csv.gz")
    plate = pd.read_csv(ROOT / "data/jump/plate.csv.gz")
    plate = plate[plate.Metadata_PlateType == "COMPOUND"]
    assert not plate.Metadata_Plate.duplicated().any(), "plate names not unique"
    m = well.merge(plate[["Metadata_Source", "Metadata_Plate", "Metadata_Batch"]],
                   on=["Metadata_Source", "Metadata_Plate"])
    plates_ws = set(m.loc[m.Metadata_JCP2022.isin(ws), "Metadata_Plate"])
    m = m[m.Metadata_Plate.isin(plates_ws) & (m.Metadata_JCP2022.isin(ws) | (m.Metadata_JCP2022 == DMSO))]
    m = m.rename(columns={"Metadata_Source": "source", "Metadata_Plate": "plate",
                          "Metadata_Well": "well", "Metadata_JCP2022": "jcp",
                          "Metadata_Batch": "batch"}).reset_index(drop=True)
    m["dataset"] = np.where(m.source == LOW_SOURCE, "low", "high")
    m["is_dmso"] = m.jcp == DMSO
    assert not m.duplicated(["source", "plate", "well"]).any()
    return m


# ----------------------------------------------------------------------------- download
def _get(url: str, tries: int = 5) -> bytes:
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            if isinstance(e, urllib.error.HTTPError) and e.code == 404:
                raise
            if k == tries - 1:
                raise
            time.sleep(2 ** k)
    raise RuntimeError("unreachable")


def fetch_plate(source: str, batch: str, plate: str, wells: list[str]) -> dict:
    out = CACHE / f"{source}__{plate}.parquet"
    if out.exists():
        return {"plate": plate, "status": "cached"}
    url = URL.format(source=source, batch=batch, plate=plate)
    t0 = time.time()
    try:
        raw = _get(url)
    except urllib.error.HTTPError as e:
        return {"plate": plate, "status": f"http_{e.code}", "url": url}
    df = pd.read_parquet(io.BytesIO(raw))
    n_all = len(df)
    df = df[df.Metadata_Well.isin(set(wells))]
    tmp = out.with_suffix(".tmp")
    df.to_parquet(tmp, index=False)
    tmp.replace(out)
    return {"plate": plate, "source": source, "status": "downloaded", "bytes": len(raw),
            "seconds": round(time.time() - t0, 2), "wells_in_file": n_all, "wells_kept": len(df)}


def download_all(meta: pd.DataFrame, record: bool = True) -> pd.DataFrame:
    CACHE.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(MANIFEST) if MANIFEST.exists() else pd.DataFrame()
    jobs = meta.groupby(["source", "batch", "plate"]).well.apply(list).reset_index()
    log(f"plates to process: {len(jobs)}")
    t0 = time.time()
    recs = []
    with ThreadPoolExecutor(N_THREADS) as ex:
        futs = [ex.submit(fetch_plate, r.source, r.batch, r.plate, r.well) for r in jobs.itertuples()]
        for i, f in enumerate(as_completed(futs), 1):
            rec = f.result()
            if rec["status"] != "cached":
                recs.append(rec)
            if i % 25 == 0 or i == len(futs):
                gb = sum(r.get("bytes", 0) for r in recs) / 1e9
                print(f"  {i}/{len(futs)} plates, {gb:.2f} GB new, {time.time() - t0:.0f}s", flush=True)
    new = pd.DataFrame(recs)
    if len(new):
        man = pd.concat([man, new], ignore_index=True)
        man = man.drop_duplicates("plate", keep="last")
        man.to_csv(MANIFEST, index=False)
        if record and len(new) == len(jobs):  # full fresh download: record its wall time
            DL_SUMMARY.write_text(json.dumps({"plates": len(new), "bytes": float(new.bytes.sum()),
                                              "wall_seconds": round(time.time() - t0, 1)}), encoding="utf-8")
    bad =man[man.status != "downloaded"] if len(man) else man
    if len(bad):
        log(f"WARNING failed plates: {bad.to_dict('records')}")
    return man


# ----------------------------------------------------------------------------- loading
def load_dataset(meta: pd.DataFrame, name: str) -> tuple[pd.DataFrame, np.ndarray, list[str], dict]:
    md = meta[meta.dataset == name]
    frames, colsets = [], {}
    for (source, plate), _ in md.groupby(["source", "plate"]):
        f = CACHE / f"{source}__{plate}.parquet"
        if not f.exists():
            continue
        d = pd.read_parquet(f)
        frames.append(d)
        colsets[(source, plate)] = [c for c in d.columns if not c.startswith("Metadata_")]
    allfeat = sorted(set().union(*colsets.values()))
    common = set(allfeat)
    by_source_missing = {}
    for (source, _), cols in colsets.items():
        miss = set(allfeat) - set(cols)
        common &= set(cols)
        by_source_missing.setdefault(source, set()).update(miss)
    order = [c for c in colsets[next(iter(colsets))] if c in common]
    info = {"n_features_union": len(allfeat), "n_features_intersection": len(order),
            "missing_by_source": {s: len(v) for s, v in sorted(by_source_missing.items()) if v}}
    prof = pd.concat([d[["Metadata_Source", "Metadata_Plate", "Metadata_Well"] + order] for d in frames],
                     ignore_index=True)
    prof = prof.rename(columns={"Metadata_Source": "source", "Metadata_Plate": "plate", "Metadata_Well": "well"})
    assert not prof.duplicated(["source", "plate", "well"]).any(), "duplicate wells in profiles"
    keyed = md.merge(prof[["source", "plate", "well"]].assign(_row=np.arange(len(prof))),
                     on=["source", "plate", "well"], how="left", indicator=True)
    missing = keyed[keyed._merge == "left_only"]
    info["expected_wells"] = len(md)
    info["missing_wells"] = len(missing)
    info["missing_wells_detail"] = (missing.groupby(["source", "plate"]).size().to_dict() if len(missing) else {})
    info["missing_cpd_wells"] = int((~missing.is_dmso).sum())
    keyed = keyed[keyed._merge == "both"].reset_index(drop=True)
    X = prof[order].to_numpy(np.float32)[keyed._row.to_numpy()]
    obs = keyed.drop(columns=["_row", "_merge"])
    return obs, X, order, info


# ----------------------------------------------------------------------------- processing
def plate_groups(obs: pd.DataFrame) -> list[tuple[np.ndarray, np.ndarray]]:
    out = []
    for _, idx in obs.groupby("plate").indices.items():
        idx = np.sort(idx)
        out.append((idx, idx[obs.is_dmso.to_numpy()[idx]]))
    return out


def normalize(obs: pd.DataFrame, X: np.ndarray, feats: list[str], info: dict):
    """Steps 2-3. Returns clipped z (float32), pre-clip max|z| per feature, kept features."""
    finite = np.isfinite(X)
    nan_frac = 1 - finite.mean(0)
    X = np.where(finite, X, np.nan).astype(np.float32)
    std = np.nanstd(X, axis=0)
    drop_nan = nan_frac > NAN_FRAC_MAX
    drop_const = (~drop_nan) & ~(std > 0)
    keep = ~(drop_nan | drop_const)
    info["drop_nan_gt1pct"] = int(drop_nan.sum())
    info["drop_constant"] = int(drop_const.sum())
    X = X[:, keep]
    feats = [f for f, k in zip(feats, keep) if k]
    info["n_after_nan_const"] = len(feats)

    groups = plate_groups(obs)
    meds = np.empty((len(groups), X.shape[1]), np.float64)
    scales = np.empty_like(meds)
    for g, (_, didx) in enumerate(groups):
        assert len(didx) >= 8, f"too few DMSO wells on plate ({len(didx)})"
        D = X[didx].astype(np.float64)
        med = np.nanmedian(D, axis=0)
        meds[g] = med
        scales[g] = MAD_SCALE * np.nanmedian(np.abs(D - med), axis=0)
    ref_scale = np.nanmedian(scales, axis=0)
    floor = np.maximum(MAD_FLOOR_REL * np.nan_to_num(ref_scale, nan=0.0), MAD_FLOOR_ABS)
    floored = scales < floor[None, :]
    info["plate_feature_scales_floored_frac"] = float(np.nanmean(floored))
    scales = np.where(np.isfinite(scales) & ~floored, scales, floor[None, :])
    Z = np.empty_like(X)
    n_imputed = 0
    for g, (idx, _) in enumerate(groups):
        x = X[idx].astype(np.float64)
        nan = np.isnan(x)
        n_imputed += int(nan.sum())
        med = meds[g]
        x = np.where(nan, med[None, :], x)              # impute with plate DMSO median
        z = (x - med[None, :]) / scales[g][None, :]
        z[~np.isfinite(z)] = 0.0                        # plate DMSO median itself NaN
        Z[idx] = z.astype(np.float32)
    info["n_values_imputed"] = n_imputed
    maxabs = np.abs(Z).max(0)
    np.clip(Z, -CLIP_Z, CLIP_Z, out=Z)
    return Z, maxabs, feats


def variance_threshold_drop(Z: np.ndarray) -> np.ndarray:
    """pycytominer.variance_threshold: drop if (2nd most common / most common) < freq_cut
    or (n_unique / n) < unique_cut."""
    n = Z.shape[0]
    drop = np.zeros(Z.shape[1], bool)
    for j in range(Z.shape[1]):
        _, cnt = np.unique(Z[:, j], return_counts=True)
        if len(cnt) < 2:
            drop[j] = True
            continue
        top2 = np.partition(cnt, len(cnt) - 2)[-2:]
        if top2[0] / top2[1] < FREQ_CUT or len(cnt) / n < UNIQUE_CUT:
            drop[j] = True
    return drop


def correlation_prune(Z: np.ndarray) -> np.ndarray:
    """Greedy |r| > thresh pruning. Features are visited in ascending order of their total
    |r| with all features (pycytominer's ranking: least redundant first) and kept if their
    |r| with every already-kept feature is <= thresh. Returns a keep mask."""
    R = np.abs(np.corrcoef(Z, rowvar=False))
    R[~np.isfinite(R)] = 0
    order = np.argsort(R.sum(0), kind="stable")
    keep = np.zeros(Z.shape[1], bool)
    kept_idx: list[int] = []
    for j in order:
        if not kept_idx or R[j, kept_idx].max() <= CORR_THRESH:
            keep[j] = True
            kept_idx.append(j)
    return keep


def select_features(Z, maxabs, feats, info):
    f = np.array(feats)
    m_block = np.isin(f, BLOCKLIST)
    info["drop_blocklist"] = int(m_block.sum())
    m_out = (~m_block) & (maxabs > OUTLIER_Z)
    info["drop_outlier_gt500"] = int(m_out.sum())
    keep = ~(m_block | m_out)
    var_drop = np.zeros(len(f), bool)
    var_drop[keep] = variance_threshold_drop(Z[:, keep])
    info["drop_low_variance"] = int(var_drop.sum())
    keep &= ~var_drop
    corr_keep = np.zeros(len(f), bool)
    corr_keep[keep] = correlation_prune(Z[:, keep])
    info["drop_correlated"] = int((keep & ~corr_keep).sum())
    keep &= corr_keep
    info["n_selected"] = int(keep.sum())
    return Z[:, keep], f[keep].tolist()


def fit_sphere(D: np.ndarray, info: dict, eps: float | None = None):
    eps = SPHERE_EPS if eps is None else eps
    mu = D.mean(0)
    C = np.cov(D - mu, rowvar=False)
    lam, U = np.linalg.eigh(C)
    lam, U = lam[::-1].clip(min=0), U[:, ::-1]
    cum = np.cumsum(lam) / lam.sum()
    k = int(np.searchsorted(cum, VAR_KEEP) + 1)
    W = U[:, :k] / np.sqrt(lam[:k] + eps)
    info["sphere_dim"] = k
    info["sphere_eig_first_last"] = [float(lam[0]), float(lam[k - 1])]
    return mu, W


# ----------------------------------------------------------------------------- similarity
def unit(V: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(V, axis=1, keepdims=True)
    return V / np.where(n > 0, n, 1)


def replicate_auc(obs: pd.DataFrame, V: np.ndarray) -> dict:
    """AUC of replicate-pair cosine vs other-compound-pair cosine (compound wells only):
    all pairs, cross-plate pairs, and pairs at different positions (not same source+well;
    the position-robust metric used to choose the sphering eps); plus mean cosines."""
    c = ~obs.is_dmso.to_numpy()
    U = unit(V[c].astype(np.float64))
    lab, pl, wl, sr = (obs[k].to_numpy()[c] for k in ("jcp", "plate", "well", "source"))
    a, b = np.triu_indices(len(U), 1)
    s = (U @ U.T)[a, b]
    same = lab[a] == lab[b]
    sameplate = pl[a] == pl[b]
    samepos = (wl[a] == wl[b]) & (sr[a] == sr[b])

    def auc(pos, neg):
        r = pd.Series(np.concatenate([pos, neg])).rank().to_numpy()
        return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))

    return {"auc": auc(s[same], s[~same]),
            "auc_cross_plate": auc(s[same & ~sameplate], s[~same & ~sameplate]),
            "auc_diff_pos": auc(s[same & ~samepos], s[~same & ~samepos]),
            "mean_cos_replicate": float(s[same].mean()), "mean_cos_other_cpd": float(s[~same].mean()),
            "mean_cos_replicate_same_pos": float(s[same & samepos].mean()),
            "mean_cos_replicate_diff_pos": float(s[same & ~samepos].mean()),
            "n_rep_pairs": int(same.sum()), "n_rep_pairs_same_plate": int((same & sameplate).sum()),
            "n_rep_pairs_same_pos": int((same & samepos).sum())}


def masked_mean_cos(U: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """U: (..., n, k) unit vectors; mask: (n, n) symmetric bool of pairs to average (no diagonal)."""
    G = U @ np.swapaxes(U, -1, -2)
    return (G * mask).sum((-1, -2)) / mask.sum()


def _pval(null: np.ndarray, x: float) -> float:
    return float((1 + (null >= x).sum()) / (1 + len(null)))


def activity(obs: pd.DataFrame, V: np.ndarray, rng) -> pd.DataFrame:
    """Replicate consistency per compound (obs must have a 0..n-1 index).

    Score = mean pairwise cosine among the compound's replicate wells (all pairs), and the
    same restricted to pairs at different positions (not same source + well).
    Primary null (design): plate-matched DMSO - each replicate well replaced by a DMSO well
    of the same plate (without replacement within plate). Slot-aligned, so the same pair
    mask applies to the null.
    Stricter null: same-source other-compound wells - each replicate well replaced by a well
    of a different working-set compound from the same source, at most one well per compound
    (so a null set is never itself a replicate set). Falls back to the whole dataset for
    sources with too few other compounds (counted)."""
    U = unit(V.astype(np.float64))
    isd = obs.is_dmso.to_numpy()
    jcp, src = obs.jcp.to_numpy(), obs.source.to_numpy()
    dmso_by_plate = {p: np.flatnonzero((obs.plate.to_numpy() == p) & isd) for p in obs.plate.unique()}
    c = obs[~isd]
    wells_of = {k: c.index.to_numpy()[v] for k, v in c.groupby("jcp").indices.items()}
    cpds_by_src = {s: sorted(set(jcp[(src == s) & ~isd])) for s in c.source.unique()}
    all_cpds = sorted(wells_of)
    rows = []
    for cid, g in c.groupby("jcp"):
        idx = g.index.to_numpy()
        n = len(idx)
        row = {"compound_id": cid, "n_wells": n}
        if n < 2:
            rows.append(row)
            continue
        # pair masks (slot-aligned)
        w_, s_ = g.well.to_numpy(), g.source.to_numpy()
        samepos = (w_[:, None] == w_[None, :]) & (s_[:, None] == s_[None, :])
        off = ~np.eye(n, dtype=bool)
        m_all, m_dp = off, off & ~samepos
        has_dp = m_dp.any()
        x_all = float(masked_mean_cos(U[idx], m_all))
        x_dp = float(masked_mean_cos(U[idx], m_dp)) if has_dp else np.nan
        # plate-matched DMSO null; slots ordered as idx
        draws = np.empty((N_NULL, n), np.int64)
        for p in np.unique(g.plate.to_numpy()):
            slots = np.flatnonzero(g.plate.to_numpy() == p)
            pool = dmso_by_plate[p]
            keys = rng.random((N_NULL, len(pool)))
            draws[:, slots] = pool[np.argsort(keys, axis=1)[:, :len(slots)]]
        Ud = U[draws]
        null_all = masked_mean_cos(Ud, m_all)
        row.update({"rep_cos": x_all, "null_mean": float(null_all.mean()), "null_sd": float(null_all.std()),
                    "null_pct": float((null_all < x_all).mean() * 100), "p": _pval(null_all, x_all),
                    "n_pairs_diffpos": int(m_dp.sum() // 2), "rep_cos_diffpos": x_dp})
        if has_dp:
            null_dp = masked_mean_cos(Ud, m_dp)
            row["p_diffpos"] = _pval(null_dp, x_dp)
        # same-source other-compound null, one well per other compound
        draws = np.empty((N_NULL, n), np.int64)
        n_fb = 0
        for s in np.unique(s_):
            slots = np.flatnonzero(s_ == s)
            m = len(slots)
            others = [k for k in cpds_by_src[s] if k != cid]
            if len(others) < max(m, 3):
                others = [k for k in all_cpds if k != cid]
                n_fb += m
                pool_w = {k: wells_of[k] for k in others}
            else:
                pool_w = {k: wells_of[k][src[wells_of[k]] == s] for k in others}
            # random well of each other compound, per draw: (N_NULL, n_others)
            tab = np.stack([pool_w[k][rng.integers(0, len(pool_w[k]), N_NULL)] for k in others], 1)
            pick = np.argsort(rng.random((N_NULL, len(others))), axis=1)[:, :m]
            draws[:, slots] = np.take_along_axis(tab, pick, 1)
        Uc = U[draws]
        null_c = masked_mean_cos(Uc, m_all)
        row.update({"cpdnull_mean": float(null_c.mean()), "p_cpdnull": _pval(null_c, x_all),
                    "cpdnull_fallback_wells": n_fb})
        if has_dp:
            row["p_cpdnull_diffpos"] = _pval(masked_mean_cos(Uc, m_dp), x_dp)
        rows.append(row)
    return pd.DataFrame(rows)


def dmso_position_effect(obs: pd.DataFrame, V: np.ndarray, rng, n_max: int = 6000) -> dict:
    """Mean cosine of DMSO pairs on different plates of the same source, same vs other well."""
    idx = np.flatnonzero(obs.is_dmso.to_numpy())
    if len(idx) > n_max:
        idx = rng.choice(idx, n_max, replace=False)
    U = unit(V[idx].astype(np.float64))
    a, b = np.triu_indices(len(idx), 1)
    s = (U @ U.T)[a, b]
    pl, wl, sr = (obs[k].to_numpy()[idx] for k in ("plate", "well", "source"))
    base = (pl[a] != pl[b]) & (sr[a] == sr[b])
    same = wl[a] == wl[b]
    return {"dmso_cos_same_pos": round(float(s[base & same].mean()), 4),
            "dmso_cos_diff_pos": round(float(s[base & ~same].mean()), 4),
            "n_same_pos_pairs": int((base & same).sum())}


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    o = np.argsort(pv)
    r = pv[o] * len(pv) / np.arange(1, len(pv) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1].clip(max=1)
    out = np.empty_like(pv)
    out[o] = r
    q[ok] = out
    return q


# ----------------------------------------------------------------------------- main
def process(meta: pd.DataFrame, name: str, rng) -> dict:
    t0 = time.time()
    log(f"\n=== {name.upper()} ({DOSE_UM[name]} uM) ===")
    obs, X, feats, info = load_dataset(meta, name)
    info["n_plates"] = int(obs.plate.nunique())
    info["n_wells"] = len(obs)
    info["n_dmso_wells"] = int(obs.is_dmso.sum())
    info["n_cpd_wells"] = int((~obs.is_dmso).sum())
    info["n_compounds"] = int(obs.loc[~obs.is_dmso, "jcp"].nunique())
    info["sources"] = sorted(obs.source.unique())
    log(f"plates {info['n_plates']}, wells {info['n_wells']} (cpd {info['n_cpd_wells']}, "
        f"DMSO {info['n_dmso_wells']}), compounds {info['n_compounds']}, "
        f"missing wells {info['missing_wells']} (cpd {info['missing_cpd_wells']})")
    log(f"features: union {info['n_features_union']}, intersection {info['n_features_intersection']}, "
        f"missing by source {info['missing_by_source']}")
    Z, maxabs, feats = normalize(obs, X, feats, info)
    del X
    log(f"after NaN/const filter: {info['n_after_nan_const']} (dropped NaN>1% {info['drop_nan_gt1pct']}, "
        f"constant {info['drop_constant']}); imputed values {info['n_values_imputed']}; "
        f"plate x feature MAD floored {info['plate_feature_scales_floored_frac']:.4f}")
    Zs, sel = select_features(Z, maxabs, feats, info)
    del Z
    info["n_image_features_selected"] = sum(f.startswith("Image_") for f in sel)
    log(f"selection: blocklist -{info['drop_blocklist']}, outlier -{info['drop_outlier_gt500']}, "
        f"low variance -{info['drop_low_variance']}, correlated -{info['drop_correlated']} "
        f"-> {info['n_selected']} features ({info['n_image_features_selected']} Image_*)")
    dm = obs.is_dmso.to_numpy()
    mu, W = fit_sphere(Zs[dm].astype(np.float64), info)
    S = ((Zs.astype(np.float64) - mu) @ W).astype(np.float32)
    log(f"sphering: PCA-whitening on {dm.sum()} DMSO wells, {info['sphere_dim']} components "
        f"({VAR_KEEP:.0%} variance), eps {SPHERE_EPS}, eigenvalue range {info['sphere_eig_first_last']}")

    info["eps_grid_auc"] = {}
    for e in EPS_GRID:
        Wg = fit_sphere(Zs[dm].astype(np.float64), {}, eps=e)[1]
        r = replicate_auc(obs, (Zs.astype(np.float64) - mu) @ Wg)
        info["eps_grid_auc"][e] = {"all": round(r["auc"], 4), "diff_pos": round(r["auc_diff_pos"], 4)}
    log(f"eps grid replicate AUC: {info['eps_grid_auc']}")
    info["rep_unsphered"] = replicate_auc(obs, Zs)
    info["rep_sphered"] = replicate_auc(obs, S)
    for k in ("rep_unsphered", "rep_sphered"):
        r = info[k]
        log(f"{k}: AUC {r['auc']:.3f} (cross-plate {r['auc_cross_plate']:.3f}, diff-pos "
            f"{r['auc_diff_pos']:.3f}), mean cos rep {r['mean_cos_replicate']:.3f} (same pos "
            f"{r['mean_cos_replicate_same_pos']:.3f}, diff pos {r['mean_cos_replicate_diff_pos']:.3f}) "
            f"vs other cpd {r['mean_cos_other_cpd']:.3f}")

    fcols = [f"f{i}" for i in range(S.shape[1])]
    wells = obs[["jcp", "source", "batch", "plate", "well", "is_dmso"]].copy()
    wells.insert(0, "compound_id", np.where(obs.is_dmso, "DMSO", obs.jcp))
    wells = pd.concat([wells.reset_index(drop=True), pd.DataFrame(S, columns=fcols)], axis=1)
    wells.to_parquet(STEP2 / f"morph_wells_{name}.parquet", index=False)

    # consensus
    cw = wells[~wells.is_dmso]
    if name == "high":
        cons = cw.groupby(["compound_id", "source"])[fcols].median().groupby("compound_id").mean()
    else:
        cons = cw.groupby("compound_id")[fcols].median()
    g = cw.groupby("compound_id")
    meta_c = pd.DataFrame({
        "n_wells": g.size(),
        "n_plates": g.plate.nunique(),
        "sources": g.source.apply(lambda s: ";".join(sorted(set(s), key=lambda x: int(x.split("_")[1])))),
        "plates": g.plate.apply(lambda s: ";".join(sorted(set(s)))),
    })
    prof = meta_c.join(cons.astype(np.float32)).reset_index()
    prof.insert(1, "dose", name)
    prof.insert(2, "dose_uM", DOSE_UM[name])

    obs0 = obs.reset_index(drop=True)
    cw0 = obs0[~obs0.is_dmso]
    info["frac_cpd_wells_shared_pos"] = float(cw0.duplicated(["jcp", "source", "well"], keep=False).mean())
    act = activity(obs0, S, rng)
    act.insert(1, "dose", name)
    act["q_bh"] = bh(act.p.to_numpy())
    act["active"] = act.p < P_ACTIVE
    act["active_diffpos"] = (act.p_diffpos < P_ACTIVE).where(act.p_diffpos.notna())  # NaN: no such pairs
    act["active_vs_cpd"] = act.p_cpdnull < P_ACTIVE
    # strict: different-position pairs vs other-compound null (all pairs if no diff-pos pairs)
    act["strict_basis"] = np.where(act.p_cpdnull_diffpos.notna(), "diffpos", "allpairs")
    act["p_strict"] = act.p_cpdnull_diffpos.fillna(act.p_cpdnull)
    act["active_strict"] = act.p_strict < P_ACTIVE
    info["n_cpd_without_diffpos_pairs"] = int(act.p_diffpos.isna().sum())
    for k in ("active_diffpos", "active_vs_cpd", "active_strict"):
        info[f"frac_{k}"] = float(act[k].dropna().astype(bool).mean())  # over compounds with a value
    info["frac_cpdnull_fallback"] = float(act.cpdnull_fallback_wells.sum() / act.n_wells.sum())
    info["dmso_position_effect"] = dmso_position_effect(obs0, S, rng)
    log(f"DMSO position effect: {info['dmso_position_effect']}")
    log(f"stricter calls: diff-pos vs DMSO {info['frac_active_diffpos']:.1%}, all pairs vs other-cpd "
        f"{info['frac_active_vs_cpd']:.1%}, strict {info['frac_active_strict']:.1%}; compounds without "
        f"diff-pos pairs {info['n_cpd_without_diffpos_pairs']}; other-cpd fallback wells "
        f"{info['frac_cpdnull_fallback']:.1%}")
    info["n_active"] = int(act.active.sum())
    info["frac_active"] = float(act.active.mean())
    info["frac_active_q05"] = float((act.q_bh < 0.05).mean())
    info["median_null_pct"] = float(act.null_pct.median())
    info["frac_rep_above_null_mean"] = float((act.rep_cos > act.null_mean).mean())
    log(f"activity: {info['n_active']}/{len(act)} active at p<{P_ACTIVE} ({info['frac_active']:.1%}); "
        f"BH q<0.05 {info['frac_active_q05']:.1%}; replicate cos > null mean for "
        f"{info['frac_rep_above_null_mean']:.1%}")
    info["selected_features"] = sel
    info["seconds"] = round(time.time() - t0, 1)
    return {"info": info, "prof": prof, "act": act}


def write_notes(meta, man, res, t_total, n_bytes, cpds):
    L = []
    a = L.append
    a("# JUMP morphology profiles: notes (written by step2/jump_profiles.py)\n")
    a("## Data access and format quirks")
    a("- Path pattern verified: `https://cellpainting-gallery.s3.amazonaws.com/cpg0016-jump/"
      "{source}/workspace/profiles/{batch}/{plate}/{plate}.parquet` (anonymous HTTPS GET).")
    a("- Each file holds only `Metadata_Source`, `Metadata_Plate`, `Metadata_Well` plus 4,762 float32 "
      "CellProfiler features (source_15: 3,805; it lacks 957 features, mostly texture/correlation/"
      "radial-distribution and Image_* ones). Well IDs are zero-padded (`A01`), matching well.csv; "
      "1536-well plates (source_1, source_9) use `AA01`-style IDs for rows beyond Z.")
    a("- Profile files can contain wells absent from well.csv; we inner-join on source/plate/well.")
    a("- Some features are mostly NaN on some plates (e.g. `Cells_AreaShape_FormFactor`); handled by the 1% NaN rule.")
    a(f"- Plates downloaded: {(man.status == 'downloaded').sum()} of {meta.plate.nunique()} needed; "
      f"failures: {int((man.status != 'downloaded').sum())}.")
    a(f"- Bandwidth: {n_bytes / 1e9:.2f} GB of plate files (+~0.17 GB of manual probes while "
      f"developing). The filtered per-plate cache `step2/jump_cache/` "
      f"({sum(f.stat().st_size for f in CACHE.glob('*.parquet')) / 1e9:.2f} GB) is deleted at the end "
      "of the run unless `--keep-cache` is passed; a re-run then re-downloads. Per-plate download log: "
      "`step2/jump_download_manifest.csv`.")
    dl = man[man.status == "downloaded"]
    if "seconds" in dl and len(dl):
        a(f"- Sum of per-plate download times {dl.seconds.sum():.0f} s over {N_THREADS} threads.")
    if DL_SUMMARY.exists():
        d = json.loads(DL_SUMMARY.read_text(encoding="utf-8"))
        a(f"- Last full fresh download: {d['plates']} plates, {d['bytes'] / 1e9:.2f} GB in "
          f"{d['wall_seconds']:.0f} s wall ({d['bytes'] / 1e6 / d['wall_seconds']:.0f} MB/s).")
    a(f"- Wall time of this run: {t_total / 60:.1f} min (processing from the cache when it was present).\n")
    a("## Parameters")
    a(f"- NaN/inf rule: drop if > {NAN_FRAC_MAX:.0%} of wells; constant features dropped; remaining NaN -> plate DMSO median.")
    a(f"- Robust z: (x - median_DMSO_plate) / max(1.4826*MAD_DMSO_plate, floor); floor = "
      f"{MAD_FLOOR_REL} x median-over-plates of 1.4826*MAD (min {MAD_FLOOR_ABS}); clip +-{CLIP_Z:g}.")
    a(f"- Selection: pycytominer blocklist (55), outlier features with pre-clip max |z| > {OUTLIER_Z:g}, "
      f"variance threshold (freq_cut {FREQ_CUT}, unique_cut {UNIQUE_CUT}), greedy |r| > {CORR_THRESH} "
      "pruning (features visited in ascending total |r|, kept if |r| <= 0.9 with all kept ones).")
    a(f"- Sphering: PCA-whitening fit on DMSO wells only, centred on DMSO mean, components covering "
      f"{VAR_KEEP:.0%} of DMSO variance, scaled by 1/sqrt(eigenvalue + {SPHERE_EPS}). eps was chosen "
      f"from the grid {EPS_GRID} by the replicate-vs-other-compound AUC over different-position pairs "
      "(position-robust), requiring it to match or beat the unsphered selected features at both doses "
      "(eps 0.1 was worse than no sphering). This is a mild use of the replicate signal for tuning.")
    a("- Consensus: HIGH = median of wells within each source, then mean across sources; LOW = median of wells.")
    a(f"- Activity: mean pairwise cosine among a compound's replicate wells (sphered) vs a plate-matched "
      f"DMSO null ({N_NULL} draws: for every replicate well one DMSO well from the same plate, without "
      f"replacement within plate); p = (1 + #null >= obs)/(1 + {N_NULL}); active = p < {P_ACTIVE}. "
      "BH q is also reported.")
    a("- The primary DMSO test is close to saturated: other-compound pairs already have mean cosine "
      "~0.03-0.04 vs ~0 for DMSO pairs, so almost any set of compound wells beats DMSO. Stricter "
      "flags are therefore given alongside `active` (the caller chooses):")
    a("  - `p_diffpos` / `active_diffpos`: replicate score over pairs at different positions only "
      "(not same source+well), against the same plate-matched DMSO null with the same pair mask "
      "(controls positional similarity).")
    a("  - `p_cpdnull` / `active_vs_cpd`: all pairs vs a same-source other-compound null (each "
      "replicate well replaced by a well of a different working-set compound of the same source, at "
      "most one well per compound; whole-dataset pool for sources with too few compounds). "
      "Conservative: the null compounds are mostly bioactive drugs too.")
    a("  - `p_strict` / `active_strict`: different-position pairs vs the other-compound null "
      "(`strict_basis` = allpairs when a compound has no different-position pairs).\n")
    for name in ("high", "low"):
        i = res[name]["info"]
        a(f"## {name.upper()} ({DOSE_UM[name]} uM)")
        a(f"- Sources: {', '.join(i['sources'])}")
        a(f"- Plates {i['n_plates']}; wells {i['n_wells']} = {i['n_cpd_wells']} compound + "
          f"{i['n_dmso_wells']} DMSO; compounds {i['n_compounds']}/175.")
        a(f"- Expected wells missing from profile files: {i['missing_wells']} "
          f"(compound wells {i['missing_cpd_wells']}) {i['missing_wells_detail'] or ''}")
        a(f"- Features: union {i['n_features_union']} -> intersection {i['n_features_intersection']} "
          f"(missing by source: {i['missing_by_source'] or 'none'}) -> NaN>1% -{i['drop_nan_gt1pct']}, "
          f"constant -{i['drop_constant']} = {i['n_after_nan_const']} -> blocklist -{i['drop_blocklist']}, "
          f"outlier -{i['drop_outlier_gt500']}, low variance -{i['drop_low_variance']}, "
          f"correlated -{i['drop_correlated']} = **{i['n_selected']} selected** "
          f"({i['n_image_features_selected']} of them Image_* features).")
        a(f"- Imputed values: {i['n_values_imputed']}; fraction of plate x feature MADs at floor: "
          f"{i['plate_feature_scales_floored_frac']:.4f}.")
        a(f"- Sphering: **{i['sphere_dim']} dimensions** (eigenvalues {i['sphere_eig_first_last'][0]:.2f} "
          f"... {i['sphere_eig_first_last'][1]:.3f}); replicate AUC by eps: {i['eps_grid_auc']}.")
        for k, lab in (("rep_unsphered", "selected features, unsphered"), ("rep_sphered", "sphered")):
            r = i[k]
            a(f"- Replicate check ({lab}): replicate-pair vs other-compound-pair cosine AUC "
              f"{r['auc']:.3f} (cross-plate {r['auc_cross_plate']:.3f}; different-position pairs "
              f"{r['auc_diff_pos']:.3f}); mean cosine replicate {r['mean_cos_replicate']:.3f} (same "
              f"position {r['mean_cos_replicate_same_pos']:.3f}, different {r['mean_cos_replicate_diff_pos']:.3f}) "
              f"vs other compound {r['mean_cos_other_cpd']:.3f} ({r['n_rep_pairs']} replicate pairs, "
              f"{r['n_rep_pairs_same_pos']} same position, {r['n_rep_pairs_same_plate']} same plate).")
        a(f"- Active (p < {P_ACTIVE}): **{i['n_active']}/{i['n_compounds']} = {i['frac_active']:.1%}**; "
          f"BH q < 0.05: {i['frac_active_q05']:.1%}; replicate cosine above null mean for "
          f"{i['frac_rep_above_null_mean']:.1%}; median null percentile {i['median_null_pct']:.1f}.")
        pe = i["dmso_position_effect"]
        a(f"- Position effect: DMSO pairs on different plates of the same source have mean cosine "
          f"{pe['dmso_cos_same_pos']:.3f} at the same well vs {pe['dmso_cos_diff_pos']:.3f} at different "
          f"wells.")
        a(f"- Stricter calls: active_diffpos {i['frac_active_diffpos']:.1%} (of the "
          f"{i['n_compounds'] - i['n_cpd_without_diffpos_pairs']} compounds with different-position "
          f"replicate pairs), active_vs_cpd "
          f"{i['frac_active_vs_cpd']:.1%}, **active_strict {i['frac_active_strict']:.1%}**; compounds "
          f"without different-position pairs: {i['n_cpd_without_diffpos_pairs']}; other-compound null "
          f"used the whole-dataset pool for {i['frac_cpdnull_fallback']:.1%} of wells.")
        a(f"- Time: {i['seconds']:.0f} s.\n")
    miss = {n: sorted(set(cpds) - set(res[n]["prof"].compound_id)) for n in ("high", "low")}
    a("## Coverage")
    a(f"- Compounds missing at HIGH: {miss['high'] or 'none'}; at LOW: {miss['low'] or 'none'}.")
    a(f"- Caveat: {res['high']['info']['frac_cpd_wells_shared_pos']:.0%} (HIGH) and "
      f"{res['low']['info']['frac_cpd_wells_shared_pos']:.0%} (LOW) of compound wells share their well "
      "position with another replicate of the same compound (same source), so part of the replicate "
      "similarity can be positional (see the DMSO position effect and the `active_diffpos` / `active_strict` flags above). "
      "`well` is kept in the wells files.")
    a("\n## Outputs")
    a("- `morph_profiles.parquet` (root): both doses, f0..fN float32, NaN in columns beyond a dose's "
      "dimensionality. Also `morph_profiles_high.parquet` and `morph_profiles_low.parquet` (root) "
      "with only that dose's columns. f-columns are NOT comparable between doses (different spaces).")
    a("- `step2/morph_activity.csv`, `step2/morph_wells_{high,low}.parquet`, `step2/jump_run.json` "
      "(full run info incl. selected feature names).")
    (STEP2 / "jump_notes.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main():
    t_start = time.time()
    rng = np.random.default_rng(SEED)
    meta = load_metadata()
    cpds = sorted(pd.read_csv(STEP2 / "compounds.csv").compound_id)
    log(f"metadata: {meta.plate.nunique()} plates, {int((~meta.is_dmso).sum())} compound wells, "
        f"{int(meta.is_dmso.sum())} DMSO wells")
    global KEEP_CACHE
    KEEP_CACHE = "--keep-cache" in sys.argv
    only = [a for a in sys.argv[1:] if not a.startswith("--")]  # smoke test: download these plates only
    if only:
        download_all(meta[meta.plate.isin(only)], record=False)
        return
    man = download_all(meta)
    n_bytes = float(man.get("bytes", pd.Series(dtype=float)).fillna(0).sum())
    log(f"download: {n_bytes / 1e9:.2f} GB total, {time.time() - t_start:.0f}s so far")
    res = {name: process(meta, name, rng) for name in ("high", "low")}

    prof = pd.concat([res["high"]["prof"], res["low"]["prof"]], ignore_index=True)
    fcols = [c for c in prof.columns if c.startswith("f") and c[1:].isdigit()]
    fcols = sorted(fcols, key=lambda c: int(c[1:]))
    head = ["compound_id", "dose", "dose_uM", "n_wells", "n_plates", "sources", "plates"]
    prof = prof[head + fcols]
    prof[fcols] = prof[fcols].astype(np.float32)
    prof.to_parquet(ROOT / "morph_profiles.parquet", index=False)
    for name in ("high", "low"):
        p = res[name]["prof"]
        p[head + [c for c in fcols if c in p.columns]].to_parquet(ROOT / f"morph_profiles_{name}.parquet", index=False)
    act = pd.concat([res["high"]["act"], res["low"]["act"]], ignore_index=True)
    act = act.rename(columns={"rep_cos": "replicate_mean_cosine", "null_pct": "null_percentile"})
    act = act[["compound_id", "dose", "n_wells", "replicate_mean_cosine", "null_mean", "null_sd",
               "null_percentile", "p", "q_bh", "active", "n_pairs_diffpos", "rep_cos_diffpos",
               "p_diffpos", "active_diffpos", "cpdnull_mean", "p_cpdnull", "active_vs_cpd",
               "p_cpdnull_diffpos", "strict_basis", "p_strict", "active_strict", "cpdnull_fallback_wells"]]
    act.to_csv(STEP2 / "morph_activity.csv", index=False)

    t_total = time.time() - t_start
    run = {"seconds_total": round(t_total, 1), "download_bytes": n_bytes,
           "high": res["high"]["info"], "low": res["low"]["info"]}
    (STEP2 / "jump_run.json").write_text(json.dumps(run, indent=1, default=str), encoding="utf-8")
    write_notes(meta, man, res, t_total, n_bytes, cpds)
    if not KEEP_CACHE:
        for f in CACHE.glob("*.parquet"):
            f.unlink()
        log(f"deleted filtered plate cache {CACHE} (pass --keep-cache to keep it)")
    log(f"\ndone in {t_total / 60:.1f} min")


if __name__ == "__main__":
    main()
