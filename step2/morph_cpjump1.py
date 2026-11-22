"""CPJUMP1 cell-line calibration: how much does compound morphology change with the cell
line alone (A549 vs U2OS) under one protocol?

Reproducible end to end:
    py -3.13 -I -X utf8 step2/morph_cpjump1.py [--keep-cache]

Follows step2/MORPH_SPEC.md (binding). Summary:
  0. Resolve the gallery accession with gallery_index.resolve('cpjump1'). Fetch (and cache in
     step2/cpjump1_cache/) the CPJUMP1 experiment metadata (2024_Chandrasekaran_NatureMethods),
     the JUMP-Target-1 compound metadata + platemap (JUMP-Target) and the NatureMethods
     additional annotations (single designed `target` gene per compound).
  1. Plates: batch 2020_11_04_CPJUMP1, compound plates only, the benchmark's own protocol filter
     (Density 100, antibiotics absent, no Cas9 compound plates) -> 4 plates per line x time.
  2. Stream <plate>.csv.gz (raw well-level CellProfiler profiles), join the platemap on well.
  3. Four datasets (A549/U2OS x 24/48 h), each normalized / feature-selected / sphered against
     its own DMSO. Raw feature columns: one intersection per cell line across both time points.
  4. Consensus per compound (median of replicate wells), replicate mAP vs DMSO with a
     plate-matched DMSO null, norm-z activity.
  5. Calibration statistics per time point and compound set (all / active / matched):
     cross-line Mantel (Spearman) with permutation, within-target permutation and bootstrap CI;
     split-half within-line reliability (Spearman-Brown corrected); attenuation factor;
     per-compound row agreement; positional diagnostics.
"""
from __future__ import annotations

import hashlib
import io
import json
import re
import shutil
import sys
import time
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gallery_index as gi  # noqa: E402

# ----------------------------------------------------------------------------- config
ROOT = Path(__file__).resolve().parent.parent
STEP2 = ROOT / "step2"
CACHE = STEP2 / "cpjump1_cache"
PLATE_CACHE = CACHE / "plates"
BATCH = "2020_11_04_CPJUMP1"
DOSE_UM = 5.0
LINES = ["A549", "U2OS"]
TIMES = [24, 48]
PRIMARY_TIME = 48

GH = "https://raw.githubusercontent.com"
META_URLS = {
    "experiment-metadata.tsv":
        f"{GH}/jump-cellpainting/2024_Chandrasekaran_NatureMethods/main/benchmark/output/experiment-metadata.tsv",
    "JUMP-Target-1_compound_metadata_additional_annotations.tsv":
        f"{GH}/jump-cellpainting/2024_Chandrasekaran_NatureMethods/main/benchmark/input/"
        "JUMP-Target-1_compound_metadata_additional_annotations.tsv",
    "JUMP-Target-1_compound_metadata.tsv":
        f"{GH}/jump-cellpainting/JUMP-Target/master/JUMP-Target-1_compound_metadata.tsv",
    "JUMP-Target-1_compound_platemap.tsv":
        f"{GH}/jump-cellpainting/JUMP-Target/master/JUMP-Target-1_compound_platemap.tsv",
}

# Feature-selection scope. "line_time": select per cell line x time point (the four datasets);
# "line": select once per cell line on both time points pooled. Raw columns are always one
# intersection per cell line. None of the calibration statistics crosses time points.
SELECTION_SCOPE = "line_time"

NAN_FRAC_MAX = 0.01
MAD_SCALE = 1.4826
MAD_FLOOR_REL = 0.01
MAD_FLOOR_ABS = 1e-12
OUTLIER_Z = 500.0
CLIP_Z = 50.0
FREQ_CUT, UNIQUE_CUT = 0.05, 0.01
CORR_THRESH = 0.9
VAR_KEEP = 0.99
SPHERE_EPS = 3.0          # fixed a priori (same value as the JUMP tier); never tuned here
DROP_BRIGHTFIELD = True   # CPJUMP1 has extra brightfield channels; drop for a 5-channel space
POS_RE = re.compile(r"AreaShape_Center_|AreaShape_BoundingBoxM(ax|in)imum_|_Location_")
BF_RE = re.compile(r"Brightfield|HighZBF|LowZBF")

N_NULL_MAP = 10000        # DMSO-null draws for replicate mAP (spec: >= 1000)
N_PERM = 9999
N_BOOT = 1000
N_SPLITS = 100            # random split halves per line x time (spec: >= 50)
SPLITS_PER_BOOT = 10      # split halves evaluated inside each bootstrap draw (rotating)
P_ACTIVE = 0.05
Z_ACTIVE = 3.0
SEED = 20201104

CHANNELS = ["AGP", "DNA", "ER", "Mito", "RNA"]
BLOCKLIST = (
    [f"Nuclei_Correlation_Manders_{a}_{b}" for a, b in permutations(CHANNELS, 2)]
    + [f"Nuclei_Correlation_RWC_{a}_{b}" for a, b in permutations(CHANNELS, 2)]
    + [f"Nuclei_Granularity_{g}_{c}" for g in (14, 15, 16) for c in CHANNELS]
)
assert len(BLOCKLIST) == 55

LOG: list[str] = []
PEAK = {"disk": 0}
OUT_FILES: list[Path] = []


def log(msg: str = "") -> None:
    print(msg, flush=True)
    LOG.append(msg)


def dir_bytes(p: Path) -> int:
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file()) if p.exists() else 0


def track_disk() -> None:
    used = dir_bytes(CACHE) + sum(f.stat().st_size for f in OUT_FILES if f.exists())
    PEAK["disk"] = max(PEAK["disk"], used)


def get(url: str, tries: int = 5) -> bytes:
    last = None
    for k in range(tries):
        try:
            r = requests.get(url, timeout=180)
            if r.status_code == 200:
                gi.Bytes.add(len(r.content))
                return r.content
            last = RuntimeError(f"HTTP {r.status_code} for {url}")
            if r.status_code == 404:
                break
        except requests.RequestException as e:  # retry transient errors
            last = e
        time.sleep(2 ** k)
    raise RuntimeError(f"source unreachable: {url}: {last}")


# ----------------------------------------------------------------------------- metadata
def fetch_metadata() -> dict[str, pd.DataFrame]:
    CACHE.mkdir(parents=True, exist_ok=True)
    out, info = {}, {}
    for name, url in META_URLS.items():
        f = CACHE / name
        if not f.exists():
            raw = get(url)
            f.write_bytes(raw)
            info[name] = "downloaded"
        else:
            info[name] = "cached"
        out[name] = pd.read_csv(f, sep="\t")
        info[name + "_sha256"] = hashlib.sha256(f.read_bytes()).hexdigest()[:16]
    out["_info"] = info
    track_disk()
    return out


def select_plates(exp: pd.DataFrame, acc: str) -> pd.DataFrame:
    e = exp[(exp.Batch == BATCH) & (exp.Perturbation == "compound")].copy()
    e_all = e.copy()
    e = e[(e.Density == 100) & (e.Antibiotics == "absent") & (e.Cell_line != "Cas9")]
    e = e[e.Plate_Map_Name == "JUMP-Target-1_compound_platemap"]
    e = e.rename(columns={"Assay_Plate_Barcode": "plate", "Cell_type": "cell_line", "Time": "time_h"})
    e["time_h"] = e.time_h.astype(int)
    assert not e.plate.duplicated().any()
    counts = e.groupby(["cell_line", "time_h"]).size()
    assert set(counts.index) == {(c, t) for c in LINES for t in TIMES}, counts
    assert (counts == 4).all(), counts
    prefixes, _ = gi.s3_list(f"{acc}/source_4/workspace/profiles/{BATCH}/", "/")
    listed = {p.rstrip("/").split("/")[-1] for p in prefixes}
    miss = set(e.plate) - listed
    if miss:
        raise RuntimeError(f"plates not in gallery listing: {sorted(miss)}")
    excluded = e_all[~e_all.Assay_Plate_Barcode.isin(e.plate)]
    return e.reset_index(drop=True), excluded


def plate_key(acc: str, plate: str) -> tuple[str, int]:
    _, keys = gi.s3_list(f"{acc}/source_4/workspace/profiles/{BATCH}/{plate}/", "/")
    want = f"{acc}/source_4/workspace/profiles/{BATCH}/{plate}/{plate}.csv.gz"
    hit = [k for k in keys if k[0] == want]
    if not hit:
        raise RuntimeError(f"raw profile file missing: {want}; have {[k for k, _ in keys]}")
    return hit[0]


def fetch_plates(acc: str, plates: pd.DataFrame) -> pd.DataFrame:
    PLATE_CACHE.mkdir(parents=True, exist_ok=True)
    recs = []
    for r in plates.itertuples():
        out = PLATE_CACHE / f"{r.plate}.parquet"
        if out.exists():
            recs.append({"plate": r.plate, "status": "cached"})
            continue
        key, size = plate_key(acc, r.plate)
        t0 = time.time()
        raw = get(gi.s3_url(key))
        assert len(raw) == size, (len(raw), size)
        d = pd.read_csv(io.BytesIO(raw), compression="gzip")
        feats = [c for c in d.columns if not c.startswith("Metadata_")]
        d[feats] = d[feats].astype(np.float32)
        d.to_parquet(out, index=False)
        del raw
        recs.append({"plate": r.plate, "status": "downloaded", "key": key, "bytes": size,
                     "seconds": round(time.time() - t0, 1), "wells": len(d), "features": len(feats)})
        track_disk()
        log(f"  {r.plate} ({r.cell_line} {r.time_h}h): {size / 1e6:.1f} MB, {len(d)} wells")
    return pd.DataFrame(recs)


def load_line(line: str, plates: pd.DataFrame, pm: pd.DataFrame, ann: pd.DataFrame):
    """All plates of one cell line, raw-column intersection, platemap join."""
    pl = plates[plates.cell_line == line]
    frames, cols = [], []
    for r in pl.itertuples():
        d = pd.read_parquet(PLATE_CACHE / f"{r.plate}.parquet")
        assert (d.Metadata_Plate == r.plate).all()
        d["time_h"] = r.time_h
        frames.append(d)
        cols.append([c for c in d.columns if not c.startswith("Metadata_") and c != "time_h"])
    union = sorted(set().union(*map(set, cols)))
    common = [c for c in cols[0] if all(c in set(x) for x in cols)]
    d = pd.concat([f[["Metadata_Plate", "Metadata_Well", "time_h"] + common] for f in frames], ignore_index=True)
    d = d.rename(columns={"Metadata_Plate": "plate", "Metadata_Well": "well"})
    m = d[["plate", "well", "time_h"]].merge(pm, left_on="well", right_on="well_position", how="left",
                                             indicator=True)
    assert (m._merge == "both").all(), "wells not in platemap"
    assert not m.duplicated(["plate", "well"]).any()
    per_plate = m.groupby("plate").size()
    assert (per_plate <= 384).all(), per_plate
    missing = {p: sorted(set(pm.well_position) - set(m.well[m.plate == p])) for p in pl.plate}
    missing = {p: w for p, w in missing.items() if w}
    if missing:
        log(f"[{line}] wells absent from profile files: {missing} "
            f"(platemap content: {[pm.set_index('well_position').broad_sample.get(w) for ws in missing.values() for w in ws]})")
    m["is_dmso"] = m.broad_sample.isna()
    n_dm = m.groupby("plate").is_dmso.sum()
    assert (n_dm >= 32).all(), n_dm  # BR00116995 (U2OS 24 h) lacks 83 wells in the gallery file
    m["broad_sample"] = m.broad_sample.fillna("DMSO")
    m = m.merge(ann, on="broad_sample", how="left")
    assert m.loc[~m.is_dmso, "pert_iname"].notna().all(), "treated wells without compound metadata"
    obs = m.drop(columns=["_merge", "well_position", "solvent"]).reset_index(drop=True)
    obs["cell_line"] = line
    X = d[common].to_numpy(np.float32)
    info = {"n_raw_union": len(union), "n_raw_intersection": len(common), "missing_wells": missing}
    return obs, X, common, info


# ----------------------------------------------------------------------------- processing
def plate_groups(obs):
    out = []
    dm = obs.is_dmso.to_numpy()
    for _, idx in obs.groupby("plate").indices.items():
        idx = np.sort(idx)
        out.append((idx, idx[dm[idx]]))
    return out


def prefilter(feats: list[str], info: dict) -> np.ndarray:
    f = np.array(feats)
    bf = np.array([bool(BF_RE.search(x)) for x in f]) if DROP_BRIGHTFIELD else np.zeros(len(f), bool)
    pos = np.array([bool(POS_RE.search(x)) for x in f]) & ~bf
    info["drop_brightfield"] = int(bf.sum())
    info["drop_positional_coords"] = int(pos.sum())
    return ~(bf | pos)


def normalize(obs, X, feats, info):
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
    meds = np.empty((len(groups), X.shape[1]))
    scales = np.empty_like(meds)
    for g, (_, didx) in enumerate(groups):
        D = X[didx].astype(np.float64)
        med = np.nanmedian(D, axis=0)
        meds[g] = med
        scales[g] = MAD_SCALE * np.nanmedian(np.abs(D - med), axis=0)
    ref = np.nanmedian(scales, axis=0)
    floor = np.maximum(MAD_FLOOR_REL * np.nan_to_num(ref, nan=0.0), MAD_FLOOR_ABS)
    floored = scales < floor[None, :]
    info["plate_feature_scales_floored_frac"] = round(float(np.nanmean(floored)), 4)
    scales = np.where(np.isfinite(scales) & ~floored, scales, floor[None, :])
    Z = np.empty_like(X)
    n_imp = 0
    for g, (idx, _) in enumerate(groups):
        x = X[idx].astype(np.float64)
        nan = np.isnan(x)
        n_imp += int(nan.sum())
        x = np.where(nan, meds[g][None, :], x)
        z = (x - meds[g][None, :]) / scales[g][None, :]
        z[~np.isfinite(z)] = 0.0
        Z[idx] = z.astype(np.float32)
    info["n_values_imputed"] = n_imp
    maxabs = np.abs(Z).max(0)
    Zraw = Z.copy()                                     # unclipped, for the variance threshold
    np.clip(Z, -CLIP_Z, CLIP_Z, out=Z)
    return Z, maxabs, feats, Zraw


def variance_threshold_drop(Z):
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


def correlation_prune(Z):
    R = np.abs(np.corrcoef(Z, rowvar=False))
    R[~np.isfinite(R)] = 0
    order = np.argsort(R.sum(0), kind="stable")
    keep = np.zeros(Z.shape[1], bool)
    kept: list[int] = []
    for j in order:
        if not kept or R[j, kept].max() <= CORR_THRESH:
            keep[j] = True
            kept.append(j)
    return keep


def select_features(Z, maxabs, feats, info, Zraw):
    """Returns a keep mask over `feats`. The variance threshold runs on the unclipped z:
    on clipped z, features where many treated wells hit +-50 (strong responders) look
    'low-variance' because the clip value becomes the most frequent value."""
    f = np.array(feats)
    m_block = np.isin(f, BLOCKLIST)
    info["drop_blocklist"] = int(m_block.sum())
    m_out = (~m_block) & (maxabs > OUTLIER_Z)
    info["drop_outlier_gt500"] = int(m_out.sum())
    keep = ~(m_block | m_out)
    info["n_after_outlier"] = int(keep.sum())
    vd = np.zeros(len(f), bool)
    vd[keep] = variance_threshold_drop(Zraw[:, keep])
    info["drop_low_variance"] = int(vd.sum())
    keep &= ~vd
    info["n_after_variance"] = int(keep.sum())
    ck = np.zeros(len(f), bool)
    ck[keep] = correlation_prune(Z[:, keep])
    info["drop_correlated"] = int((keep & ~ck).sum())
    keep &= ck
    info["n_selected"] = int(keep.sum())
    return keep


def fit_sphere(D, info):
    mu = D.mean(0)
    C = np.cov(D - mu, rowvar=False)
    lam, U = np.linalg.eigh(C)
    lam, U = lam[::-1].clip(min=0), U[:, ::-1]
    cum = np.cumsum(lam) / lam.sum()
    k = int(np.searchsorted(cum, VAR_KEEP) + 1)
    W = U[:, :k] / np.sqrt(lam[:k] + SPHERE_EPS)
    info["sphere_dim"] = k
    info["sphere_eig_first_last"] = [round(float(lam[0]), 3), round(float(lam[k - 1]), 4)]
    info["sphere_n_dmso_fit"] = int(D.shape[0])
    return mu, W


def unit_rows(V):
    n = np.linalg.norm(V, axis=1, keepdims=True)
    return V / np.where(n > 0, n, 1)


# ----------------------------------------------------------------------------- replicate mAP
def avg_precision(sims: np.ndarray, pos: np.ndarray) -> float:
    order = np.argsort(-sims, kind="stable")
    hit = pos[order]
    ranks = np.flatnonzero(hit) + 1
    return float(np.mean(np.arange(1, len(ranks) + 1) / ranks))


def map_of_set(U, rep_idx, neg_idx):
    """Mean AP: each replicate queries (other replicates + negatives) by cosine."""
    aps = []
    for q in rep_idx:
        others = rep_idx[rep_idx != q]
        pool = np.concatenate([others, neg_idx])
        s = U[pool] @ U[q]
        pos = np.zeros(len(pool), bool)
        pos[:len(others)] = True
        aps.append(avg_precision(s, pos))
    return float(np.mean(aps))


def replicate_map(obs, U, rng):
    """Per compound: replicate mAP vs DMSO wells on the same plates; p from a plate-matched
    null of random DMSO sets of the same size and plate composition."""
    dm = obs.is_dmso.to_numpy()
    plate = obs.plate.to_numpy()
    dmso_by_plate = {p: np.flatnonzero(dm & (plate == p)) for p in np.unique(plate)}
    rows, null_cache = [], {}
    trt = obs[~obs.is_dmso]
    for bs, idx in trt.groupby("broad_sample").indices.items():
        idx = trt.index.to_numpy()[idx]
        pl_counts = tuple(sorted(pd.Series(plate[idx]).value_counts().items()))
        neg = np.concatenate([dmso_by_plate[p] for p, _ in pl_counts])
        obs_map = map_of_set(U, idx, neg)
        if pl_counts not in null_cache:
            null = np.empty(N_NULL_MAP)
            for b in range(N_NULL_MAP):
                fake = np.concatenate([rng.choice(dmso_by_plate[p], c, replace=False) for p, c in pl_counts])
                rest = np.setdiff1d(neg, fake, assume_unique=True)
                null[b] = map_of_set(U, fake, rest)
            null_cache[pl_counts] = null
        null = null_cache[pl_counts]
        p = (1 + np.sum(null >= obs_map)) / (1 + N_NULL_MAP)
        rows.append({"broad_sample": bs, "rep_mAP": obs_map, "rep_p": p,
                     "null_mAP_mean": float(null.mean())})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- one dataset
def process_line(line, plates, pm, ann, rng):
    obs_all, X_all, feats0, linfo = load_line(line, plates, pm, ann)
    keep0 = prefilter(feats0, linfo)
    feats0 = [f for f, k in zip(feats0, keep0) if k]
    X_all = X_all[:, keep0]
    linfo["n_after_prefilter"] = len(feats0)
    log(f"[{line}] raw features union {linfo['n_raw_union']}, intersection {linfo['n_raw_intersection']}, "
        f"-brightfield {linfo['drop_brightfield']}, -positional {linfo['drop_positional_coords']} "
        f"= {len(feats0)}")
    # normalization per dataset (line x time); the per-plate robust z is the same either way
    parts = {}
    for t in TIMES:
        m = (obs_all.time_h == t).to_numpy()
        obs = obs_all[m].reset_index(drop=True)
        info = {}
        Z, maxabs, feats, Zraw = normalize(obs, X_all[m], feats0, info)
        parts[t] = dict(obs=obs, Z=Z, maxabs=maxabs, feats=feats, info=info, Zraw=Zraw)
    if SELECTION_SCOPE == "line":
        common = sorted(set.intersection(*[set(p["feats"]) for p in parts.values()]), key=feats0.index)
        Zs, mx = [], []
        for p in parts.values():
            ix = [p["feats"].index(f) for f in common]
            p["Z"], p["maxabs"], p["feats"] = p["Z"][:, ix], p["maxabs"][ix], common
            p["Zraw"] = p["Zraw"][:, ix]
            Zs.append(p["Z"])
            mx.append(p["maxabs"])
        sinfo = {}
        keep = select_features(np.vstack(Zs), np.max(mx, axis=0), common, sinfo,
                               np.vstack([p["Zraw"] for p in parts.values()]))
        for p in parts.values():
            p["info"].update(sinfo)
            p["keep"] = keep
    else:
        for p in parts.values():
            p["keep"] = select_features(p["Z"], p["maxabs"], p["feats"], p["info"], p["Zraw"])
    out = {}
    for t, p in parts.items():
        obs, info = p["obs"], p["info"]
        p.pop("Zraw", None)
        info.update(linfo)
        Zs = p["Z"][:, p["keep"]]
        sel = [f for f, k in zip(p["feats"], p["keep"]) if k]
        dm = obs.is_dmso.to_numpy()
        mu, W = fit_sphere(Zs[dm].astype(np.float64), info)
        V = ((Zs - mu) @ W).astype(np.float32)
        info["n_wells"] = len(obs)
        info["n_dmso_wells"] = int(dm.sum())
        info["n_compound_wells"] = int((~dm).sum())
        info["plates"] = sorted(obs.plate.unique().tolist())
        log(f"[{line} {t}h] plates {info['plates']}; wells {len(obs)} ({info['n_dmso_wells']} DMSO); "
            f"features {info['n_after_nan_const']} after NaN/const (-{info['drop_nan_gt1pct']}/-{info['drop_constant']}) "
            f"-> blocklist -{info['drop_blocklist']}, outlier -{info['drop_outlier_gt500']}, "
            f"low-var -{info['drop_low_variance']}, corr -{info['drop_correlated']} = {info['n_selected']} selected; "
            f"sphere dim {info['sphere_dim']}")
        U = unit_rows(V.astype(np.float64))
        t0 = time.time()
        rep = replicate_map(obs, U, rng)
        log(f"[{line} {t}h] replicate mAP done ({time.time() - t0:.0f}s): "
            f"reproducible {np.mean(rep.rep_p < P_ACTIVE):.3f} of {len(rep)}")
        out[t] = dict(obs=obs, V=V, feats=sel, info=info, rep=rep)
    return out


def consensus(obs, V):
    rows, prof = [], []
    for bs, idx in obs.groupby("broad_sample").indices.items():
        o = obs.iloc[idx]
        prof.append(np.median(V[idx], axis=0))
        rows.append({"broad_sample": bs, "n_wells": len(idx), "n_plates": o.plate.nunique(),
                     "plates": ";".join(sorted(o.plate.unique())), "sites": "source_4",
                     "wells": ";".join(sorted(o.well.unique()))})
    return pd.DataFrame(rows), np.vstack(prof)


# ----------------------------------------------------------------------------- statistics
def cos_mat(P):
    U = unit_rows(P.astype(np.float64))
    return U @ U.T


def spearman(a, b):
    ra, rb = rankdata(a), rankdata(b)
    ra -= ra.mean()
    rb -= rb.mean()
    d = np.sqrt((ra @ ra) * (rb @ rb))
    return float(ra @ rb / d) if d > 0 else np.nan


def pair_index(n, group=None):
    """Upper-triangle pairs, minus pairs within the same group (e.g. two IDs of one unit)."""
    ia, ib = np.triu_indices(n, 1)
    if group is not None:
        g = np.asarray(group)
        m = g[ia] != g[ib]
        ia, ib = ia[m], ib[m]
    return ia, ib


def mantel(SA, SB, group, rng, classes=None, n_perm=N_PERM):
    """Spearman Mantel r; joint row/col permutation of SB on ranks; optional within-class null."""
    n = len(SA)
    ia, ib = pair_index(n, group)
    ra = rankdata(SA[ia, ib])
    ra = (ra - ra.mean()) / ra.std()
    RB = np.zeros((n, n))
    rb = rankdata(SB[ia, ib])
    RB[ia, ib] = rb
    RB[ib, ia] = rb
    r_obs = spearman(SA[ia, ib], SB[ia, ib])

    def perm_r(perm):
        v = RB[perm[ia], perm[ib]]
        return float(np.corrcoef(ra, v)[0, 1])

    nulls = {}
    null = np.array([perm_r(rng.permutation(n)) for _ in range(n_perm)])
    nulls["perm"] = null
    if classes is not None:
        cls = np.asarray(classes)
        blocks = [np.flatnonzero(cls == c) for c in pd.unique(cls)]
        blocks = [b for b in blocks if len(b) > 1]
        nullc = np.empty(n_perm)
        for k in range(n_perm):
            perm = np.arange(n)
            for b in blocks:
                perm[b] = rng.permutation(b)
            nullc[k] = perm_r(perm)
        nulls["within_target"] = nullc
        nulls["n_target_blocks"] = len(blocks)
        from math import factorial, prod
        nulls["n_distinct_perms"] = prod(factorial(len(b)) for b in blocks)
    return r_obs, len(ia), nulls


def pval(null, x):
    return float((1 + np.sum(null >= x)) / (1 + len(null)))


def split_half_mats(obs, V, names, rng, n_splits=N_SPLITS):
    """For each split: per compound, random halves of its replicate wells -> two consensus
    matrices -> two cosine matrices (float32), rows in `names` order."""
    idx_by = obs.groupby("broad_sample").indices
    groups = [np.asarray(idx_by[n]) for n in names]
    H1 = np.empty((n_splits, len(names), len(names)), np.float32)
    H2 = np.empty_like(H1)
    for s in range(n_splits):
        P1, P2 = [], []
        for g in groups:
            g = rng.permutation(g)
            h = len(g) // 2
            P1.append(np.median(V[g[:h]], axis=0))
            P2.append(np.median(V[g[h:]], axis=0))
        H1[s] = cos_mat(np.vstack(P1))
        H2[s] = cos_mat(np.vstack(P2))
    return H1, H2


def sb(r):
    return 2 * r / (1 + r)


def rel_on(H1, H2, sub, ia, ib, splits):
    vals = [spearman(H1[s][sub[ia], sub[ib]], H2[s][sub[ia], sub[ib]]) for s in splits]
    return float(np.mean(vals))


def calibrate_set(tag, t, names, SA, SU, HA, HU, group, classes, rng):
    """All calibration statistics for one compound set at one time point."""
    n = len(names)
    rows = []

    def add(stat, value, lo=np.nan, hi=np.nan, p=np.nan, **kw):
        rows.append({"statistic": stat, "time_h": t, "compound_set": tag, "n_compounds": n,
                     "value": value, "ci_low": lo, "ci_high": hi, "p": p, **kw})

    r, npairs, nulls = mantel(SA, SU, group, rng, classes)
    ia, ib = pair_index(n, group)
    sub = np.arange(n)
    allsplits = range(HA[0].shape[0])
    rel_half = {L: [spearman(H[0][s][ia, ib], H[1][s][ia, ib]) for s in allsplits]
                for L, H in (("A549", HA), ("U2OS", HU))}
    rh = {L: float(np.mean(v)) for L, v in rel_half.items()}
    rsb = {L: sb(v) for L, v in rh.items()}
    att = r / np.sqrt(rsb["A549"] * rsb["U2OS"]) if min(rsb.values()) > 0 else np.nan
    att_half = r / np.sqrt(rh["A549"] * rh["U2OS"]) if min(rh.values()) > 0 else np.nan

    # paired bootstrap over compounds
    grp = np.asarray(group)
    B = {k: [] for k in ("r_cross", "relh_A549", "relh_U2OS", "relsb_A549", "relsb_U2OS",
                         "att", "att_half", "att_A549", "att_U2OS")}
    n_nonpos = 0
    ta, tb = np.triu_indices(n, 1)
    for b in range(N_BOOT):
        bi = rng.integers(0, n, n)
        m = grp[bi[ta]] != grp[bi[tb]]          # drops self-pairs from duplicates + within-unit
        a_, b_ = bi[ta[m]], bi[tb[m]]
        if len(a_) < 3:
            continue
        rc = spearman(SA[a_, b_], SU[a_, b_])
        spl = [(b * SPLITS_PER_BOOT + j) % len(allsplits) for j in range(SPLITS_PER_BOOT)]
        rhb = {L: float(np.mean([spearman(H[0][s][a_, b_], H[1][s][a_, b_]) for s in spl]))
               for L, H in (("A549", HA), ("U2OS", HU))}
        rsbb = {L: sb(v) for L, v in rhb.items()}
        B["r_cross"].append(rc)
        B["relh_A549"].append(rhb["A549"])
        B["relh_U2OS"].append(rhb["U2OS"])
        B["relsb_A549"].append(rsbb["A549"])
        B["relsb_U2OS"].append(rsbb["U2OS"])
        ok = min(rsbb.values()) > 0
        okh = min(rhb.values()) > 0
        n_nonpos += int(not ok)
        B["att"].append(rc / np.sqrt(rsbb["A549"] * rsbb["U2OS"]) if ok else np.nan)
        B["att_half"].append(rc / np.sqrt(rhb["A549"] * rhb["U2OS"]) if okh else np.nan)
        B["att_A549"].append(rc / rsbb["A549"] if rsbb["A549"] > 0 else np.nan)
        B["att_U2OS"].append(rc / rsbb["U2OS"] if rsbb["U2OS"] > 0 else np.nan)

    def ci(k):
        v = np.asarray(B[k], float)
        v = v[np.isfinite(v)]
        return (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) if len(v) > 10 else (np.nan, np.nan)

    nb = len(B["r_cross"])
    add("mantel_r_cross", r, *ci("r_cross"), p=pval(nulls["perm"], r), n_pairs=npairs,
        null_mean=float(nulls["perm"].mean()), n_boot=nb,
        note=f"Spearman Mantel A549 vs U2OS; p from {N_PERM} joint row/col permutations")
    p_wt = pval(nulls["within_target"], r) if nulls["n_distinct_perms"] >= 1000 else np.nan  # untestable
    add("mantel_r_cross_within_target_null", r, *ci("r_cross"), p=p_wt,
        n_pairs=npairs, null_mean=float(nulls["within_target"].mean()), n_boot=nb,
        n_perm_classes=nulls["n_target_blocks"], n_distinct_perms=float(nulls["n_distinct_perms"]),
        note=f"p from {N_PERM} permutations within JUMP-Target `target` classes "
             f"({nulls['n_target_blocks']} classes with >1 member; singletons fixed; "
             f"{nulls['n_distinct_perms']:.3g} distinct permutations)"
             + (f"; NOT TESTABLE: only {nulls['n_distinct_perms']} distinct permutations, so p >= "
                f"{1 / nulls['n_distinct_perms']:.3g} even under perfect agreement"
                if nulls["n_distinct_perms"] < 1000 else ""))
    for L in LINES:
        add(f"rel_half_{L}", rh[L], *ci(f"relh_{L}"), n_pairs=npairs, n_boot=nb,
            note=f"split-half Mantel Spearman within {L}, mean over {len(allsplits)} random splits "
                 f"(sd {np.std(rel_half[L]):.3f})")
        add(f"rel_SB_{L}", rsb[L], *ci(f"relsb_{L}"), n_pairs=npairs, n_boot=nb,
            note="Spearman-Brown corrected to full replicate count: 2r/(1+r)")
    add("attenuation_SB", att, *ci("att"), n_pairs=npairs, n_boot=nb,
        note=f"r_cross / sqrt(rel_SB_A549 * rel_SB_U2OS); not clipped; bootstrap draws with a "
             f"non-positive reliability: {n_nonpos}")
    add("attenuation_halfrel", att_half, *ci("att_half"), n_pairs=npairs, n_boot=nb,
        note="r_cross / sqrt(rel_half_A549 * rel_half_U2OS) (uncorrected reliabilities)")
    for L in LINES:
        add(f"r_cross_over_rel_SB_{L}", r / rsb[L] if rsb[L] > 0 else np.nan, *ci(f"att_{L}"),
            n_pairs=npairs, n_boot=nb, note=f"raw ratio r_cross / rel_SB_{L}")

    # per-compound row agreement
    G = np.asarray(group)
    rho = np.full(n, np.nan)
    for i in range(n):
        m = G != G[i]
        if m.sum() >= 3:
            rho[i] = spearman(SA[i, m], SU[i, m])
    v = rho[np.isfinite(rho)]
    bmed = [np.median(rng.choice(v, len(v))) for _ in range(N_BOOT)]
    add("row_agreement_median", float(np.median(v)), float(np.percentile(bmed, 2.5)),
        float(np.percentile(bmed, 97.5)), note="median per-compound Spearman(A549 row, U2OS row); CI by bootstrap")
    add("row_agreement_mean", float(np.mean(v)))
    add("row_agreement_q25", float(np.percentile(v, 25)))
    add("row_agreement_q75", float(np.percentile(v, 75)))
    add("row_agreement_min", float(v.min()))
    add("row_agreement_max", float(v.max()))
    add("row_agreement_frac_positive", float(np.mean(v > 0)))
    return rows, rho


# ----------------------------------------------------------------------------- diagnostics
def position_only_diag(ds, cons_lt, line, t, rng):
    """(1) DMSO well positions as pseudo-compounds through the replicate-mAP test (their wells
    as replicates, the other DMSO wells of the same plates as negatives, same null): the
    position-only false-positive rate. (2) Norm of the DMSO-position consensus relative to the
    treated consensus norm (positional share of a cosine scales roughly with its square)."""
    o, V = ds["obs"], ds["V"]
    U = unit_rows(V.astype(np.float64))
    dm = o.is_dmso.to_numpy()
    plate = o.plate.to_numpy()
    dmso_by_plate = {p: np.flatnonzero(dm & (plate == p)) for p in np.unique(plate)}
    null_cache, ps, maps, norms = {}, [], [], []
    for w, idx in o[o.is_dmso].groupby("well").indices.items():
        idx = o[o.is_dmso].index.to_numpy()[idx]
        if len(idx) < 2:
            continue
        pl_counts = tuple(sorted(pd.Series(plate[idx]).value_counts().items()))
        pool = np.concatenate([dmso_by_plate[p] for p, _ in pl_counts])
        m = map_of_set(U, idx, np.setdiff1d(pool, idx))
        if pl_counts not in null_cache:
            null = np.empty(N_NULL_MAP)
            for b in range(N_NULL_MAP):
                fake = np.concatenate([rng.choice(dmso_by_plate[p], c, replace=False) for p, c in pl_counts])
                null[b] = map_of_set(U, fake, np.setdiff1d(pool, fake))
            null_cache[pl_counts] = null
        ps.append(pval(null_cache[pl_counts], m))
        maps.append(m)
        norms.append(np.linalg.norm(np.median(V[idx], axis=0)))
    mc, P = cons_lt
    trt = (mc.broad_sample != "DMSO").to_numpy()
    tn = np.median(np.linalg.norm(P[trt], axis=1))
    base = {"time_h": t, "compound_set": f"{line}_dmso_positions", "n_compounds": len(ps),
            "ci_low": np.nan, "ci_high": np.nan, "p": np.nan}
    return [
        {**base, "statistic": f"diag_dmso_position_frac_reproducible_{line}", "value": float(np.mean(np.array(ps) < P_ACTIVE)),
         "note": f"fraction of DMSO well positions (pseudo-compounds) with replicate-mAP p < {P_ACTIVE}; "
                 f"mean mAP {np.mean(maps):.3f}"},
        {**base, "statistic": f"diag_norm_ratio_dmso_position_over_treated_{line}", "value": float(np.median(norms) / tn),
         "note": f"median DMSO-position consensus norm {np.median(norms):.3f} / median treated consensus norm {tn:.3f}"},
    ]


def dmso_position_mantel(dsA, dsU, t, rng):
    """Positional-only baseline: DMSO well positions as pseudo-compounds."""
    P = {}
    for L, ds in (("A549", dsA), ("U2OS", dsU)):
        o, V = ds["obs"], ds["V"]
        d = o[o.is_dmso]
        wells = sorted(d.well.unique())
        P[L] = (wells, np.vstack([np.median(V[d.index[d.well == w]], axis=0) for w in wells]))
    assert P["A549"][0] == P["U2OS"][0]
    n = len(P["A549"][0])
    r, npairs, nulls = mantel(cos_mat(P["A549"][1]), cos_mat(P["U2OS"][1]), np.arange(n), rng)
    return {"statistic": "diag_dmso_position_mantel", "time_h": t, "compound_set": "dmso_positions",
            "n_compounds": n, "value": r, "ci_low": np.nan, "ci_high": np.nan,
            "p": pval(nulls["perm"], r), "n_pairs": npairs, "null_mean": float(nulls["perm"].mean()),
            "note": "cross-line Mantel over DMSO well positions (median over the 4 plates); "
                    "positional-only baseline for r_cross"}


def same_vs_diff_well(ds, line, t):
    """Compounds plated at two positions per plate: replicate cosine for pairs on different
    plates at the same well vs at different wells."""
    o, U = ds["obs"], unit_rows(ds["V"].astype(np.float64))
    same, diff = [], []
    for bs, idx in o[~o.is_dmso].groupby("broad_sample").indices.items():
        idx = o[~o.is_dmso].index.to_numpy()[idx]
        if o.well[idx].nunique() < 2:
            continue
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                i, j = idx[a], idx[b]
                if o.plate[i] == o.plate[j]:
                    continue
                c = float(U[i] @ U[j])
                (same if o.well[i] == o.well[j] else diff).append(c)
    base = {"time_h": t, "compound_set": f"{line}_twice_plated", "ci_low": np.nan, "ci_high": np.nan, "p": np.nan}
    return [
        {**base, "statistic": f"diag_rep_cos_same_well_{line}", "value": float(np.mean(same)),
         "n_pairs": len(same), "note": "mean cosine, replicate pairs on different plates at the same well"},
        {**base, "statistic": f"diag_rep_cos_diff_well_{line}", "value": float(np.mean(diff)),
         "n_pairs": len(diff), "note": "mean cosine, replicate pairs on different plates at different wells"},
    ]


# ----------------------------------------------------------------------------- main
def main():
    keep_cache = "--keep-cache" in sys.argv
    t_start = time.time()
    rng = np.random.default_rng(SEED)
    acc = gi.resolve("cpjump1")
    log(f"accession: {acc}")
    meta = fetch_metadata()
    exp = meta["experiment-metadata.tsv"]
    cm = meta["JUMP-Target-1_compound_metadata.tsv"]
    pm = meta["JUMP-Target-1_compound_platemap.tsv"]
    add = meta["JUMP-Target-1_compound_metadata_additional_annotations.tsv"]
    plates, excluded = select_plates(exp, acc)
    log(f"plates used: {len(plates)}; compound plates excluded by protocol filter: "
        f"{excluded[['Assay_Plate_Barcode', 'Cell_type', 'Time', 'Density', 'Cell_line']].values.tolist()}")

    # compound annotation (keyed by full broad_sample)
    ann = cm[cm.broad_sample.notna()][["broad_sample", "InChIKey", "pert_iname", "target_list",
                                         "pert_type", "control_type"]]
    ann = ann.merge(add[["broad_sample", "target"]], on="broad_sample", how="left")
    assert ann.broad_sample.is_unique and len(ann) == 306, len(ann)
    assert ann.target.notna().all(), "compounds without designed target"
    assert set(pm.broad_sample.dropna()) == set(ann.broad_sample)
    units = pd.read_csv(STEP2 / "units.csv")
    cal = units[units.tier == "CAL_CPJUMP1"]
    id2unit = {i: u for u, ids in zip(cal.unit, cal.dataset_ids) for i in ids.split(";")}
    miss = set(id2unit) - set(ann.broad_sample)
    if miss:
        raise RuntimeError(f"CAL_CPJUMP1 dataset_ids not in JUMP-Target-1: {miss}")
    ann["unit"] = ann.broad_sample.map(id2unit)
    ann["broad_id"] = ann.broad_sample.str[:13]
    assert ann.broad_id.is_unique
    log(f"compounds: {len(ann)} treated; matched CAL_CPJUMP1 units {cal.unit.nunique()} "
        f"-> {ann.unit.notna().sum()} broad_samples")

    man = fetch_plates(acc, plates)
    track_disk()

    # ---- datasets
    ds = {}
    for L in LINES:
        for t, d in process_line(L, plates, pm, ann, rng).items():
            ds[(L, t)] = d
        track_disk()

    # ---- consensus, activity, outputs
    prof_rows, prof_mats, act_rows, well_frames, cons = [], [], [], [], {}
    annk = ann.set_index("broad_sample")
    for (L, t), d in ds.items():
        meta_c, P = consensus(d["obs"], d["V"])
        trt = meta_c.broad_sample != "DMSO"
        norms = np.linalg.norm(P, axis=1)
        nt = norms[trt.to_numpy()]
        med, mad = np.median(nt), MAD_SCALE * np.median(np.abs(nt - np.median(nt)))
        meta_c["norm_z"] = (norms - med) / mad
        meta_c.loc[~trt, "norm_z"] = np.nan
        meta_c = meta_c.merge(d["rep"], on="broad_sample", how="left")
        meta_c["active_rep"] = meta_c.rep_p < P_ACTIVE
        meta_c["active_z"] = meta_c.norm_z > Z_ACTIVE
        # position-robust call: consensus norm beyond the 95th percentile of the consensus norms of
        # the 64 DMSO well positions (DMSO wells of one position across the 4 plates = pure position + noise)
        dmo = d["obs"][d["obs"].is_dmso]
        pos_norms = np.array([np.linalg.norm(np.median(d["V"][dmo.index.to_numpy()[ix]], axis=0))
                              for ix in dmo.groupby("well").indices.values()])
        q95 = float(np.percentile(pos_norms, 95))
        meta_c["active_pos"] = norms > q95
        meta_c.loc[~trt, ["active_rep", "active_z", "active_pos"]] = False
        d["info"]["dmso_position_norm_q50_q95"] = [round(float(np.median(pos_norms)), 3), round(q95, 3)]
        d["info"]["frac_active_pos"] = float(meta_c.loc[trt, "active_pos"].mean())
        meta_c["cell_line"], meta_c["time_h"] = L, t
        meta_c["dose_uM"] = np.where(trt, DOSE_UM, 0.0)
        for c in ["unit", "broad_id", "pert_iname", "target", "control_type", "InChIKey"]:
            meta_c[c] = meta_c.broad_sample.map(annk[c])
        meta_c.loc[~trt, ["broad_id", "pert_iname", "unit"]] = ["DMSO", "DMSO", "DMSO"]
        cons[(L, t)] = (meta_c, P)
        d["info"]["frac_reproducible"] = float(meta_c.loc[trt, "active_rep"].mean())
        d["info"]["frac_active_z"] = float(meta_c.loc[trt, "active_z"].mean())
        d["info"]["frac_reproducible_matched"] = float(meta_c.loc[trt & meta_c.unit.notna(), "active_rep"].mean())
        prof_rows.append(meta_c)
        prof_mats.append(P)
        w = d["obs"][["plate", "well", "time_h", "cell_line", "broad_sample", "is_dmso"]].copy()
        w["unit"] = w.broad_sample.map(annk["unit"])
        w["pert_iname"] = w.broad_sample.map(annk["pert_iname"]).fillna("DMSO")
        w = pd.concat([w.reset_index(drop=True),
                       pd.DataFrame(d["V"], columns=[f"f{i}" for i in range(d["V"].shape[1])])], axis=1)
        well_frames.append(w)

    maxdim = max(P.shape[1] for P in prof_mats)
    cols_meta = ["unit", "broad_id", "broad_sample", "pert_iname", "target", "control_type", "dose_uM",
                 "cell_line", "time_h", "n_wells", "n_plates", "plates", "sites"]
    blocks = []
    for mc, P in zip(prof_rows, prof_mats):
        F = np.full((len(P), maxdim), np.nan, np.float32)
        F[:, :P.shape[1]] = P
        blocks.append(pd.concat([mc[cols_meta].reset_index(drop=True),
                                 pd.DataFrame(F, columns=[f"f{i}" for i in range(maxdim)])], axis=1))
    morph = pd.concat(blocks, ignore_index=True)
    f_out = STEP2 / "morph_cpjump1.parquet"
    morph.to_parquet(f_out, index=False)
    OUT_FILES.append(f_out)
    act = pd.concat(prof_rows, ignore_index=True)[
        ["unit", "broad_id", "broad_sample", "pert_iname", "target", "dose_uM", "cell_line", "time_h",
         "n_wells", "rep_mAP", "rep_p", "null_mAP_mean", "active_rep", "norm_z", "active_z", "active_pos"]]
    f_act = STEP2 / "morph_cpjump1_activity.csv"
    act.to_csv(f_act, index=False)
    OUT_FILES.append(f_act)
    wells = pd.concat(well_frames, ignore_index=True)
    f_w = STEP2 / "morph_cpjump1_wells.parquet"
    wells.to_parquet(f_w, index=False)
    OUT_FILES.append(f_w)
    track_disk()
    log(f"wrote {f_out.name} ({len(morph)} rows), {f_act.name}, {f_w.name} "
        f"({f_w.stat().st_size / 1e6:.1f} MB)")

    # ---- calibration
    cal_rows, row_frames = [], []
    for t in TIMES:
        (mA, PA), (mU, PU) = cons[("A549", t)], cons[("U2OS", t)]
        names = sorted(set(mA.broad_sample[mA.broad_sample != "DMSO"])
                       & set(mU.broad_sample[mU.broad_sample != "DMSO"]))
        ixA = mA.reset_index(drop=True).reset_index().set_index("broad_sample")["index"]
        ixU = mU.reset_index(drop=True).reset_index().set_index("broad_sample")["index"]
        SA_all = cos_mat(PA[ixA[names].to_numpy()])
        SU_all = cos_mat(PU[ixU[names].to_numpy()])
        t0 = time.time()
        HA = split_half_mats(ds[("A549", t)]["obs"], ds[("A549", t)]["V"], names, rng)
        HU = split_half_mats(ds[("U2OS", t)]["obs"], ds[("U2OS", t)]["V"], names, rng)
        log(f"[{t}h] split halves built ({time.time() - t0:.0f}s)")
        actA = mA.set_index("broad_sample").active_rep
        actU = mU.set_index("broad_sample").active_rep
        unit_of = [annk.unit.get(nm) if pd.notna(annk.unit.get(nm)) else nm for nm in names]
        targets = [annk.target[nm] for nm in names]
        sets = {
            "all": np.arange(len(names)),
            "active": np.array([i for i, nm in enumerate(names) if actA[nm] and actU[nm]]),
            "matched": np.array([i for i, nm in enumerate(names) if pd.notna(annk.unit.get(nm))]),
        }
        # strength strata fixed a priori: mean over lines of the consensus-norm percentile
        nA = np.linalg.norm(PA[ixA[names].to_numpy()], axis=1)
        nU = np.linalg.norm(PU[ixU[names].to_numpy()], axis=1)
        strength = (rankdata(nA) + rankdata(nU)) / (2 * len(names))
        lo_q, hi_q = np.quantile(strength, [1 / 3, 2 / 3])
        posA = mA.set_index("broad_sample").active_pos
        posU = mU.set_index("broad_sample").active_pos
        sets["beyond_position"] = np.array([i for i, nm in enumerate(names) if posA[nm] and posU[nm]])
        sets["weak_tercile"] = np.flatnonzero(strength <= lo_q)
        sets["strong_tercile"] = np.flatnonzero(strength >= hi_q)
        for tag, sel in sets.items():
            t0 = time.time()
            sub_names = [names[i] for i in sel]
            rows, rho = calibrate_set(
                tag, t, sub_names, SA_all[np.ix_(sel, sel)], SU_all[np.ix_(sel, sel)],
                (HA[0][:, sel][:, :, sel], HA[1][:, sel][:, :, sel]),
                (HU[0][:, sel][:, :, sel], HU[1][:, sel][:, :, sel]),
                [unit_of[i] for i in sel], [targets[i] for i in sel], rng)
            for rr in rows:
                if tag == "matched":
                    rr["note"] = (rr.get("note", "") + f"; {len(set(unit_of[i] for i in sel))} units / "
                                  f"{len(sel)} broad_samples, within-unit pairs dropped").lstrip("; ")
            cal_rows += rows
            row_frames.append(pd.DataFrame({"broad_sample": sub_names, "time_h": t, "compound_set": tag,
                                            "row_spearman": rho}))
            r0 = rows[0]
            log(f"[{t}h {tag}] n={len(sel)} r_cross={r0['value']:.3f} [{r0['ci_low']:.3f},{r0['ci_high']:.3f}] "
                f"p={r0['p']:.4g}, p_within_target={rows[1]['p']:.4g}; "
                + ", ".join(f"{x['statistic']}={x['value']:.3f}" for x in rows[2:7])
                + f" ({time.time() - t0:.0f}s)")
        cal_rows.append(dmso_position_mantel(ds[("A549", t)], ds[("U2OS", t)], t, rng))
        for L in LINES:
            cal_rows += same_vs_diff_well(ds[(L, t)], L, t)
            cal_rows += position_only_diag(ds[(L, t)], cons[(L, t)], L, t, rng)
            log(f"[{t}h {L}] position-only: " + ", ".join(f"{x['statistic']}={x['value']:.3f}" for x in cal_rows[-2:]))
        del HA, HU
    calib = pd.DataFrame(cal_rows)
    calib["primary"] = calib.time_h == PRIMARY_TIME
    first = ["statistic", "time_h", "compound_set", "n_compounds", "value", "ci_low", "ci_high", "p"]
    calib = calib[first + [c for c in calib.columns if c not in first]]
    f_cal = STEP2 / "cpjump1_calibration.csv"
    calib.to_csv(f_cal, index=False)
    OUT_FILES.append(f_cal)
    rowag = pd.concat(row_frames, ignore_index=True)
    for c in ["pert_iname", "target", "unit"]:
        rowag[c] = rowag.broad_sample.map(annk[c])
    f_row = STEP2 / "cpjump1_row_agreement.csv"
    rowag.to_csv(f_row, index=False)
    OUT_FILES.append(f_row)
    track_disk()

    # ---- run info, cleanup, notes
    n_bytes = gi.Bytes.n
    peak = PEAK["disk"]
    plate_bytes = int(man.get("bytes", pd.Series(dtype=float)).fillna(0).sum()) if len(man) else 0
    if not keep_cache:
        shutil.rmtree(PLATE_CACHE, ignore_errors=True)
    secs = round(time.time() - t_start, 1)
    summary = {"bytes": n_bytes, "peak_disk_bytes": peak, "seconds": secs}
    f_sum = STEP2 / "cpjump1_download_summary.json"
    f_sum.write_text(json.dumps(summary), encoding="utf-8")
    run = {"accession": acc, "batch": BATCH, "meta_files": meta["_info"], "selection_scope": SELECTION_SCOPE,
           "plates": plates.to_dict("records"), "excluded_compound_plates": excluded.to_dict("records"),
           "plate_downloads": man.to_dict("records"), "plate_bytes": plate_bytes,
           "datasets": {f"{L}_{t}h": {**{k: v for k, v in d["info"].items()}, "selected_features": d["feats"]}
                        for (L, t), d in ds.items()},
           "summary": summary, "log": LOG}
    (STEP2 / "cpjump1_run.json").write_text(json.dumps(run, indent=1, default=str), encoding="utf-8")
    write_notes(acc, meta, plates, excluded, ds, calib, man, summary, keep_cache, ann)
    log(f"done in {secs}s; bytes {n_bytes / 1e6:.1f} MB; peak disk {peak / 1e6:.1f} MB")


def fmt_row(c, stat, t, tag):
    r = c[(c.statistic == stat) & (c.time_h == t) & (c.compound_set == tag)]
    if not len(r):
        return "n/a"
    r = r.iloc[0]
    s = f"{r.value:.3f}"
    if np.isfinite(r.ci_low):
        s += f" [{r.ci_low:.3f}, {r.ci_high:.3f}]"
    if np.isfinite(r.p):
        s += f" (p<={r.p:.0e})" if r.p <= 1 / (N_PERM + 1) + 1e-12 else f" (p={r.p:.2g})"
    return s


def write_notes(acc, meta, plates, excluded, ds, calib, man, summary, keep_cache, ann):
    import fsspec, pyarrow, scipy, sklearn  # noqa: E401  (versions only)
    L = ["# CPJUMP1 calibration notes (written by step2/morph_cpjump1.py)", ""]
    L += ["## Sources",
          f"- Gallery accession `{acc}` (resolved with `gallery_index.resolve('cpjump1')`), batch `{BATCH}`.",
          f"- Profiles: `{acc}/source_4/workspace/profiles/{BATCH}/<plate>/<plate>.csv.gz` (verified; raw "
          "well-level CellProfiler aggregates, only `Metadata_Plate`/`Metadata_Well` as metadata, 384 wells, "
          "5,792 features incl. brightfield channels). Other files in the folder (`_augmented`, `_normalized*`, "
          "`_spherized`) are not used.",
          "- Metadata (cached in `step2/cpjump1_cache/`, sha256 prefix): "
          + "; ".join(f"`{k}` {v}" for k, v in meta["_info"].items() if k.endswith("_sha256")),
          "  - experiment metadata and additional annotations: GitHub jump-cellpainting/2024_Chandrasekaran_NatureMethods "
          "(`benchmark/output/experiment-metadata.tsv`, `benchmark/input/JUMP-Target-1_compound_metadata_additional_annotations.tsv`, main);",
          "  - compound metadata and platemap: GitHub jump-cellpainting/JUMP-Target (master).",
          "- `target` (single designed target gene per compound; 146 genes x 2 compounds + 14 singletons) comes from the "
          "NatureMethods additional annotations: the JUMP-Target compound_metadata.tsv itself only has `target_list` "
          "(multi-gene). The compound key for everything is the full `broad_sample`; units map through `units.csv` "
          "`dataset_ids` (not InChIKey: the two metadata files disagree on some InChIKeys, e.g. amlodipine).",
          ""]
    L += ["## Plates (one protocol)",
          "- Filter = the NatureMethods benchmark's own: batch == 2020_11_04_CPJUMP1 (not `_DL`), Perturbation == compound, "
          "Density == 100, Antibiotics absent, Cas9 compound plates dropped. ORF/CRISPR plates excluded.",
          ]
    for (cl, t), g in plates.groupby(["cell_line", "time_h"]):
        L.append(f"- {cl} {t} h: {', '.join(g.plate)} (anomaly flags: {', '.join(g.Anomaly.astype(str))}; "
                 f"images: {', '.join(g.Number_of_images.astype(str))})")
    L.append("- Excluded compound plates of the batch (density/Cas9): "
             + ", ".join(f"{r.Assay_Plate_Barcode} ({r.Cell_type} {r.Time}h, density {r.Density}, {r.Cell_line})"
                         for r in excluded.itertuples()))
    L += ["- Caveats: all 24 h plates carry the `WGA` anomaly flag (both lines, balanced); 2 of 4 U2OS 48 h plates carry "
          "`Mitotracker`; the A549 48 h plates were imaged with 16 instead of 9 sites per well (~49k images). "
          "Kept as published by the benchmark.",
          "- Every plate uses the same JUMP-Target-1 platemap (306 compounds, 292 once and 14 twice per plate, "
          "64 DMSO wells), so each compound sits at the same well in all replicates and in both cell lines. "
          "Positional effects can therefore inflate replicate mAP, the split-half ceiling and r_cross alike; "
          "see the diagnostics below.", ""]
    L += ["## Processing (per cell line x time point, against its own DMSO)",
          f"- Raw columns: one intersection per cell line across both time points. Feature-selection scope: "
          f"`{SELECTION_SCOPE}` (constant `SELECTION_SCOPE`; 'line' would select once per cell line). Sphering is "
          "always per line x time.",
          f"- Pre-filter (additions to the spec): brightfield-channel features dropped (Brightfield/HighZBF/LowZBF, "
          f"to keep the standard 5-channel Cell Painting space of the other tiers) and image-coordinate features "
          f"(AreaShape_Center_*, BoundingBoxMin/Max, Location_*).",
          f"- NaN/inf > 1% or constant dropped; NaN -> plate DMSO median; robust z vs plate DMSO with MAD floor "
          f"(0.01 x median plate scale); clip +-50; pycytominer 55-feature blocklist; outlier features (pre-clip "
          f"max|z| > {OUTLIER_Z:g}); variance threshold (freq_cut {FREQ_CUT}, unique_cut {UNIQUE_CUT}) evaluated on the "
          f"UNCLIPPED z (on clipped z, strong-responder features where many treated wells sit at +-50 were flagged "
          f"as low-variance because the clip value became the modal value: 426 features in U2OS 48 h, 187 in U2OS 24 h); "
          f"greedy |r| > "
          f"{CORR_THRESH} pruning.",
          f"- Sphering: PCA whitening fit on DMSO wells only, 99% DMSO variance, eigenvalue + {SPHERE_EPS} "
          "(fixed a priori, same as the JUMP tier; not tuned on any replicate or cross-line statistic).",
          "- Consensus: median of replicate wells (single site, source_4). DMSO reference row = median of all DMSO wells.",
          f"- Replicate mAP: each replicate queries its other replicates + all DMSO wells on the same plates by cosine "
          f"(sphered); null = {N_NULL_MAP} random DMSO sets with the same per-plate composition (queried against the "
          f"remaining DMSO wells); p = (1+#null>=obs)/(1+N); reproducible / active_rep if p < {P_ACTIVE}. "
          f"active_z: robust z of the consensus L2 norm vs all treated compounds > {Z_ACTIVE:g}.", ""]
    L += ["| dataset | plates | wells (DMSO) | raw int. | after prefilter | NaN/const drop | after | blocklist | outlier | low var | corr | selected | sphere dim | reproducible | active_z | active_pos | reproducible (matched) |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for (cl, t), d in ds.items():
        i = d["info"]
        L.append(f"| {cl} {t} h | {len(i['plates'])} | {i['n_wells']} ({i['n_dmso_wells']}) | {i['n_raw_intersection']} | "
                 f"{i['n_after_prefilter']} (-{i['drop_brightfield']} bf, -{i['drop_positional_coords']} pos) | "
                 f"-{i['drop_nan_gt1pct']}/-{i['drop_constant']} | {i['n_after_nan_const']} | -{i['drop_blocklist']} | "
                 f"-{i['drop_outlier_gt500']} | -{i['drop_low_variance']} | -{i['drop_correlated']} | {i['n_selected']} | "
                 f"{i['sphere_dim']} | {i['frac_reproducible']:.3f} | {i['frac_active_z']:.3f} | {i['frac_active_pos']:.3f} | {i['frac_reproducible_matched']:.3f} |")
    L += ["", f"- Imputed values / floored plate x feature MAD fraction: " + "; ".join(
        f"{cl} {t}h {d['info']['n_values_imputed']} / {d['info']['plate_feature_scales_floored_frac']}"
        for (cl, t), d in ds.items()),
          f"- Compounds: {len(ann)} treated (all 5 uM), 4 wells per line x time (8 for the 14 twice-plated), except where "
          f"wells are absent from the gallery profile files: "
          + "; ".join(f"{p}: {len(w)} wells ({', '.join(w[:6])}{', ...' if len(w) > 6 else ''})"
                      for (cl, t), d in ds.items() for p, w in d['info']['missing_wells'].items()
                      if p in d['info']['plates'])
          + ". BR00116995 (U2OS 24 h) is a partial plate (83 of 384 wells missing in odd rows A-O, cols 1-14); kept, so "
          f"those compounds have 3 U2OS 24 h replicates. "
          f"Matched subset: {ann.unit.nunique()} CAL_CPJUMP1 units = {ann.unit.notna().sum()} broad_samples "
          "(dexamethasone and quinidine/quinine map 2 IDs to one unit; within-unit pairs are dropped from the matched "
          "Mantel vectors and bootstrap).", ""]
    L += ["## Calibration (cpjump1_calibration.csv; 48 h primary)",
          f"- Mantel: Spearman over the upper triangle of the compound x compound cosine matrices of the consensus "
          f"profiles (DMSO excluded; compounds present in both lines). p from {N_PERM} joint row/column permutations "
          f"of the U2OS matrix (on ranks), and from {N_PERM} permutations within JUMP-Target `target` classes "
          f"(singletons fixed). CI: {N_BOOT} paired bootstrap draws over compounds (self-pairs from duplicated "
          "compounds dropped); the same draw gives r_cross, both reliabilities and the ratios.",
          f"- Reliability: {N_SPLITS} random split halves per line x time (each compound's wells split independently, "
          f"2+2 or 4+4), Mantel Spearman between the half matrices; mean over splits; Spearman-Brown 2r/(1+r). "
          f"Within each bootstrap draw {SPLITS_PER_BOOT} rotating splits are used. Reliabilities are recomputed on "
          "each compound set.",
          "- Active set: replicate-mAP p < 0.05 in both lines at that time point. CAUTION: this call is position-"
          "confounded (all replicates share one well; 36-70% of DMSO well positions pass the same test, see diagnostics), "
          "so `active` is close to `all`.",
          "- `beyond_position` set: consensus norm > 95th percentile of the 64 DMSO-position consensus norms in BOTH "
          "lines (`active_pos` in the activity file); a position-robust activity call.", ""]
    L += ["| time | set | n | r_cross [95% CI] (perm p) | p within-target | rel_SB A549 | rel_SB U2OS | attenuation_SB | attenuation (half rel) | row agreement median |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for t in TIMES:
        for tag in ("all", "active", "beyond_position", "strong_tercile", "weak_tercile", "matched"):
            r = calib[(calib.statistic == "mantel_r_cross") & (calib.time_h == t) & (calib.compound_set == tag)]
            if not len(r):
                continue
            wr = calib[(calib.statistic == "mantel_r_cross_within_target_null") & (calib.time_h == t)
                       & (calib.compound_set == tag)].iloc[0]
            if wr.n_distinct_perms < 1000:
                pw = f"not testable ({int(wr.n_perm_classes)} swappable classes, min p {1 / wr.n_distinct_perms:.2g})"
            elif wr.p <= 1 / (N_PERM + 1) + 1e-12:
                pw = f"<={wr.p:.0e}"
            else:
                pw = f"{wr.p:.2g}"
            L.append(f"| {t} h | {tag} | {int(r.n_compounds.iloc[0])} | {fmt_row(calib, 'mantel_r_cross', t, tag)} | "
                     f"{pw} | {fmt_row(calib, 'rel_SB_A549', t, tag)} | {fmt_row(calib, 'rel_SB_U2OS', t, tag)} | "
                     f"{fmt_row(calib, 'attenuation_SB', t, tag)} | {fmt_row(calib, 'attenuation_halfrel', t, tag)} | "
                     f"{fmt_row(calib, 'row_agreement_median', t, tag)} |")
    L += ["", "Diagnostics (positional confound):"]
    for t in TIMES:
        L.append(f"- {t} h: DMSO-position cross-line Mantel {fmt_row(calib, 'diag_dmso_position_mantel', t, 'dmso_positions')}; "
                 + "; ".join(f"{cl} twice-plated replicate cosine same well "
                             f"{fmt_row(calib, f'diag_rep_cos_same_well_{cl}', t, f'{cl}_twice_plated')} vs different well "
                             f"{fmt_row(calib, f'diag_rep_cos_diff_well_{cl}', t, f'{cl}_twice_plated')}" for cl in LINES)
                 + "; " + "; ".join(
                     f"{cl} DMSO positions passing the replicate test "
                     f"{fmt_row(calib, f'diag_dmso_position_frac_reproducible_{cl}', t, f'{cl}_dmso_positions')}, "
                     f"DMSO-position/treated norm ratio "
                     f"{fmt_row(calib, f'diag_norm_ratio_dmso_position_over_treated_{cl}', t, f'{cl}_dmso_positions')}"
                     for cl in LINES))
    a = pd.read_csv(STEP2 / "morph_cpjump1_activity.csv")
    a = a[a.broad_sample != "DMSO"]
    for t in TIMES:
        w = a[a.time_h == t].pivot_table(index="broad_sample", columns="cell_line", values="active_pos",
                                         aggfunc="first").astype(bool)
        A, U = w["A549"], w["U2OS"]
        L.append(f"- {t} h `active_pos` (beyond position): both lines {(A & U).sum()}, A549 only {(A & ~U).sum()}, "
                 f"U2OS only {(~A & U).sum()}, neither {(~A & ~U).sum()}. Sets requiring activity in both lines "
                 "(beyond_position, strong_tercile) exclude line-specific activity, itself a cell-line effect, so they "
                 "bracket the attenuation from above; the A549 48 h plates were imaged at 16 sites/well, so part of the "
                 "asymmetry may be sensitivity.")
    L += ["- Strength strata (`weak_tercile`, `strong_tercile`): compounds ranked by the mean over lines of their "
          "consensus-norm percentile, fixed before any calibration statistic was looked at in those strata. If r_cross "
          "in the weak tercile drifts toward the DMSO-position Mantel, the all-compound number is position-contaminated."]
    L += ["", "## Outputs",
          "- `step2/morph_cpjump1.parquet`: consensus per broad_sample x cell line x time (all 306 compounds + a DMSO "
          "reference row), keyed by `broad_id` (broad_sample prefix) and `unit` where matched; f0..fN float32, "
          "NaN beyond a dataset's dimensionality. f-columns are NOT comparable across the four datasets.",
          "- `step2/morph_cpjump1_activity.csv`, `step2/morph_cpjump1_wells.parquet` (sphered wells), "
          "`step2/cpjump1_calibration.csv`, `step2/cpjump1_row_agreement.csv` (per-compound row Spearman), "
          "`step2/cpjump1_run.json` (selected feature names, per-plate download log, full log).", ""]
    L += ["## Bytes, disk, time",
          f"- Downloaded this run: {summary['bytes'] / 1e6:.1f} MB (plates {pd.to_numeric(man.get('bytes', pd.Series(dtype=float)), errors='coerce').fillna(0).sum() / 1e6:.1f} MB newly downloaded, "
          "plus S3 listings and GitHub metadata). Development runs before this one: ~128 MB (a first run that downloaded all 16 plate files, ~113 MB, "
          "before stopping at a missing-well assertion; two probe downloads of one plate file, ~14 MB; listings and "
          "GitHub files). The final run started from an empty `step2/cpjump1_cache/`.",
          f"- Peak disk (cache + outputs): {summary['peak_disk_bytes'] / 1e6:.1f} MB. Plate cache "
          f"`step2/cpjump1_cache/plates/` {'kept (--keep-cache)' if keep_cache else 'deleted at the end'}; "
          "metadata cache kept.",
          f"- Wall time: {summary['seconds']:.0f} s.",
          f"- Packages: python {sys.version.split()[0]}, pandas {pd.__version__}, numpy {np.__version__}, scipy {scipy.__version__}, "
          f"scikit-learn {sklearn.__version__}, pyarrow {pyarrow.__version__}, fsspec {fsspec.__version__}, requests {requests.__version__}.",
          "", "## Deviations from MORPH_SPEC",
          "- Outputs cover ALL 306 compounds, not only matched units (task-mandated: the calibration uses all compounds).",
          "- Extra pre-filter: brightfield features and image-coordinate features dropped (counts above).",
          "- Feature space: normalization/selection/sphering per line x time; the 'shared within a cell line' requirement "
          "is met at the raw-column level (one intersection per line). No statistic crosses time points.",
          "- Variance threshold on unclipped z (see above); the JUMP tier ran it on clipped z.",
          "- Plates with missing wells (BR00116995 partial, BR00116993 one well) kept; DMSO-count check relaxed.",
          "- Within-unit pairs (dexamethasone x2, quinidine/quinine) are dropped from every Mantel vector, not only the "
          "matched set (2 of ~46k pairs in the all set).",
          "- Single site/source, so the consensus is a plain median."]
    f = STEP2 / "cpjump1_notes.md"
    f.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
