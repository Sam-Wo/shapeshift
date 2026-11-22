"""T1 morphology profiles: LINCS Cell Painting (cpg0004), A549, 48 h, 6-dose compound screen.

Run:  py -3.13 -I -X utf8 step2/morph_lincs.py [--plate-cache DIR]

Follows step2/MORPH_SPEC.md plus the coordinator addendum on activity calls:
  1. Per-plate well-level *_augmented.csv.gz (raw aggregated CellProfiler + platemap metadata)
     streamed from the gallery into memory (no raw file is written unless --plate-cache is given,
     a development option). Metadata_* dropped, plus non-morphological bookkeeping columns
     (Location_*, Number_Object_Number, Parent_*, Neighbors_*ObjectNumber*). Features NaN/inf in
     > 1% of wells or constant are dropped; remaining NaNs imputed with the plate's DMSO median.
  2. Per-plate robust z vs DMSO: (x - med) / max(1.4826 MAD, floor), clip +-50.
  3. Selection (whole batch): pycytominer blocklist, outlier (pre-clip max|z| > 500),
     variance threshold, greedy |r| > 0.9 pruning.
  4. PCA whitening fitted on DMSO wells only, 99% variance, eigenvalue + eps (eps fixed = 3.0,
     the JUMP T3 value; not tuned here).
  5. Consensus per unit x nominal dose: median over replicate wells (single site).
  6. Replicate mAP (cosine) vs DMSO on the same plates; null = plate-matched random DMSO sets.
  7. Activity: active_rep (6), active_z (robust z of consensus norm vs all treated compounds at
     that dose), active_norm_dmso (consensus norm vs DMSO pseudo-consensus norms),
     active_vs_cpd (replicate mAP vs sets of other compounds' wells from the same plates, same dose).
"""
from __future__ import annotations

import argparse
import io
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from importlib.metadata import version
from pathlib import Path

STEP2 = Path(__file__).resolve().parent
sys.path.insert(0, str(STEP2))  # -I drops the script dir from sys.path

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import requests  # noqa: E402

from gallery_index import Bytes, resolve, s3_list, s3_url  # noqa: E402

ROOT = STEP2.parent
CACHE = STEP2 / "lincs_cache"
BATCH = "2016_04_01_a549_48hr_batch1"
OUT_CONS = STEP2 / "morph_lincs.parquet"
OUT_ACT = STEP2 / "morph_lincs_activity.csv"
OUT_WELLS = STEP2 / "morph_lincs_wells.parquet"
OUT_LOW = STEP2 / "morph_lincs_lowest_active.csv"
OUT_LIB = STEP2 / "morph_lincs_library_activity.csv"
OUT_NOTES = STEP2 / "lincs_notes.md"
OUT_DL = STEP2 / "lincs_download_summary.json"
TIER = "T1_LINCS_A549"
SITE = "broad"

NOMINAL_UM = np.array([0.0412, 0.123, 0.37, 1.11, 3.33, 10.0])  # units.csv labels
DOSE_TOL = 1.5            # max ratio actual/nominal (either way) to assign a nominal dose
NAN_FRAC_MAX = 0.01
MAD_SCALE = 1.4826
MAD_FLOOR_REL = 0.01      # floor = this x median-over-plates of the plate DMSO scale
MAD_FLOOR_ABS = 1e-12
OUTLIER_Z = 500.0
CLIP_Z = 50.0
FREQ_CUT, UNIQUE_CUT = 0.05, 0.01
CORR_THRESH = 0.9
VAR_KEEP = 0.99
SPHERE_EPS = 3.0
N_NULL = 2000
P_ACTIVE = 0.05
Z_ACTIVE = 3.0
SEED = 0
N_THREADS = 8
CHANNELS = ["AGP", "DNA", "ER", "Mito", "RNA"]
BLOCKLIST = set(
    [f"Nuclei_Correlation_Manders_{a}_{b}" for a in CHANNELS for b in CHANNELS if a != b]
    + [f"Nuclei_Correlation_RWC_{a}_{b}" for a in CHANNELS for b in CHANNELS if a != b]
    + [f"Nuclei_Granularity_{g}_{c}" for g in (14, 15, 16) for c in CHANNELS])
assert len(BLOCKLIST) == 55

LOG: list[str] = []
ASSEMBLY: dict = {}
DISK = {"peak": 0}
TRACKED: set[Path] = set()


def log(msg=""):
    print(msg, flush=True)
    LOG.append(str(msg))


def dir_bytes(p: Path) -> int:
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file()) if p.exists() else 0


def disk_check(extra_dirs=()):
    n = dir_bytes(CACHE) + sum(dir_bytes(Path(d)) for d in extra_dirs if d)
    n += sum(p.stat().st_size for p in TRACKED if p.exists())
    DISK["peak"] = max(DISK["peak"], n)
    if n > 1.5e9:
        raise RuntimeError(f"disk budget exceeded: {n / 1e9:.2f} GB")
    return n


def get(url: str, tries: int = 5) -> bytes:
    last = None
    for k in range(tries):
        try:
            r = requests.get(url, timeout=180)
            if r.status_code == 200:
                Bytes.add(len(r.content))
                return r.content
            last = RuntimeError(f"HTTP {r.status_code} for {url}")
            if r.status_code in (403, 404):
                break
        except requests.RequestException as e:
            last = e
        time.sleep(2 ** k)
    raise RuntimeError(f"source unreachable: {url}: {last}")


# ----------------------------------------------------------------------------- metadata
def fetch_metadata(acc: str):
    """Cache barcode_platemap + platemaps under lincs_cache/; return (barcodes, platemap)."""
    base = f"{acc}/broad/workspace/metadata/platemaps/{BATCH}/"
    (CACHE / "platemap").mkdir(parents=True, exist_ok=True)
    bp_path = CACHE / "barcode_platemap.csv"
    if not bp_path.exists():
        bp_path.write_bytes(get(s3_url(base + "barcode_platemap.csv")))
    _, keys = s3_list(base + "platemap/")
    assert keys, "no platemaps listed"
    for k, _ in keys:
        p = CACHE / "platemap" / Path(k).name
        if not p.exists():
            p.write_bytes(get(s3_url(k)))
    bp = pd.read_csv(bp_path)
    pm = pd.concat([pd.read_csv(f, sep="\t") for f in sorted((CACHE / "platemap").glob("*.txt"))],
                   ignore_index=True)
    return bp, pm


def list_plates(acc: str) -> dict:
    prefixes, _ = s3_list(f"{acc}/broad/workspace/profiles/{BATCH}/")
    plates = {}
    for p in prefixes:
        plate = p.rstrip("/").split("/")[-1]
        plates[plate] = f"{p}{plate}_augmented.csv.gz"
    return plates


def load_units():
    u = pd.read_csv(STEP2 / "units.csv")
    u = u[u.tier == TIER]
    ids = u.assign(bs=u.dataset_ids.str.split(";")).explode("bs")[["unit", "bs"]]
    ids["bs"] = ids.bs.str.strip()
    return u, ids


# ----------------------------------------------------------------------------- download
def fetch_plate(plate: str, key: str, plate_cache: Path | None):
    if plate_cache is not None:
        f = plate_cache / f"{plate}.parquet"
        if f.exists():
            return plate, pd.read_parquet(f), 0
    raw = get(s3_url(key))
    df = pd.read_csv(io.BytesIO(raw), compression="gzip", low_memory=False)
    meta = [c for c in df.columns if c.startswith("Metadata_")]
    feat = [c for c in df.columns if not c.startswith("Metadata_")]
    df = pd.concat([df[meta], df[feat].apply(pd.to_numeric, errors="coerce").astype(np.float32)], axis=1)
    if plate_cache is not None:
        df.to_parquet(plate_cache / f"{plate}.parquet", index=False)
    return plate, df, len(raw)


def download(plates: dict, plate_cache):
    t0 = time.time()
    frames, nbytes = {}, {}
    with ThreadPoolExecutor(N_THREADS) as ex:
        futs = [ex.submit(fetch_plate, p, k, plate_cache) for p, k in sorted(plates.items())]
        for i, f in enumerate(futs, 1):
            p, df, nb = f.result()  # raises loudly on failure
            frames[p], nbytes[p] = df, nb
            if i % 20 == 0 or i == len(futs):
                print(f"  {i}/{len(futs)} plates, {sum(nbytes.values()) / 1e6:.0f} MB, "
                      f"{time.time() - t0:.0f}s", flush=True)
                disk_check([plate_cache])
    return frames, nbytes, time.time() - t0


def assemble(frames: dict, pm: pd.DataFrame):
    featsets = {p: [c for c in d.columns if not c.startswith("Metadata_")] for p, d in frames.items()}
    common = set.intersection(*(set(f) for f in featsets.values()))
    union = sorted(set().union(*featsets.values()))
    first = [c for c in featsets[sorted(featsets)[0]] if c in common]
    lacking = {p: len(set(union) - set(f)) for p, f in featsets.items() if set(f) != set(union)}
    ASSEMBLY.update(n_features_union=len(union), n_features_intersection=len(first),
                    plates_lacking_features=lacking,
                    features_not_on_all_plates=sorted(set(union) - common))
    log(f"feature columns: union {len(union)}, intersection {len(first)}; "
        f"{len(lacking)} plates lack {sorted(set(lacking.values()))} features")
    obs_cols = ["Metadata_Plate", "Metadata_Well", "Metadata_Plate_Map_Name", "Metadata_broad_sample",
                "Metadata_mmoles_per_liter", "Metadata_pert_type", "Metadata_pert_iname"]
    obs = pd.concat([d[obs_cols] for _, d in sorted(frames.items())], ignore_index=True)
    X = np.concatenate([d[first].to_numpy(np.float32) for _, d in sorted(frames.items())])
    obs.columns = ["plate", "well", "plate_map", "broad_sample", "dose_actual_uM", "pert_type", "pert_iname"]
    assert not obs.duplicated(["plate", "well"]).any()
    # cross-check against the cached platemaps
    chk = obs.merge(pm.rename(columns={"plate_map_name": "plate_map", "well_position": "well"}),
                    on=["plate_map", "well"], how="left", validate="m:1")
    bs_pm = chk.broad_sample_y.fillna("DMSO")
    assert (bs_pm == chk.broad_sample_x).all(), "broad_sample mismatch vs platemap"
    dpm = chk.mmoles_per_liter.fillna(0.0).to_numpy()
    assert np.allclose(dpm, chk.dose_actual_uM.to_numpy(), rtol=1e-3, atol=1e-6), "dose mismatch vs platemap"
    obs["is_dmso"] = obs.broad_sample == "DMSO"
    assert (obs.loc[obs.is_dmso, "pert_type"] == "control").all()
    obs["b13"] = np.where(obs.is_dmso, "DMSO", obs.broad_sample.str[:13])
    return obs, X, first


def nominal_dose(actual: np.ndarray):
    a = np.asarray(actual, float)
    out = np.full(len(a), np.nan)
    pos = a > 0
    la = np.log(a[pos])[:, None]
    d = np.abs(la - np.log(NOMINAL_UM)[None, :])
    j = d.argmin(1)
    ok = d[np.arange(len(j)), j] <= np.log(DOSE_TOL)
    vals = NOMINAL_UM[j]
    vals[~ok] = np.nan
    out[pos] = vals
    return out


# ----------------------------------------------------------------------------- processing
NONMORPH = ("_Location_", "_Number_Object_Number", "_Parent_", "ObjectNumber")


def normalize(obs, X, feats, info):
    f = np.array(feats)
    nonmorph = np.array([any(s in c for s in NONMORPH) for c in feats])
    info["n_features_raw"] = len(feats)
    info["drop_nonmorph_bookkeeping"] = int(nonmorph.sum())
    X = X[:, ~nonmorph]
    f = f[~nonmorph]
    finite = np.isfinite(X)
    nan_frac = 1 - finite.mean(0)
    X = np.where(finite, X, np.nan).astype(np.float32)
    std = np.nanstd(X, axis=0)
    drop_nan = nan_frac > NAN_FRAC_MAX
    drop_const = (~drop_nan) & ~(std > 0)
    keep = ~(drop_nan | drop_const)
    info["drop_nan_gt1pct"] = int(drop_nan.sum())
    info["drop_constant"] = int(drop_const.sum())
    X, f = X[:, keep], f[keep]
    info["n_after_nan_const"] = len(f)
    plates = obs.plate.to_numpy()
    dm = obs.is_dmso.to_numpy()
    groups = [np.sort(ix) for ix in obs.groupby("plate").indices.values()]
    meds = np.empty((len(groups), X.shape[1]))
    scales = np.empty_like(meds)
    for g, idx in enumerate(groups):
        D = X[idx[dm[idx]]].astype(np.float64)
        assert len(D) >= 8, f"too few DMSO wells on {plates[idx[0]]}"
        meds[g] = np.nanmedian(D, 0)
        scales[g] = MAD_SCALE * np.nanmedian(np.abs(D - meds[g]), 0)
    ref = np.nanmedian(scales, 0)
    floor = np.maximum(MAD_FLOOR_REL * np.nan_to_num(ref), MAD_FLOOR_ABS)
    floored = ~(scales >= floor[None])
    info["plate_feature_scales_floored_frac"] = float(floored.mean())
    scales = np.where(floored, floor[None], scales)
    Z = np.empty_like(X)
    n_imp = 0
    for g, idx in enumerate(groups):
        x = X[idx].astype(np.float64)
        nan = np.isnan(x)
        n_imp += int(nan.sum())
        x = np.where(nan, meds[g][None], x)
        z = (x - meds[g][None]) / scales[g][None]
        z[~np.isfinite(z)] = 0.0
        Z[idx] = z
    info["n_values_imputed"] = n_imp
    maxabs = np.abs(Z).max(0)
    np.clip(Z, -CLIP_Z, CLIP_Z, out=Z)
    return Z, maxabs, f


def variance_threshold_drop(Z):
    n = Z.shape[0]
    drop = np.zeros(Z.shape[1], bool)
    for j in range(Z.shape[1]):
        _, cnt = np.unique(Z[:, j], return_counts=True)
        if len(cnt) < 2:
            drop[j] = True
            continue
        top2 = np.partition(cnt, len(cnt) - 2)[-2:]
        drop[j] = (top2[0] / top2[1] < FREQ_CUT) or (len(cnt) / n < UNIQUE_CUT)
    return drop


def correlation_prune(Z):
    mu = Z.mean(0, dtype=np.float64)
    sd = Z.std(0, dtype=np.float64)
    S = ((Z - mu) / sd).astype(np.float32)
    R = np.abs((S.T @ S) / len(S)).astype(np.float64)
    R[~np.isfinite(R)] = 0
    order = np.argsort(R.sum(0), kind="stable")
    keep = np.zeros(Z.shape[1], bool)
    kept = []
    for j in order:
        if not kept or R[j, kept].max() <= CORR_THRESH:
            keep[j] = True
            kept.append(j)
    return keep


def select_features(Z, maxabs, f, info):
    m_block = np.isin(f, list(BLOCKLIST))
    info["drop_blocklist"] = int(m_block.sum())
    m_out = (~m_block) & (maxabs > OUTLIER_Z)
    info["drop_outlier_gt500"] = int(m_out.sum())
    keep = ~(m_block | m_out)
    vd = np.zeros(len(f), bool)
    vd[keep] = variance_threshold_drop(Z[:, keep])
    info["drop_low_variance"] = int(vd.sum())
    keep &= ~vd
    ck = np.zeros(len(f), bool)
    ck[keep] = correlation_prune(Z[:, keep])
    info["drop_correlated"] = int((keep & ~ck).sum())
    keep &= ck
    info["n_selected"] = int(keep.sum())
    return Z[:, keep], f[keep].tolist()


def fit_sphere(D, info):
    mu = D.mean(0)
    lam, U = np.linalg.eigh(np.cov(D - mu, rowvar=False))
    lam, U = lam[::-1].clip(min=0), U[:, ::-1]
    k = int(np.searchsorted(np.cumsum(lam) / lam.sum(), VAR_KEEP) + 1)
    info["sphere_dim"] = k
    info["sphere_eig_first_last"] = [round(float(lam[0]), 3), round(float(lam[k - 1]), 4)]
    info["sphere_n_eig_below_eps"] = int((lam[:k] < SPHERE_EPS).sum())
    return mu, U[:, :k] / np.sqrt(lam[:k] + SPHERE_EPS)


def unitv(V):
    n = np.linalg.norm(V, axis=-1, keepdims=True)
    return V / np.where(n > 0, n, 1)


# ----------------------------------------------------------------------------- mAP machinery
def batch_map(G, Q, neg, neg_valid=None):
    """Mean average precision for sets of queries.
    G: (M, M) cosine matrix of a local pool; Q: (D, n) indices into G (each row one replicate
    set); neg: (N,) indices of negatives; neg_valid: (D, N) bool or None.
    Each replicate queries {other replicates} + negatives ranked by cosine. Returns (D,) mAP."""
    D, n = Q.shape
    pos = G[Q[:, :, None], Q[:, None, :]]                        # (D, n, n)
    off = ~np.eye(n, dtype=bool)
    pos = pos[:, off].reshape(D, n, n - 1)
    pos = -np.sort(-pos, axis=-1)                                # descending
    ng = G[Q[:, :, None], neg[None, None, :]]                    # (D, n, N)
    if neg_valid is not None:
        ng = np.where(neg_valid[:, None, :], ng, -np.inf)
    above = (ng[:, :, None, :] > pos[:, :, :, None]).sum(-1)     # (D, n, n-1)
    k = np.arange(1, n, dtype=float)
    ap = (k / (k + above)).mean(-1)
    return ap.mean(-1)


def draw_slots(rng, slot_cands: list[np.ndarray], n_draw: int):
    """One candidate per slot, without replacement among slots that share a candidate list
    (slots on the same plate share a list). Returns (n_draw, n_slots) indices."""
    out = np.empty((n_draw, len(slot_cands)), np.int64)
    groups = {}
    for s, c in enumerate(slot_cands):
        groups.setdefault(id(c), (c, []))[1].append(s)
    for c, slots in groups.values():
        keys = rng.random((n_draw, len(c)))
        out[:, slots] = c[np.argsort(keys, axis=1)[:, :len(slots)]]
    return out


def pval(null, x):
    return float((1 + (null >= x).sum()) / (1 + len(null)))


# ----------------------------------------------------------------------------- main analysis
def run(args):
    t_start = time.time()
    info: dict = {}
    CACHE.mkdir(exist_ok=True)
    plate_cache = Path(args.plate_cache) if args.plate_cache else None
    if plate_cache:
        plate_cache.mkdir(parents=True, exist_ok=True)

    acc = resolve("lincs")
    log(f"accession: {acc}; batch {BATCH}")
    bp, pm = fetch_metadata(acc)
    plates = list_plates(acc)
    _, flat = s3_list(f"{acc}/broad/workspace/profiles/{BATCH}/", delimiter="")
    keyset = {k for k, _ in flat}
    for p, k in plates.items():
        if k not in keyset:
            raise RuntimeError(f"augmented profile missing for plate {p}: {k}")
    info["plate_file_sizes_bytes"] = int(sum(s for k, s in flat if k.endswith("_augmented.csv.gz")))
    bp_set = set(bp.Assay_Plate_Barcode)
    missing = sorted(bp_set - set(plates))
    extra = sorted(set(plates) - bp_set)
    assert not extra, f"profiled plates not in barcode_platemap: {extra}"
    bpm = bp.set_index("Assay_Plate_Barcode").Plate_Map_Name
    info["n_plates_barcode_map"] = len(bp_set)
    info["n_plates_profiled"] = len(plates)
    info["plates_missing"] = {p: bpm[p] for p in missing}
    per_map = bpm[list(plates)].value_counts()
    assert per_map.min() >= 2, per_map
    info["plates_per_platemap"] = per_map.value_counts().sort_index().to_dict()
    log(f"plates: {len(plates)} profiled of {len(bp_set)} in barcode_platemap; missing {info['plates_missing']}")
    disk_check([plate_cache])

    frames, nbytes, t_dl = download(plates, plate_cache)
    info["bytes_profiles"] = int(sum(nbytes.values()))
    info["download_seconds"] = round(t_dl, 1)
    obs, X, feats = assemble(frames, pm)
    have = set(zip(obs.plate, obs.well))
    info["missing_wells"] = sorted(f"{p} {r}{c:02d}" for p in sorted(set(obs.plate)) for r in "ABCDEFGHIJKLMNOP"
                                   for c in range(1, 25) if (p, f"{r}{c:02d}") not in have)
    info["n_wells_total"] = len(obs)
    del frames
    info.update(ASSEMBLY)
    log(f"wells: {len(obs)}, raw feature columns {len(feats)}, DMSO wells {int(obs.is_dmso.sum())}")

    # dose mapping
    obs["dose_uM"] = nominal_dose(obs.dose_actual_uM.to_numpy())
    trt = ~obs.is_dmso
    off = trt & obs.dose_uM.isna()
    info["offladder_wells"] = int(off.sum())
    info["offladder_detail"] = (obs[off].groupby(["broad_sample", "pert_iname"], dropna=False)
                                .dose_actual_uM.apply(lambda s: ";".join(f"{v:.4g}" for v in sorted(set(s))))
                                .reset_index().to_dict("records"))
    act_doses = (obs[trt & obs.dose_uM.notna()].groupby("dose_uM").dose_actual_uM
                 .agg(lambda s: sorted({round(v, 4) for v in s})))
    info["actual_doses_by_nominal"] = {float(k): v for k, v in act_doses.items()}
    # duplicate wells of one compound at one nominal dose on one plate would mean bin collision
    coll = obs[trt & obs.dose_uM.notna()].groupby(["plate", "b13", "dose_uM"]).size()
    info["bin_collisions"] = int((coll > 1).sum())
    log(f"off-ladder treated wells: {info['offladder_wells']}; bin collisions {info['bin_collisions']}")

    # units
    units, ids = load_units()
    info["n_units_tier"] = int(units.unit.nunique())
    shared = ids.groupby("bs").unit.nunique()
    info["ids_shared_by_units"] = {b: sorted(ids[ids.bs == b].unit) for b in shared[shared > 1].index}
    miss_ids = sorted(set(ids.bs) - set(obs.broad_sample))
    info["unit_ids_not_in_batch"] = miss_ids

    # 1-4 processing
    Z, maxabs, f = normalize(obs, X, feats, info)
    del X
    Zs, fsel = select_features(Z, maxabs, f, info)
    del Z
    log(f"features: raw {info['n_features_raw']}, bookkeeping -{info['drop_nonmorph_bookkeeping']}, "
        f"nan>1% -{info['drop_nan_gt1pct']}, const -{info['drop_constant']} -> {info['n_after_nan_const']}; "
        f"blocklist -{info['drop_blocklist']}, outlier -{info['drop_outlier_gt500']}, "
        f"lowvar -{info['drop_low_variance']}, corr -{info['drop_correlated']} -> {info['n_selected']}")
    dm = obs.is_dmso.to_numpy()
    info["n_dmso"] = int(dm.sum())
    mu, W = fit_sphere(Zs[dm].astype(np.float64), info)
    V = ((Zs.astype(np.float64) - mu) @ W).astype(np.float32)
    del Zs
    U = unitv(V.astype(np.float64))
    if args.dump:
        Path(args.dump).mkdir(parents=True, exist_ok=True)
        obs.to_parquet(Path(args.dump) / "obs.parquet")
        np.save(Path(args.dump) / "V.npy", V)
    log(f"sphering: {info['sphere_dim']} dims, eig first/last {info['sphere_eig_first_last']}, eps {SPHERE_EPS}")

    rng = np.random.default_rng(SEED)
    plate_arr = obs.plate.to_numpy()
    well_arr = obs.well.to_numpy()
    plate_idx = {p: np.sort(ix) for p, ix in obs.groupby("plate").indices.items()}
    dmso_by_plate = {p: ix[dm[ix]] for p, ix in plate_idx.items()}
    onl = trt.to_numpy() & obs.dose_uM.notna().to_numpy()
    dose_arr = obs.dose_uM.to_numpy()
    b13_arr = obs.b13.to_numpy()

    # ---- conditions: matched units x dose, and all on-ladder compounds x dose (norm_z reference)
    id2units = ids.groupby("bs").unit.apply(list).to_dict()
    unit_ids = ids.groupby("unit").bs.apply(set).to_dict()
    cond_rows = []
    wsel = onl & obs.broad_sample.isin(id2units).to_numpy()
    for i in np.flatnonzero(wsel):
        for un in id2units[obs.broad_sample.iat[i]]:
            cond_rows.append((un, dose_arr[i], i))
    cw = pd.DataFrame(cond_rows, columns=["unit", "dose_uM", "row"])
    conds = cw.groupby(["unit", "dose_uM"]).row.apply(lambda s: np.sort(s.to_numpy())).reset_index()
    log(f"matched conditions: {len(conds)} ({conds.unit.nunique()} units)")

    # ---- shared structures for the tests
    bpm_batch = bp.set_index("Assay_Plate_Barcode").Batch_Number
    batch_arr = obs.plate.map(bpm_batch).to_numpy()
    assert not pd.isna(batch_arr).any()
    is_b1 = batch_arr == 1
    map_arr = obs.plate_map.to_numpy()
    plates_sorted = sorted(plate_idx)
    pid = {p: j for j, p in enumerate(plates_sorted)}
    plate_id_arr = obs.plate.map(pid).to_numpy()
    DM = np.flatnonzero(dm)
    dloc = np.full(len(obs), -1, np.int64)
    dloc[DM] = np.arange(len(DM))
    for p in plates_sorted:
        assert len(dmso_by_plate[p]) == 24, f"plate {p} has {len(dmso_by_plate[p])} DMSO wells"
    plate_dmso = np.stack([dloc[dmso_by_plate[p]] for p in plates_sorted])  # (n_plates, 24) local
    U32 = U.astype(np.float32)
    UD = U32[DM]
    V64 = V

    def is_edge(w: str) -> bool:
        return w[0] in "AP" or int(w[1:]) in (1, 24)

    def score_sets(Q, exclude_self):
        """Q: (D, n) global rows, one well per distinct plate. Each well queries the other
        n-1 plus all DMSO wells on the set's plates (minus the set itself if exclude_self).
        Returns (mAP (D,), norm of median profile (D,))."""
        D, n = Q.shape
        UQ = U32[Q]
        pos = np.einsum("dik,djk->dij", UQ, UQ)
        negs = plate_dmso[plate_id_arr[Q]].reshape(D, n * 24)
        uq, inv = np.unique(Q, return_inverse=True)
        S = U32[uq] @ UD.T
        ng = S[inv.reshape(D, n)[:, :, None], negs[:, None, :]]
        valid = None
        if exclude_self:
            valid = ~(negs[:, None, :] == dloc[Q][:, :, None]).any(1)
        off_ = ~np.eye(n, dtype=bool)
        pp = -np.sort(-pos[:, off_].reshape(D, n, n - 1), axis=-1)
        if valid is not None:
            ng = np.where(valid[:, None, :], ng, -np.inf)
        above = (ng[:, :, None, :] > pp[:, :, :, None]).sum(-1)
        k = np.arange(1, n, dtype=float)
        mAP = (k / (k + above)).mean(-1).mean(-1)
        norms = np.linalg.norm(np.median(V64[Q].astype(np.float64), axis=1), axis=-1)
        return mAP, norms

    # (1) plate-matched random DMSO null (spec step 6): cached per plate set
    null_cache = {}

    def dmso_null(rows):
        key = tuple(sorted(plate_arr[rows]))
        if key not in null_cache:
            slot_c = [dmso_by_plate[p] for p in key]
            Q = draw_slots(rng, slot_c, N_NULL)
            null_cache[key] = score_sets(Q, exclude_self=True)
        return null_cache[key]

    # (2) position-matched DMSO null: whole other-map same-position DMSO sets
    dmso_wells = sorted(set(well_arr[dm]))
    maps_all = sorted(set(map_arr))
    mapset = {}
    for j in DM:
        mapset.setdefault((map_arr[j], well_arr[j]), []).append(j)
    mapset = {k: np.array(sorted(v, key=lambda r: plate_arr[r])) for k, v in mapset.items()}
    set_score_cache = {}

    def comp_of(rows):
        return int(is_b1[rows].sum()), int((~is_b1[rows]).sum())

    def trimmed(M, w, c):
        rows = mapset[(M, w)]
        b1, bk = rows[is_b1[rows]], rows[~is_b1[rows]]
        if len(b1) < c[0] or len(bk) < c[1]:
            return None
        return np.sort(np.concatenate([b1[:c[0]], bk[:c[1]]]))

    def set_score(M, w, c):
        key = (M, w, c)
        if key not in set_score_cache:
            r = trimmed(M, w, c)
            set_score_cache[key] = None if r is None else tuple(float(x[0]) for x in score_sets(r[None], True))
        return set_score_cache[key]

    def groups_of(rows):
        out = []
        for M in sorted(set(map_arr[rows])):
            g = rows[map_arr[rows] == M]
            ws = set(well_arr[g])
            assert len(ws) == 1, f"replicates of one map at several positions: {ws}"
            out.append((M, g, comp_of(g), is_edge(next(iter(ws)))))
        return out

    posnull_cache = {}

    def pos_dmso_null(rows):
        grp = groups_of(rows)
        own = frozenset(g[0] for g in grp)
        sig = (own, tuple((g[2], g[3]) for g in grp))
        if sig in posnull_cache:
            return posnull_cache[sig]
        cand = []
        for _, _, c, e in grp:
            cand.append([(M, w) for M in maps_all if M not in own for w in dmso_wells
                         if is_edge(w) == e and trimmed(M, w, c) is not None])
        if len(grp) == 1:
            sc = [set_score(M, w, grp[0][2]) for M, w in cand[0]]
            res = (np.array([x[0] for x in sc]), np.array([x[1] for x in sc]))
        else:
            Qs = []
            for _ in range(N_NULL):
                while True:
                    pick = [cand[i][rng.integers(len(cand[i]))] for i in range(len(grp))]
                    if len({m for m, _ in pick}) == len(pick):
                        break
                Qs.append(np.concatenate([trimmed(M, w, grp[i][2]) for i, (M, w) in enumerate(pick)]))
            res = score_sets(np.array(Qs), exclude_self=True)
        posnull_cache[sig] = res
        return res

    # (3) other compounds, same plates, same dose (addendum b)
    def cpd_null(rows, dose, own_b13):
        ps = sorted(set(plate_arr[rows]))
        cand = np.concatenate([plate_idx[p] for p in ps])
        cand = cand[onl[cand] & (dose_arr[cand] == dose) & ~np.isin(b13_arr[cand], list(own_b13))]
        cpds = np.unique(b13_arr[cand])
        cidx = {c: j for j, c in enumerate(cpds)}
        tab = {p: np.full(len(cpds), -1, np.int64) for p in ps}
        for g in cand:
            tab[plate_arr[g]][cidx[b13_arr[g]]] = g
        keys = rng.random((N_NULL, len(cpds)))
        used = np.zeros((N_NULL, len(cpds)), bool)
        Q = np.empty((N_NULL, len(rows)), np.int64)
        ar = np.arange(N_NULL)
        for s, p in enumerate(plate_arr[rows]):
            kk = np.where((tab[p] >= 0)[None] & ~used, keys, np.inf)
            c = kk.argmin(1)
            assert np.isfinite(kk[ar, c]).all(), "not enough null compounds"
            used[ar, c] = True
            Q[:, s] = tab[p][c]
        return score_sets(Q, exclude_self=False)[0]

    # (4) position-matched other compounds: same well, same batch, distinct other maps per slot
    rows_at = {}
    for j in np.flatnonzero(onl):
        rows_at.setdefault((well_arr[j], batch_arr[j], map_arr[j]), []).append(j)

    def cpd_pos_null(rows, own_b13):
        own = set(map_arr[rows])
        n = len(rows)
        tab = np.full((n, len(maps_all), 4), -1, np.int64)
        cnt = np.zeros((n, len(maps_all)), np.int64)
        for s, r in enumerate(rows):
            for mi, M in enumerate(maps_all):
                if M in own:
                    continue
                c = [g for g in rows_at.get((well_arr[r], batch_arr[r], M), []) if b13_arr[g] not in own_b13]
                cnt[s, mi] = len(c)
                tab[s, mi, :len(c)] = c
        keys = rng.random((N_NULL, len(maps_all)))
        used = np.zeros((N_NULL, len(maps_all)), bool)
        Q = np.empty((N_NULL, n), np.int64)
        ar = np.arange(N_NULL)
        for s in range(n):
            kk = np.where((cnt[s] > 0)[None] & ~used, keys, np.inf)
            c = kk.argmin(1)
            if not np.isfinite(kk[ar, c]).all():
                return None
            used[ar, c] = True
            pick = (rng.random(N_NULL) * cnt[s, c]).astype(np.int64)
            Q[:, s] = tab[s, c, pick]
        return score_sets(Q, exclude_self=False)[0]

    def evaluate(rows, dose, own_b13, do_cpd=True):
        rows = np.sort(rows)
        assert len(set(plate_arr[rows])) == len(rows), "two wells of one condition on one plate"
        rec = {"n_wells": len(rows)}
        cons = np.median(V64[rows].astype(np.float64), 0)
        rec["norm"] = float(np.linalg.norm(cons))
        rec["_cons"] = cons.astype(np.float32)
        nanf = dict(rep_mAP=np.nan, rep_p=np.nan, rep_null_mean=np.nan, norm_dmso_p=np.nan,
                    norm_dmso_null_mean=np.nan, rep_pos_p=np.nan, norm_pos_p=np.nan, n_pos_null=0,
                    vs_cpd_p=np.nan, cpd_null_mean=np.nan, vs_cpd_pos_p=np.nan, cpd_pos_null_mean=np.nan)
        rec.update(nanf)
        if len(rows) < 2:
            return rec
        x = float(score_sets(rows[None], exclude_self=False)[0][0])
        mn, nn = dmso_null(rows)
        pm_, pn_ = pos_dmso_null(rows)
        rec.update(rep_mAP=x, rep_p=pval(mn, x), rep_null_mean=float(mn.mean()),
                   norm_dmso_p=pval(nn, rec["norm"]), norm_dmso_null_mean=float(nn.mean()),
                   rep_pos_p=pval(pm_, x), norm_pos_p=pval(pn_, rec["norm"]), n_pos_null=len(pm_),
                   rep_pos_null_mean=float(pm_.mean()), norm_pos_null_mean=float(pn_.mean()))
        if do_cpd:
            cn = cpd_null(rows, dose, own_b13)
            rec.update(vs_cpd_p=pval(cn, x), cpd_null_mean=float(cn.mean()))
            cp = cpd_pos_null(rows, own_b13)
            if cp is not None:
                rec.update(vs_cpd_pos_p=pval(cp, x), cpd_pos_null_mean=float(cp.mean()))
            else:
                rec.update(vs_cpd_pos_p=np.nan, cpd_pos_null_mean=np.nan)
        return rec

    # ---- matched conditions
    recs = []
    t0 = time.time()
    for r in conds.itertuples():
        rows = r.row
        ps = sorted(set(plate_arr[rows]))
        rec = {"unit": r.unit, "dose_uM": r.dose_uM, "n_plates": len(ps), "plates": ";".join(ps),
               "sites": SITE, "well_position": ";".join(sorted(set(well_arr[rows]))),
               "plate_maps": ";".join(sorted(set(map_arr[rows]))),
               "dose_actual_uM": ";".join(f"{v:.4g}" for v in sorted(set(obs.dose_actual_uM.to_numpy()[rows])))}
        rec.update(evaluate(rows, r.dose_uM, {s_[:13] for s_ in unit_ids[r.unit]}))
        wells = well_arr[rows]
        a_, b_ = np.triu_indices(len(rows), 1)
        rec["frac_rep_pairs_same_well"] = float((wells[a_] == wells[b_]).mean()) if len(a_) else np.nan
        recs.append(rec)
    log(f"matched-condition tests: {time.time() - t0:.0f}s")
    act = pd.DataFrame(recs)

    # ---- library: every on-ladder treated compound x dose (norm_z reference + calibration)
    t0 = time.time()
    lib = (pd.DataFrame({"b13": b13_arr[onl], "dose_uM": dose_arr[onl], "row": np.flatnonzero(onl)})
           .groupby(["b13", "dose_uM"]).row.apply(np.array).reset_index())
    lrecs = []
    for r in lib.itertuples():
        e = evaluate(r.row, r.dose_uM, {r.b13}, do_cpd=args.library_cpd_nulls)
        e.pop("_cons")
        e.update(b13=r.b13, dose_uM=r.dose_uM, well_position=";".join(sorted(set(well_arr[r.row]))),
                 plate_maps=";".join(sorted(set(map_arr[r.row]))))
        lrecs.append(e)
    lib = pd.DataFrame(lrecs)
    log(f"library tests ({len(lib)} compound x dose): {time.time() - t0:.0f}s")
    info["norm_ref_n_per_dose"] = lib.groupby("dose_uM").size().to_dict()
    zstats = lib.groupby("dose_uM").norm.agg(med="median", mad=lambda s_: np.median(np.abs(s_ - s_.median())))

    def add_calls(df):
        df = df.join(zstats, on="dose_uM")
        df["norm_z"] = (df.norm - df.med) / (MAD_SCALE * df.mad)
        df = df.drop(columns=["med", "mad"])
        def call(p, thr, greater=False):
            c = (p > thr) if greater else (p < thr)
            return c.astype("boolean").mask(p.isna())
        df["active_rep"] = call(df.rep_p, P_ACTIVE)
        df["active_z"] = call(df.norm_z, Z_ACTIVE, greater=True)
        df["active_norm_dmso"] = call(df.norm_dmso_p, P_ACTIVE)
        df["active_vs_cpd"] = call(df.vs_cpd_p, P_ACTIVE)
        df["active_rep_pos"] = call(df.rep_pos_p, P_ACTIVE)
        df["active_norm_pos"] = call(df.norm_pos_p, P_ACTIVE)
        df["active_vs_cpd_pos"] = call(df.vs_cpd_pos_p, P_ACTIVE)
        return df

    act = add_calls(act)
    lib = add_calls(lib)
    CALLS = ["active_rep", "active_z", "active_norm_dmso", "active_vs_cpd", "active_rep_pos",
             "active_norm_pos", "active_vs_cpd_pos"]
    lib[CALLS] = lib[CALLS].astype(float)
    lib_rate = lib.groupby("dose_uM")[CALLS].mean()
    info["library_active_frac_by_dose"] = lib_rate.round(4).reset_index().to_dict("records")
    log("library active fraction by dose:\n" + lib_rate.round(3).to_string())

    # ---- position quality. Layout (identical on all 28 maps): each compound fills one row x
    # 6-column block (doses 10 -> 0.04 uM left to right). DMSO sits only in A/B01-06, F07-12, J13-18,
    # so DMSO nulls cannot sample columns 19-24 or the outer rows. Reference region chosen from the
    # per-well library heatmap: rows C-N, columns 1-18 (interior, near DMSO positions).
    from scipy.stats import binom
    single = lib[~lib.well_position.str.contains(";")].copy()
    single["row_"] = single.well_position.str[0]
    single["col_"] = single.well_position.str[1:].astype(int)
    single["block"] = single.row_ + "_" + ((single.col_ - 1) // 6 + 1).astype(str)
    ref_m = single.row_.isin(list("CDEFGHIJKLMN")) & (single.col_ <= 18)
    POSCALLS = ["active_rep_pos", "active_norm_pos"]
    ref_rate = single[ref_m].groupby("dose_uM")[POSCALLS].mean()
    info["ref_region_rate_by_dose"] = ref_rate.round(4).reset_index().to_dict("records")
    flags = {}
    for c_ in POSCALLS:
        g = single.groupby(["well_position", "dose_uM"])[c_].agg(["sum", "count"]).reset_index()
        g["rate"] = g.dose_uM.map(ref_rate[c_])
        g["p"] = binom.sf(g["sum"] - 1, g["count"], g["rate"])
        fl = g[g.p < 0.05 / len(g)]
        flags[c_] = set(zip(fl.well_position, fl.dose_uM))
        info[f"n_flagged_wells_{c_}"] = len(fl)
        info[f"position_rate_spread_{c_}"] = (g.assign(frac=g["sum"] / g["count"]).groupby("dose_uM").frac
                                              .agg(["min", "median", "max"]).round(3).reset_index()
                                              .to_dict("records"))
    log(f"flagged wells (vs reference region, Bonferroni): rep_pos {len(flags['active_rep_pos'])}, "
        f"norm_pos {len(flags['active_norm_pos'])}")
    # block reliability: library rate at the block's lowest-dose (0.0412 uM) well(s), all maps
    lowd = single[single.dose_uM == NOMINAL_UM[0]]
    blk = lowd.groupby("block")[POSCALLS].agg(["sum", "count"])
    blk.columns = [f"{a_}_{b_}" for a_, b_ in blk.columns]
    ref0 = ref_rate.loc[NOMINAL_UM[0]]
    pb = {}
    for c_ in POSCALLS:
        blk[f"{c_}_rate"] = blk[f"{c_}_sum"] / blk[f"{c_}_count"]
        pb[c_] = binom.sf(blk[f"{c_}_sum"] - 1, blk[f"{c_}_count"], ref0[c_])
    blk["block_pass_binom"] = (pb["active_rep_pos"] >= 0.01) & (pb["active_norm_pos"] >= 0.01)
    # columns 19-24 hold no DMSO well, so no DMSO null samples them: excluded regardless of the test
    # (rule read from the per-well heatmap; the binomial test alone lacks power with 28 compounds)
    blk["block_in_cols_1_18"] = blk.index.str.split("_").str[1].astype(int) <= 3
    blk["block_reliable"] = blk.block_pass_binom & blk.block_in_cols_1_18
    info["n_blocks_pass_binom"] = int(blk.block_pass_binom.sum())
    info["n_blocks"] = len(blk)
    info["n_blocks_reliable"] = int(blk.block_reliable.sum())
    info["unreliable_blocks"] = sorted(blk.index[~blk.block_reliable])
    log(f"blocks: {len(blk)}, reliable {int(blk.block_reliable.sum())}; unreliable {info['unreliable_blocks']}")
    well_block = lambda w: f"{w[0]}_{(int(w[1:]) - 1) // 6 + 1}"  # noqa: E731

    def block_info(well_positions):
        bl = sorted({well_block(w) for w in well_positions})
        ok = all(blk.block_reliable.get(b_, False) for b_ in bl)
        r1 = float(np.mean([blk[f"active_rep_pos_rate"].get(b_, np.nan) for b_ in bl]))
        r2 = float(np.mean([blk[f"active_norm_pos_rate"].get(b_, np.nan) for b_ in bl]))
        return ";".join(bl), ok, r1, r2

    for df in (act, lib):
        bi = [block_info(wp.split(";")) for wp in df.well_position]
        df["block"] = [x[0] for x in bi]
        df["block_reliable"] = [x[1] for x in bi]
        df["block_lowdose_rate_rep_pos"] = [x[2] for x in bi]
        df["block_lowdose_rate_norm_pos"] = [x[3] for x in bi]
    lib_rel_rate = lib[lib.block_reliable].groupby("dose_uM")[CALLS].mean()
    info["library_active_frac_by_dose_reliable_blocks"] = lib_rel_rate.round(4).reset_index().to_dict("records")
    log("library active fraction by dose, reliable blocks only:\n" + lib_rel_rate.round(3).to_string())
    # region rates at the lowest dose (all three regions)
    reg = np.where(ref_m, "rows C-N cols 1-18", np.where(single.col_ >= 19, "cols 19-24",
                   "outer rows (A,B,O,P) cols 1-18"))
    info["lowdose_rate_by_region"] = (single.assign(region=reg)[single.dose_uM == NOMINAL_UM[0]]
                                      .groupby("region")[POSCALLS].mean().round(4).to_dict("index"))
    info["rate_by_region_and_dose"] = (single.assign(region=reg).groupby(["region", "dose_uM"])[POSCALLS]
                                       .mean().round(3).reset_index().to_dict("records"))
    # library sustained-dose bound: fraction of reliable-block compounds whose sustained dose is 0.0412
    sb = {}
    for c_ in POSCALLS:
        n0 = ntot = 0
        for _, g in lib[lib.block_reliable].sort_values("dose_uM").groupby("b13"):
            g = g[g[c_].notna()]
            v = g[c_].astype(bool).to_numpy()
            if len(v) < 6:
                continue
            ntot += 1
            n0 += bool(v.all())
        sb[c_] = {"n_compounds": ntot, "frac_sustained_at_0.0412": round(n0 / max(ntot, 1), 4)}
    info["library_sustained_lowest_bound"] = sb
    log(f"library (reliable blocks) sustained-from-0.0412 fraction: {sb}")
    # persistence 0.04 -> 0.12 uM inside reliable blocks
    wp_ = lib[lib.block_reliable].pivot_table(index="b13", columns="dose_uM", values="active_rep_pos",
                                              aggfunc="first")
    a0, a1 = wp_.get(NOMINAL_UM[0]), wp_.get(NOMINAL_UM[1])
    if a0 is not None and a1 is not None:
        m_ = a0.notna() & a1.notna()
        info["persistence_rep_pos_reliable"] = {
            "P(act0.123|act0.0412)": round(float(a1[m_ & (a0 == 1)].astype(float).mean()), 3),
            "P(act0.123|not0.0412)": round(float(a1[m_ & (a0 == 0)].astype(float).mean()), 3)}

    # ---- positive controls: the two 20 uM compounds on every plate (one same-position set per map)
    pcs = []
    hi = off & (obs.dose_actual_uM > 15)
    for (bs, mp, w), gp in obs[hi].groupby(["broad_sample", "plate_map", "well"]):
        if len(gp) < 2:
            continue
        e = evaluate(gp.index.to_numpy(), np.nan, {bs[:13]}, do_cpd=False)
        pcs.append({"broad_sample": bs, "pert_iname": gp.pert_iname.iat[0], "plate_map": mp, "well": w,
                    **{k: e[k] for k in ("rep_mAP", "rep_p", "norm", "norm_dmso_p", "rep_pos_p", "norm_pos_p")}})
    pcs = pd.DataFrame(pcs)
    info["posctrl"] = (pcs.groupby(["broad_sample", "pert_iname"])
                       .agg(n_sets=("plate_map", "size"), median_rep_mAP=("rep_mAP", "median"),
                            frac_rep_p05=("rep_p", lambda s_: float((s_ < 0.05).mean())),
                            frac_rep_pos_p05=("rep_pos_p", lambda s_: float((s_ < 0.05).mean())),
                            median_norm=("norm", "median"),
                            frac_norm_dmso_p05=("norm_dmso_p", lambda s_: float((s_ < 0.05).mean())),
                            frac_norm_pos_p05=("norm_pos_p", lambda s_: float((s_ < 0.05).mean())))
                       .reset_index().round(4).to_dict("records"))
    log(f"positive controls: {info['posctrl']}")

    # ---- calibration: within-map same-position DMSO sets through every DMSO-based test
    cal = []
    for (M, w), rows in mapset.items():
        if len(rows) < 2:
            continue
        x, nm = (float(v[0]) for v in score_sets(rows[None], exclude_self=True))
        mn, nn = dmso_null(rows)
        pm_, pn_ = pos_dmso_null(rows)
        cal.append({"plate_map": M, "well": w, "n": len(rows), "edge": is_edge(w), "map": x,
                    "rep_p": pval(mn, x), "norm_dmso_p": pval(nn, nm),
                    "rep_pos_p": pval(pm_, x), "norm_pos_p": pval(pn_, nm), "n_pos_null": len(pm_)})
    cal = pd.DataFrame(cal)
    cinfo = {"n_sets": len(cal), "n_edge_sets": int(cal.edge.sum()),
             "pos_null_sizes": [int(cal.n_pos_null.min()), int(cal.n_pos_null.max())]}
    for col in ["rep_p", "norm_dmso_p", "rep_pos_p", "norm_pos_p"]:
        cinfo[f"{col}<0.05"] = round(float((cal[col] < 0.05).mean()), 4)
        cinfo[f"{col}<0.05_edge"] = round(float((cal[col][cal.edge] < 0.05).mean()), 4)
        cinfo[f"{col}<0.05_interior"] = round(float((cal[col][~cal.edge] < 0.05).mean()), 4)
    info["calib_same_pos_dmso"] = cinfo
    log(f"same-position DMSO calibration: {cinfo}")

    # ---- position diagnostics (addendum c)
    dI = np.flatnonzero(dm)
    Ud = U[dI]
    Gd = Ud @ Ud.T
    a, b = np.triu_indices(len(dI), 1)
    diffplate = plate_arr[dI][a] != plate_arr[dI][b]
    samew = well_arr[dI][a] == well_arr[dI][b]
    gd = Gd[a, b]
    info["dmso_cos_same_pos_diff_plate"] = float(gd[diffplate & samew].mean())
    info["dmso_cos_diff_pos_diff_plate"] = float(gd[diffplate & ~samew].mean())
    info["dmso_cos_diff_pos_same_plate"] = float(gd[~diffplate].mean())
    samemap = map_arr[dI][a] == map_arr[dI][b]
    samebatch = batch_arr[dI][a] == batch_arr[dI][b]
    split = {}
    for nm_, m_ in [("same map, same batch", samemap & samebatch), ("same map, other batch", samemap & ~samebatch),
                    ("other map, same batch", ~samemap & samebatch),
                    ("other map, other batch", ~samemap & ~samebatch)]:
        mm = diffplate & samew & m_
        split[nm_] = {"n_pairs": int(mm.sum()), "cos": round(float(gd[mm].mean()), 4)}
    info["dmso_cos_same_pos_split"] = split
    # replicate and other-compound pair cosines among matched wells (cross-plate)
    mrows = np.unique(cw.row.to_numpy())
    Um = U[mrows]
    Gm = Um @ Um.T
    a, b = np.triu_indices(len(mrows), 1)
    lab = b13_arr[mrows]
    dsm = dose_arr[mrows]
    dp = plate_arr[mrows][a] != plate_arr[mrows][b]
    same = (lab[a] == lab[b]) & (dsm[a] == dsm[b])
    sw = well_arr[mrows][a] == well_arr[mrows][b]
    g = Gm[a, b]
    info["matched_rep_pairs_cross_plate"] = int((same & dp).sum())
    info["matched_rep_pairs_same_well_frac"] = float(sw[same & dp].mean())
    info["matched_rep_cos_cross_plate"] = float(g[same & dp].mean())
    info["other_cpd_cos_cross_plate_same_dose"] = float(g[~same & dp & (dsm[a] == dsm[b]) & (lab[a] != lab[b])].mean())
    info["other_cpd_cos_cross_plate_same_well"] = float(g[~same & dp & sw & (lab[a] != lab[b])].mean()) \
        if (~same & dp & sw & (lab[a] != lab[b])).any() else None
    log(f"DMSO cos same-pos/diff-plate {info['dmso_cos_same_pos_diff_plate']:.4f}, diff-pos/diff-plate "
        f"{info['dmso_cos_diff_pos_diff_plate']:.4f}; replicate cos {info['matched_rep_cos_cross_plate']:.4f} "
        f"(same well {info['matched_rep_pairs_same_well_frac']:.3f}); other-cpd same-dose cos "
        f"{info['other_cpd_cos_cross_plate_same_dose']:.4f}")

    # ---- outputs
    act = act.sort_values(["unit", "dose_uM"]).reset_index(drop=True)
    k = V.shape[1]
    fcols = [f"f{i}" for i in range(k)]
    cons = pd.DataFrame(np.stack(act._cons.to_numpy()), columns=fcols)
    dm_cons = np.median(V[dm].astype(np.float64), 0).astype(np.float32)
    meta = act[["unit", "dose_uM", "n_wells", "n_plates", "plates", "sites", "dose_actual_uM"]].copy()
    meta.insert(2, "cell_line", "A549")
    meta.insert(3, "time_h", 48)
    dmrow = pd.DataFrame([{"unit": "DMSO", "dose_uM": 0.0, "cell_line": "A549", "time_h": 48,
                           "n_wells": int(dm.sum()), "n_plates": len(plate_idx),
                           "plates": ";".join(sorted(plate_idx)), "sites": SITE, "dose_actual_uM": "0"}])
    cons_df = pd.concat([pd.concat([meta, cons], axis=1),
                         pd.concat([dmrow, pd.DataFrame(dm_cons[None], columns=fcols)], axis=1)],
                        ignore_index=True)
    cons_df[fcols] = cons_df[fcols].astype(np.float32)
    cons_df.to_parquet(OUT_CONS, index=False)
    TRACKED.add(OUT_CONS)

    acols = ["unit", "dose_uM", "cell_line", "time_h", "n_wells", "rep_mAP", "rep_p", "active_rep",
             "norm_z", "active_z", "norm", "norm_dmso_p", "active_norm_dmso", "vs_cpd_p", "active_vs_cpd",
             "rep_pos_p", "active_rep_pos", "norm_pos_p", "active_norm_pos", "vs_cpd_pos_p",
             "active_vs_cpd_pos", "pos_flag_rep", "pos_flag_norm",
             "rep_null_mean", "rep_pos_null_mean", "norm_dmso_null_mean", "norm_pos_null_mean",
             "cpd_null_mean", "cpd_pos_null_mean", "n_pos_null", "n_plates", "plate_maps", "well_position",
             "frac_rep_pairs_same_well", "dose_actual_uM"]
    act["cell_line"], act["time_h"] = "A549", 48
    keyset_ = list(zip(act.well_position, act.dose_uM))
    act["pos_flag_rep"] = [k_ in flags["active_rep_pos"] for k_ in keyset_]
    act["pos_flag_norm"] = [k_ in flags["active_norm_pos"] for k_ in keyset_]
    acols += ["block", "block_reliable", "block_lowdose_rate_rep_pos", "block_lowdose_rate_norm_pos"]
    act[acols].to_csv(OUT_ACT, index=False, float_format="%.6g")
    TRACKED.add(OUT_ACT)

    # wells: matched-unit wells (one row per unit; a well shared by two units appears twice) + DMSO
    wl = cw[["unit", "row"]].copy()
    wl = pd.concat([wl, pd.DataFrame({"unit": "DMSO", "row": np.flatnonzero(dm)})], ignore_index=True)
    wmeta = obs.iloc[wl.row.to_numpy()][["plate", "well", "plate_map", "broad_sample", "dose_actual_uM",
                                          "dose_uM", "is_dmso"]].reset_index(drop=True)
    wmeta["dose_uM"] = wmeta.dose_uM.fillna(0.0)
    wmeta.insert(0, "unit", wl.unit.to_numpy())
    wmeta.insert(1, "cell_line", "A549")
    wmeta.insert(2, "time_h", 48)
    wmeta["site"] = SITE
    wdf = pd.concat([wmeta, pd.DataFrame(V[wl.row.to_numpy()], columns=fcols)], axis=1)
    wdf.to_parquet(OUT_WELLS, index=False)
    TRACKED.add(OUT_WELLS)

    SHORT = {"active_rep": "rep", "active_z": "z", "active_norm_dmso": "norm_dmso",
             "active_vs_cpd": "vs_cpd", "active_rep_pos": "rep_pos", "active_norm_pos": "norm_pos",
             "active_vs_cpd_pos": "vs_cpd_pos"}

    def lowest(col):
        return act[act[col].fillna(False).astype(bool)].groupby("unit").dose_uM.min()

    def sustained(col):
        """Lowest dose from which every higher tested dose is also active (NaN if top tested
        dose inactive); tested = call not NaN."""
        out = {}
        for u_, g in act.sort_values("dose_uM").groupby("unit"):
            g = g[g[col].notna()]
            v = g[col].astype(bool).to_numpy()
            if not len(v) or not v[-1]:
                out[u_] = np.nan
                continue
            j = len(v) - 1
            while j > 0 and v[j - 1]:
                j -= 1
            out[u_] = g.dose_uM.to_numpy()[j]
        return pd.Series(out)

    low = pd.DataFrame({"unit": sorted(act.unit.unique())}).set_index("unit")
    low["lowest_active_dose_rep"] = lowest("active_rep")
    low["lowest_active_dose_z"] = lowest("active_z")
    low["n_doses_tested"] = act[act.rep_p.notna()].groupby("unit").size()
    low["n_doses_tested"] = low.n_doses_tested.fillna(0).astype(int)
    for c_ in CALLS[2:]:
        low[f"lowest_active_dose_{SHORT[c_]}"] = lowest(c_)
    for c_ in CALLS:
        low[f"sustained_active_dose_{SHORT[c_]}"] = sustained(c_)
    for col, fl in [("rep_pos", "pos_flag_rep"), ("norm_pos", "pos_flag_norm")]:
        m_ = act.merge(low[f"lowest_active_dose_{col}"].rename("lad"), left_on="unit", right_index=True)
        m_ = m_[m_.dose_uM == m_.lad]
        low[f"lowest_{col}_on_flagged_well"] = (m_.set_index("unit")[fl].reindex(low.index)
                                                .astype("boolean").fillna(False))
    ub = act.groupby("unit").agg(block=("block", "first"), block_reliable=("block_reliable", "all"),
                                 block_lowdose_rate_rep_pos=("block_lowdose_rate_rep_pos", "first"),
                                 block_lowdose_rate_norm_pos=("block_lowdose_rate_norm_pos", "first"))
    low = low.join(ub)
    low.reset_index().to_csv(OUT_LOW, index=False)
    TRACKED.add(OUT_LOW)
    info["units_in_reliable_blocks"] = int(low.block_reliable.sum())

    # non-monotonic calls: active at a dose but inactive at some higher tested dose
    nonmono = {}
    for col in CALLS:
        n = 0
        for u_, g in act.sort_values("dose_uM").groupby("unit"):
            v = g[col].dropna().astype(bool).to_numpy()
            if v.any() and not v[np.argmax(v):].all():
                n += 1
        nonmono[col] = n
    info["n_units_nonmonotonic"] = nonmono

    for c_ in CALLS:
        act[c_ + "_f"] = act[c_].astype(float)
    info["n_matched_vs_cpd_pos_nan"] = int(act.vs_cpd_pos_p.isna().sum())
    per_dose = act.assign(**{c_: act[c_ + "_f"] for c_ in CALLS}).groupby("dose_uM").agg(
        n_units=("unit", "nunique"), median_n_wells=("n_wells", "median"), min_n_wells=("n_wells", "min"),
        n_wells=("n_wells", "sum"),
        frac_reproducible=("active_rep", "mean"), frac_active_z=("active_z", "mean"),
        frac_active_norm_dmso=("active_norm_dmso", "mean"), frac_active_vs_cpd=("active_vs_cpd", "mean"),
        frac_active_rep_pos=("active_rep_pos", "mean"), frac_active_norm_pos=("active_norm_pos", "mean"),
        frac_active_vs_cpd_pos=("active_vs_cpd_pos", "mean"), median_rep_mAP=("rep_mAP", "median"),
        n_units_reliable_block=("block_reliable", "sum"))
    rel = act[act.block_reliable].assign(**{c_: act[c_ + "_f"] for c_ in CALLS}).groupby("dose_uM")
    per_dose = per_dose.join(rel[["active_rep_pos", "active_norm_pos"]].mean()
                             .rename(columns=lambda c_: "frac_" + c_ + "_reliable_blocks"))
    wells_per_dose = (obs[onl].groupby("dose_uM").agg(wells_all_cpds=("well", "size"),
                                                      plates=("plate", "nunique")))
    per_dose = per_dose.join(wells_per_dose)
    uw = cw.drop_duplicates("row").groupby("dose_uM").agg(unique_wells=("row", "size"))
    uw["plates_matched"] = cw.drop_duplicates("row").groupby("dose_uM").row.apply(
        lambda r_: len(set(plate_arr[r_.to_numpy()])))
    per_dose = per_dose.join(uw)
    info["per_dose"] = per_dose.round(4).reset_index().to_dict("records")
    log(per_dose.round(3).to_string())

    lib_out = lib[["b13", "dose_uM", "n_wells", "well_position", "plate_maps", "rep_mAP", "rep_p", "norm",
                   "norm_z", "norm_dmso_p", "rep_pos_p", "norm_pos_p", "vs_cpd_p", "vs_cpd_pos_p"] + CALLS]
    lib_out.to_csv(OUT_LIB, index=False, float_format="%.6g")
    TRACKED.add(OUT_LIB)
    info["units_covered"] = int(act.unit.nunique())
    info["units_all6_doses"] = int((low.n_doses_tested == 6).sum())
    info["wells_out"] = len(wdf)
    info["wells_out_bytes"] = OUT_WELLS.stat().st_size
    info["seconds_total"] = round(time.time() - t_start, 1)
    disk_check([plate_cache])
    if plate_cache and not args.keep_plate_cache:
        for fp in plate_cache.glob("*.parquet"):
            fp.unlink()
    info["bytes_total"] = int(Bytes.n)
    info["peak_disk_bytes"] = int(DISK["peak"])
    OUT_DL.write_text(json.dumps({"bytes": int(Bytes.n), "peak_disk_bytes": int(DISK["peak"]),
                                  "seconds": info["seconds_total"]}, indent=1), encoding="utf-8")
    TRACKED.add(OUT_DL)
    (CACHE / "run_info.json").write_text(json.dumps(info, indent=1, default=str), encoding="utf-8")
    write_notes(acc, info, low, per_dose, plate_cache)
    log(f"done in {info['seconds_total']}s; bytes {Bytes.n / 1e6:.1f} MB; peak disk {DISK['peak'] / 1e6:.1f} MB")


def write_notes(acc, i, low, per_dose, plate_cache):
    pk = ["pandas", "numpy", "scipy", "scikit-learn", "pyarrow", "requests", "fsspec"]
    vers = {p: version(p) for p in pk}
    c = i["calib_same_pos_dmso"]
    L = []
    a = L.append
    a("# T1 LINCS Cell Painting (A549, 48 h, 6 doses): morphology profiles\n")
    a("Script: `step2/morph_lincs.py` (`py -3.13 -I -X utf8 step2/morph_lincs.py`). Metadata cache: "
      "`step2/lincs_cache/` (barcode_platemap.csv, 28 platemaps, run_info.json).\n")
    a("## Sources")
    a(f"- Accession `{acc}` (resolved via `gallery_index.resolve('lincs')`), batch `{BATCH}` "
      "(the Batch2 pilot `2017_12_05_Batch2` is not used).")
    a(f"- Profiles: `{acc}/broad/workspace/profiles/{BATCH}/<plate>/<plate>_augmented.csv.gz` "
      "(well-level, non-normalized aggregated CellProfiler features + platemap metadata; CSV, gzip). "
      "Streamed into memory, never written to disk.")
    a(f"- Platemaps: `{acc}/broad/workspace/metadata/platemaps/{BATCH}/` (barcode_platemap.csv + "
      "platemap/*.txt); cross-checked well by well against the augmented `Metadata_broad_sample` and "
      "`Metadata_mmoles_per_liter` (all agree).")
    a(f"- Plates: {i['n_plates_profiled']} profiled of {i['n_plates_barcode_map']} in barcode_platemap; "
      f"missing (no profiles in the gallery): {i['plates_missing']}. Plates per platemap: "
      f"{i['plates_per_platemap']} (count of plates -> number of platemaps).")
    a("- Compound key: `units.csv` tier `T1_LINCS_A549` `dataset_ids` are full 22-char broad_sample IDs "
      "(e.g. BRD-K70401845-003-09-6); all match `Metadata_broad_sample` exactly, so wells are mapped by full "
      "ID (each prefix has a single batch ID in this batch). The 13-char prefix keys compounds for the "
      "`norm_z` reference and the other-compound null. "
      f"IDs shared by two units (wells emitted under both): {i['ids_shared_by_units']}. "
      f"Unit IDs absent from the batch: {i['unit_ids_not_in_batch'] or 'none'}.")
    a("")
    a("## Doses")
    a("- `Metadata_mmoles_per_liter` is in fact micromolar (10 = 10 uM top dose). Each well is mapped to the "
      f"nearest nominal level in log space, {[float(x) for x in NOMINAL_UM]} uM (the `units.csv` labels), if within "
      f"{DOSE_TOL}x; otherwise the well is off-ladder (kept for normalization/feature selection/sphering, "
      "excluded from dose conditions and from the norm_z / other-compound references).")
    a(f"- Actual values per nominal dose: {i['actual_doses_by_nominal']}.")
    a(f"- Off-ladder treated wells: {i['offladder_wells']} (mostly the two 20 uM positive controls on every "
      f"plate); bin collisions (same compound, two actual doses on one nominal level on one plate): "
      f"{i['bin_collisions']}. All matched units sit exactly on the standard ladder.")
    a("")
    a("## Processing (one feature space for the whole batch)")
    a(f"- Feature columns: union {i['n_features_union']}; {len(i['plates_lacking_features'])} plates "
      f"({sorted(i['plates_lacking_features'])}) lack {len(i['features_not_on_all_plates'])} features "
      f"(all DNA RadialDistribution), so the plate intersection is used: {i['n_features_raw']} columns.")
    a(f"- Raw feature columns used: {i['n_features_raw']}. Dropped non-morphological bookkeeping columns "
      f"(`*_Location_*`, `*_Number_Object_Number`, `*_Parent_*`, `*ObjectNumber*`): "
      f"-{i['drop_nonmorph_bookkeeping']}; NaN/inf > 1% of wells: -{i['drop_nan_gt1pct']}; constant: "
      f"-{i['drop_constant']} -> {i['n_after_nan_const']}. Imputed values (plate DMSO median): "
      f"{i['n_values_imputed']}.")
    a(f"- Robust z per plate vs its 24 DMSO wells, 1.4826 MAD with floor {MAD_FLOOR_REL} x median-over-plates "
      f"scale; plate x feature scales at floor: {i['plate_feature_scales_floored_frac']:.4f}; clip +-{CLIP_Z:g}.")
    a(f"- Selection: pycytominer blocklist -{i['drop_blocklist']}, outlier (pre-clip max|z| > {OUTLIER_Z:g}) "
      f"-{i['drop_outlier_gt500']}, variance threshold (freq_cut {FREQ_CUT}, unique_cut {UNIQUE_CUT}) "
      f"-{i['drop_low_variance']}, greedy |r| > {CORR_THRESH} -{i['drop_correlated']} -> "
      f"**{i['n_selected']} features**.")
    a(f"- Sphering: PCA whitening on {i['n_dmso']} DMSO wells, {VAR_KEEP:.0%} variance -> "
      f"**{i['sphere_dim']} dims**, scale 1/sqrt(eigenvalue + {SPHERE_EPS}) (eps fixed a priori to the JUMP T3 "
      f"value, not tuned); eigenvalues first/last kept {i['sphere_eig_first_last']}, "
      f"{i['sphere_n_eig_below_eps']} kept eigenvalues below eps.")
    a("- Consensus: median of replicate wells per unit x nominal dose (single site `broad`; units with two "
      "IDs pool the wells of both).")
    a("")
    lr = pd.DataFrame(i["library_active_frac_by_dose"]).set_index("dose_uM")
    lrr = pd.DataFrame(i["library_active_frac_by_dose_reliable_blocks"]).set_index("dose_uM")
    rr = pd.DataFrame(i["ref_region_rate_by_dose"]).set_index("dose_uM")
    d0 = NOMINAL_UM[0]
    a("## Read this first: well position dominates the DMSO-based calls")
    a("- All replicate plates of a platemap share one layout and every dose of a compound sits at a fixed "
      f"well ({i['matched_rep_pairs_same_well_frac']:.1%} of {i['matched_rep_pairs_cross_plate']} cross-plate "
      "replicate pairs of matched units share a well; the rest are the two-ID units). Each compound fills one "
      "row x 6-column block (10 -> 0.04 uM left to right); DMSO occupies only A01-06, B01-06, F07-12, J13-18.")
    a(f"- DMSO cosine (sphered), different plates: same position {i['dmso_cos_same_pos_diff_plate']:.4f} vs "
      f"different position {i['dmso_cos_diff_pos_diff_plate']:.4f} (same plate, different position "
      f"{i['dmso_cos_diff_pos_same_plate']:.4f}). Same position, by map and batch: "
      + "; ".join(f"{k} {v['cos']:.3f} (n={v['n_pairs']})" for k, v in i["dmso_cos_same_pos_split"].items())
      + ". The position effect is batch- and map-specific. Matched replicate "
      f"pairs {i['matched_rep_cos_cross_plate']:.4f}; different compounds at the same dose "
      f"{i['other_cpd_cos_cross_plate_same_dose']:.4f}.")
    a(f"- Same-position DMSO sets ({c['n_sets']}, each DMSO position across the plates of one map) through the "
      f"spec tests: `active_rep` fires in {c['rep_p<0.05']:.1%} (edge {c['rep_p<0.05_edge']:.1%}, interior "
      f"{c['rep_p<0.05_interior']:.1%}), `active_norm_dmso` in {c['norm_dmso_p<0.05']:.1%} (edge "
      f"{c['norm_dmso_p<0.05_edge']:.1%}, interior {c['norm_dmso_p<0.05_interior']:.1%}). **The spec calls "
      "`active_rep` and `active_norm_dmso` (and `active_vs_cpd`) are therefore not usable as 'departs from "
      "control' in this dataset**; they are kept because the spec mandates them.")
    a("- Position-matched DMSO null (added): `rep_pos_p` / `norm_pos_p` compare the condition with whole "
      "same-position DMSO sets of other platemaps (same edge/interior class, trimmed to the condition's batch "
      f"composition: batch-1 plate + later-batch plates); null sizes {c['pos_null_sizes']} sets (fewer than "
      "1000; deviation). Leave-one-map-out on the DMSO sets gives "
      f"{c['rep_pos_p<0.05']:.1%} / {c['norm_pos_p<0.05']:.1%} (edge {c['rep_pos_p<0.05_edge']:.1%} / "
      f"{c['norm_pos_p<0.05_edge']:.1%}, interior {c['rep_pos_p<0.05_interior']:.1%} / "
      f"{c['norm_pos_p<0.05_interior']:.1%}); near 5% by construction, so not the false-positive "
      "rate at compound wells.")
    a("- Operational bound: the library (every on-ladder compound x dose, "
      f"{sum(i['norm_ref_n_per_dose'].values())} conditions) at the lowest dose. Rate of `active_rep_pos` / "
      f"`active_norm_pos` at 0.0412 uM: all wells {lr.loc[d0, 'active_rep_pos']:.1%} / "
      f"{lr.loc[d0, 'active_norm_pos']:.1%}; reference region (rows C-N, columns 1-18, chosen from "
      f"the per-well heatmap) {rr.loc[d0, 'active_rep_pos']:.1%} / {rr.loc[d0, 'active_norm_pos']:.1%}; "
      f"reliable blocks {lrr.loc[d0, 'active_rep_pos']:.1%} / {lrr.loc[d0, 'active_norm_pos']:.1%}. "
      f"Lowest-dose rate by region (rep_pos / norm_pos): "
      + "; ".join(f"{k} {v['active_rep_pos']:.1%} / {v['active_norm_pos']:.1%}"
                  for k, v in i["lowdose_rate_by_region"].items())
      + ". No DMSO well samples columns 19-24 or rows M-P, so no "
      "DMSO null can calibrate them, and a treated-well null there would be dose-relative.")
    a(f"- Block reliability: a block (row x 6 columns) is unreliable if its lowest-dose well's library rate "
      f"(over all 28 maps) of `active_rep_pos` or `active_norm_pos` exceeds the reference-region rate "
      f"(one-sided binomial p < 0.01; {i['n_blocks_pass_binom']} of {i['n_blocks']} pass) AND lies in columns "
      "1-18 (columns 19-24 have no DMSO well; rule read from the heatmap, where the binomial test lacks power "
      f"with 28 compounds per well). {i['n_blocks_reliable']} of {i['n_blocks']} blocks reliable; "
      f"unreliable: {i['unreliable_blocks']}. **Matched units in reliable blocks: "
      f"{i['units_in_reliable_blocks']} of {i['units_covered']}.** Wells flagged individually (Bonferroni, vs "
      f"reference region at the same dose): rep_pos {i['n_flagged_wells_active_rep_pos']}, norm_pos "
      f"{i['n_flagged_wells_active_norm_pos']}.")
    a(f"- Bound for the recommended column: fraction of library compounds in reliable blocks (all 6 doses "
      f"tested) whose sustained active dose is 0.0412 uM: {i['library_sustained_lowest_bound']}. This includes "
      "real low-dose activity, so it is an upper bound on the positional false-positive rate, not an estimate.")
    a(f"- Persistence inside reliable blocks: {i.get('persistence_rep_pos_reliable')} (adjacent wells of one "
      "block, so not a clean artifact test).")
    a("- `active_vs_cpd_pos` (other maps' compounds at the same well and batch, distinct maps per slot) lacks "
      "the same-map component of the position effect (0.106 vs 0.071) and is anticonservative; "
      f"{i['n_matched_vs_cpd_pos_nan']} matched conditions have no feasible draw (NaN).")
    a(f"- Positive controls (bortezomib, MG-132 at 20 uM; one same-position set per map and well): "
      f"{i['posctrl']}.")
    a("")
    a("**Recommended for the lowest-dose comparison:** `sustained_active_dose_rep_pos` (or `_norm_pos`) "
      "restricted to `block_reliable` units, i.e. the lowest dose from which every higher tested dose is "
      "active under the position-matched DMSO null; each dose is a different well, so isolated positional "
      "false positives drop out. Treat `lowest_active_dose_rep_pos` as an optimistic bound.")
    a("")
    a("## Calls (columns)")
    a(f"- `rep_mAP`: each replicate queries its other replicates + all DMSO wells on the condition's plates "
      f"(cosine). `rep_p`: {N_NULL} plate-matched random DMSO sets (one DMSO well per replicate plate; the "
      "drawn wells are removed from the negatives). `active_rep` = p < 0.05 (spec).")
    a(f"- `norm_z`: robust z of the consensus L2 norm vs all on-ladder compounds at that dose (n per dose "
      f"{i['norm_ref_n_per_dose']}); `active_z` = z > {Z_ACTIVE:g} (spec).")
    a(f"- `norm_dmso_p` / `active_norm_dmso` (addendum a): consensus norm vs {N_NULL} DMSO pseudo-consensus "
      "norms (median of one random DMSO well per replicate plate).")
    a(f"- `vs_cpd_p` / `active_vs_cpd` (addendum b): replicate mAP vs {N_NULL} sets of other compounds' wells on "
      "the same plates at the same dose (one distinct compound per slot).")
    a("- `rep_pos_p`, `norm_pos_p`, `vs_cpd_pos_p` and their `active_*_pos` calls: position-matched variants "
      "(above). `pos_flag_*`, `block`, `block_reliable`, `block_lowdose_rate_*`: position quality.")
    a("- Calls are NaN where the test could not be run (n_wells < 2 or no feasible null).")
    a("")
    a("### Per dose (matched units)")
    a("```")
    a(per_dose.round(3).T.to_string())
    a("```")
    a("")
    a("### Library active fraction by dose (all on-ladder compounds; compound nulls not run for the library)")
    a("```")
    a(lr.round(3).to_string())
    a("```")
    a("Reliable blocks only:")
    a("```")
    a(lrr.round(3).to_string())
    a("```")
    a(f"Units covered: {i['units_covered']} of {i['n_units_tier']}; with all 6 doses tested: "
      f"{i['units_all6_doses']}. Units with a non-monotonic call (active at a dose, inactive at a higher "
      f"tested dose): {i['n_units_nonmonotonic']}.")
    a("")
    a("## Outputs")
    a("- `step2/morph_lincs.parquet`: unit, dose_uM, cell_line, time_h, n_wells, n_plates, plates, sites, "
      "dose_actual_uM, f0..fN (float32); plus a `DMSO` reference row (dose 0, all DMSO wells).")
    a("- `step2/morph_lincs_library_activity.csv`: DMSO-based calls for every on-ladder compound x dose "
      "(13-char prefix), used for the norm_z reference and position diagnostics.")
    a("- `step2/morph_lincs_activity.csv`, `step2/morph_lincs_lowest_active.csv`, "
      f"`step2/morph_lincs_wells.parquet` ({i['wells_out']} rows, {i['wells_out_bytes'] / 1e6:.1f} MB; "
      "matched-unit wells + all DMSO wells; sphered).")
    a("")
    a("## Resources")
    a(f"- Bytes downloaded: {i['bytes_total']:,} (profiles {i['bytes_profiles']:,}); peak disk "
      f"{i['peak_disk_bytes']:,} bytes (metadata cache + outputs; profiles are held in memory only"
      f"{'; dev plate cache counted' if plate_cache else ''}); download {i['download_seconds']} s, total "
      f"{i['seconds_total']} s.")
    a(f"- Python {sys.version.split()[0]}; " + ", ".join(f"{k} {v}" for k, v in vers.items()))
    a("")
    a("## Deviations from MORPH_SPEC")
    a("- Non-morphological bookkeeping columns (Location, object numbers, parent IDs) dropped before step 1; "
      "they encode image coordinates, not morphology.")
    a("- Features absent on some plates (36 DNA RadialDistribution columns missing on 10 plates) are dropped "
      "for the whole batch (plate intersection) to keep one feature space.")
    a("- pycytominer blocklist applied in feature selection, as in the JUMP T3 run.")
    a("- Off-ladder wells (beyond 1.5x of a nominal level, chiefly the 20 uM positive controls) are excluded "
      "from dose conditions instead of being forced to the nearest level (20 uM would otherwise land on 10 uM).")
    a("- dose_uM uses the `units.csv` labels (0.0412, 0.123, 0.37, 1.11, 3.33, 10) rather than rounded values.")
    a("- In the DMSO null of `rep_mAP` the drawn DMSO wells are removed from the negatives, so null pools are "
      "smaller by n than the observed pools (slightly conservative).")
    a("- Added `active_norm_dmso` and `active_vs_cpd` per the coordinator addendum, plus position-matched "
      "variants, block reliability and sustained lowest doses (position effect, see top). The position-"
      f"matched DMSO null has {c['pos_null_sizes'][0]}-{c['pos_null_sizes'][1]} sets, fewer than the 1000 "
      "draws the spec asks for (only 24 DMSO "
      "positions x 27 other maps exist).")
    a(f"- Wells absent from the profile files: {i['missing_wells']} ({i['n_wells_total']:,} of "
      f"{i['n_plates_profiled'] * 384:,} wells).")
    a("")
    a("## Log")
    a("```")
    L.extend(LOG)
    a("```")
    OUT_NOTES.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--plate-cache", default=None, help="dev only: cache parsed plates here")
    ap.add_argument("--keep-plate-cache", action="store_true")
    ap.add_argument("--library-cpd-nulls", action="store_true",
                    help="also run the (slow, ~40 min) compound nulls for every library compound")
    ap.add_argument("--dump", default=None, help="dev only: dump obs + sphered wells here")
    run(ap.parse_args())
