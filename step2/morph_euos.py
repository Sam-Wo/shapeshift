"""T2 morphology profiles: EU-OPENSCREEN bioactives Cell Painting, HepG2, 24 h, 10 uM (controls 5 uM).

Follows step2/MORPH_SPEC.md (+ coordinator addendum on activity calls).

Source of the well-level CellProfiler profiles
----------------------------------------------
The Cell Painting Gallery accession (resolved via gallery_index.resolve('euos')) holds images,
load_data and platemaps only (no workspace/profiles or backend prefix). The authors' per-plate,
well-level aggregated CellProfiler profiles are on Zenodo (concept record 13309565; the latest
version is resolved through the REST API at run time). They sit inside one large zip
(`Aggregated_Profiles.zip`, ~1.35 GB) as one CSV per plate. We never store the zip: its central
directory is read with HTTP range requests and only the HepG2 plate members are fetched (again by
range request), inflated in memory, converted to float32 and cached per plate under
step2/euos_cache/plates/ (deleted at the end unless --keep-cache).

Run:  py -3.13 -I -X utf8 step2/morph_euos.py [--keep-cache]
"""
from __future__ import annotations

import io
import json
import re
import shutil
import struct
import sys
import time
import zlib
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import requests

STEP2 = Path(__file__).resolve().parent
ROOT = STEP2.parent
sys.path.insert(0, str(STEP2))  # -I drops the script dir from sys.path
from gallery_index import Bytes, resolve, s3_list, s3_url  # noqa: E402

CACHE = STEP2 / "euos_cache"
PLATE_CACHE = CACHE / "plates"
META_CACHE = CACHE / "metadata"
OUT_PROFILES = STEP2 / "morph_euos.parquet"
OUT_ACTIVITY = STEP2 / "morph_euos_activity.csv"
OUT_CONTROLS = STEP2 / "morph_euos_controls.csv"
OUT_WELLS = STEP2 / "morph_euos_wells.parquet"
OUT_NOTES = STEP2 / "euos_notes.md"
OUT_SUMMARY = STEP2 / "euos_download_summary.json"
OUT_RUN = CACHE / "euos_run.json"

ZENODO_CONCEPT = "13309565"  # concept DOI 10.5281/zenodo.13309565 (EU-OS bioactives, iScience 2025)
ZIP_NAME = "Aggregated_Profiles.zip"
TIER = "T2_EUOS_HepG2"
MATCH_DATASET = "EU-OS_cpg0036_HepG2+U2OS"
SITES = ["FMP", "IMTM", "MEDINA", "USC"]
DMSO = "DMSO"
CONTROLS = ["Nocodazole", "Tetrandrine"]
DOSE_UM = 10.0
CELL_LINE, TIME_H = "HepG2", 24
N_THREADS = 4

NAN_FRAC_MAX = 0.01
MAD_SCALE = 1.4826
MAD_FLOOR_REL = 0.01
MAD_FLOOR_ABS = 1e-12
OUTLIER_Z = 500.0
CLIP_Z = 50.0
FREQ_CUT, UNIQUE_CUT = 0.05, 0.01
CORR_THRESH = 0.9
CORR_SUBSAMPLE = 20000   # wells used to estimate the feature correlation matrix
VAR_KEEP = 0.99
SPHERE_EPS = 3.0         # same regularization as T3 (jump_profiles.py); DMSO features ~unit variance
N_NULL = 1000
N_NULL_CTRL = 1000
P_ACTIVE = 0.05
Z_ACTIVE = 3.0
SEED = 0

CHANNELS = ["AGP", "DNA", "ER", "Mito", "RNA"]
# pycytominer blocklist, translated to this dataset's compartment prefix ("Nuc_" not "Nuclei_")
BLOCKLIST = (
    [f"Nuc_Correlation_Manders_{a}_{b}" for a in CHANNELS for b in CHANNELS if a != b]
    + [f"Nuc_Correlation_RWC_{a}_{b}" for a in CHANNELS for b in CHANNELS if a != b]
    + [f"Nuc_Granularity_{g}_{c}" for g in (14, 15, 16) for c in CHANNELS]
)
# object bookkeeping / image-position columns (not morphology)
BOOKKEEPING = re.compile(r"(BoundingBoxM(ax|in)imum_[XY]|_Center_[XY]$|_Location_|ClosestObjectNumber"
                         r"|_ObjectNumber|_Parent_|_Children_|_Number_Object_Number)")

LOG: list[str] = []
STATE = {"peak_disk": 0, "zenodo_bytes": 0}
KEEP_CACHE = "--keep-cache" in sys.argv


def log(msg: str = "") -> None:
    print(msg, flush=True)
    LOG.append(msg)


def dir_bytes(p: Path) -> int:
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file()) if p.exists() else 0


def track_disk() -> None:
    STATE["peak_disk"] = max(STATE["peak_disk"], dir_bytes(CACHE))


# ----------------------------------------------------------------------------- HTTP helpers
def http_get(url: str, headers: dict | None = None, expect: int = 200, tries: int = 5, **kw) -> bytes:
    for k in range(tries):
        try:
            r = requests.get(url, headers=headers or {}, timeout=300, **kw)
            if r.status_code == expect:
                Bytes.add(len(r.content))
                return r.content
            if r.status_code in (404, 410, 416):
                break
            err = RuntimeError(f"HTTP {r.status_code} for {url} (expected {expect})")
        except requests.RequestException as e:  # retried, then raised
            err = e
        time.sleep(2 ** k)
    raise RuntimeError(f"source unreachable: {url} ({headers}) -> {locals().get('err', r.status_code)}")


def zenodo_latest() -> dict:
    """Latest version of the concept record (cached JSON) + its file list."""
    f = META_CACHE / "zenodo_record.json"
    if f.exists():
        return json.loads(f.read_text(encoding="utf-8"))
    raw = http_get(f"https://zenodo.org/api/records/{ZENODO_CONCEPT}")
    rec = json.loads(raw)
    vers = json.loads(http_get("https://zenodo.org/api/records",
                               params={"q": f"conceptrecid:{ZENODO_CONCEPT}", "all_versions": "true"}))
    out = {"concept": ZENODO_CONCEPT, "record_id": rec["id"], "version": rec["metadata"].get("version"),
           "title": rec["metadata"]["title"], "doi": rec.get("doi"),
           "publication_date": rec["metadata"].get("publication_date"),
           "files": [{"key": x["key"], "size": x["size"], "url": x["links"]["self"]} for x in rec["files"]],
           "all_versions": [{"id": h["id"], "version": h["metadata"].get("version"),
                             "date": h["metadata"].get("publication_date"),
                             "files": {x["key"]: x["size"] for x in h.get("files", [])}}
                            for h in vers["hits"]["hits"]]}
    META_CACHE.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


def zip_directory(url: str, size: int) -> pd.DataFrame:
    """Central directory of a remote zip via range requests (ZIP64 aware). Cached as CSV."""
    f = META_CACHE / "aggregated_profiles_zip_directory.csv"
    if f.exists():
        return pd.read_csv(f)
    a = size - 65536 - 22
    tail = http_get(url, {"Range": f"bytes={a}-{size - 1}"}, expect=206)
    i = tail.rfind(b"PK\x05\x06")
    if i < 0:
        raise RuntimeError("zip end-of-central-directory not found")
    _, _, _, _, n, cdsize, cdoff, _ = struct.unpack("<IHHHHIIH", tail[i:i + 22])
    j = tail.rfind(b"PK\x06\x06")
    if j >= 0:
        n, cdsize, cdoff = struct.unpack("<QQQ", tail[j + 24:j + 48])
    cd = http_get(url, {"Range": f"bytes={cdoff}-{cdoff + cdsize - 1}"}, expect=206)
    rows, p = [], 0
    while p < len(cd) and cd[p:p + 4] == b"PK\x01\x02":
        method, = struct.unpack("<H", cd[p + 10:p + 12])
        csz, usz, fnl, exl, cml = struct.unpack("<IIHHH", cd[p + 20:p + 34])
        lho, = struct.unpack("<I", cd[p + 42:p + 46])
        fn = cd[p + 46:p + 46 + fnl].decode("utf8")
        ex = cd[p + 46 + fnl:p + 46 + fnl + exl]
        q = 0
        while q + 4 <= len(ex):
            hid, hl = struct.unpack("<HH", ex[q:q + 4])
            if hid == 1:
                vals = list(struct.unpack("<" + "Q" * (hl // 8), ex[q + 4:q + 4 + 8 * (hl // 8)]))
                if usz == 0xFFFFFFFF:
                    usz = vals.pop(0)
                if csz == 0xFFFFFFFF:
                    csz = vals.pop(0)
                if lho == 0xFFFFFFFF:
                    lho = vals.pop(0)
            q += 4 + hl
        rows.append(dict(name=fn, method=method, csize=csz, usize=usz, offset=lho))
        p += 46 + fnl + exl + cml
    if len(rows) != n:
        raise RuntimeError(f"zip directory parse: {len(rows)} of {n} entries")
    df = pd.DataFrame(rows)
    df.to_csv(f, index=False)
    return df


def zip_member(url: str, row) -> bytes:
    """One zip member via a single range request (local header + data), inflated in memory."""
    end = int(row.offset) + 30 + 4096 + int(row.csize)
    b = http_get(url, {"Range": f"bytes={int(row.offset)}-{end}"}, expect=206)
    STATE["zenodo_bytes"] += len(b)
    if b[:4] != b"PK\x03\x04":
        raise RuntimeError(f"bad local header for {row.name}")
    fnl, exl = struct.unpack("<HH", b[26:30])
    data = b[30 + fnl + exl:30 + fnl + exl + int(row.csize)]
    if int(row.method) == 8:
        out = zlib.decompress(data, -15)
    elif int(row.method) == 0:
        out = data
    else:
        raise RuntimeError(f"unsupported zip method {row.method}")
    if len(out) != int(row.usize):
        raise RuntimeError(f"size mismatch for {row.name}: {len(out)} vs {row.usize}")
    return out


# ----------------------------------------------------------------------------- metadata
def platemaps(acc: str) -> pd.DataFrame:
    """HepG2 platemaps of all four sites from the gallery (cached under euos_cache/metadata)."""
    frames = []
    for site in SITES:
        d = META_CACHE / "platemaps" / site
        if len(list(d.glob("*__*.txt"))) != 28:
            shutil.rmtree(d, ignore_errors=True)
            d.mkdir(parents=True, exist_ok=True)
            _, keys = s3_list(f"{acc}/{site}/workspace/metadata/platemaps/", delimiter="")
            keys = [k for k, _ in keys if "HepG2" in k and "/platemap/" in k and k.endswith(".txt")]
            if len(keys) != 28:
                raise RuntimeError(f"expected 28 HepG2 platemaps for {site}, found {len(keys)}")
            with ThreadPoolExecutor(8) as ex:
                blobs = list(ex.map(lambda k: http_get(s3_url(k)), keys))
            for k, b in zip(keys, blobs):  # keep the batch folder name: <batch>__<plate>.txt
                (d / f"{k.split('/')[-3]}__{Path(k).name}").write_bytes(b)
        for f in sorted(d.glob("*__*.txt")):
            frames.append(pd.read_csv(f, sep="\t").assign(site=site, batch=f.name.split("__")[0]))
    pm = pd.concat(frames, ignore_index=True)
    pm = pm.rename(columns={"plate_map_name": "plate", "well_position": "well", "Metadata_EOS": "eos",
                            "Metadata_Concentration": "conc", "Metadata_Plate": "layout",
                            "Metadata_Replicate": "rep"})
    assert not pm.duplicated(["site", "plate", "well"]).any()
    if not (pm.groupby("site").plate.nunique().reindex(SITES) == 28).all():
        raise RuntimeError("platemaps: expected 28 HepG2 plates per site")
    return pm[["site", "batch", "plate", "layout", "rep", "well", "eos", "conc"]]


def unit_map() -> tuple[pd.DataFrame, dict]:
    u = pd.read_csv(STEP2 / "units.csv")
    u = u[u.tier == TIER]
    assert len(u) == 172, f"expected 172 {TIER} units, got {len(u)}"
    m = pd.read_csv(ROOT / "cellcheck/cellline_matches.csv")
    m = m[m.dataset == MATCH_DATASET]
    ids_units = {i for s in u.dataset_ids for i in str(s).split(";")}
    ids_match = {i for s in m.dataset_cpd_ids.dropna() for i in str(s).split(";")}
    eos2unit = {}
    for r in u.itertuples():
        for i in str(r.dataset_ids).split(";"):
            assert i not in eos2unit, f"EOS id {i} in two units"
            eos2unit[i] = r.unit
    info = {"n_units": len(u), "n_eos_ids": len(eos2unit),
            "ids_in_units_not_matches": sorted(ids_units - ids_match),
            "ids_in_matches_not_units": sorted(ids_match - ids_units),
            "units_multi_id": int((u.dataset_ids.str.count(";") > 0).sum())}
    return u, eos2unit, info


# ----------------------------------------------------------------------------- download
def site_of(name: str) -> str | None:
    m = re.match(r"Aggregated_Profiles/aggregated_data/([A-Za-z]+)_HepG2/.+\.csv$", name)
    return m.group(1) if m else None


def fetch_plate(url: str, row) -> dict:
    out = PLATE_CACHE / (Path(row.name).stem + ".parquet")
    if out.exists():
        return {"file": row.name, "status": "cached"}
    t0 = time.time()
    raw = zip_member(url, row)
    df = pd.read_csv(io.BytesIO(raw))
    meta = [c for c in df.columns if c.startswith("Metadata_")]
    feats = [c for c in df.columns if not c.startswith("Metadata_")]
    df[feats] = df[feats].astype(np.float32)
    tmp = out.with_suffix(".tmp")
    df.to_parquet(tmp, index=False)
    tmp.replace(out)
    return {"file": row.name, "status": "downloaded", "bytes": int(row.csize), "seconds": round(time.time() - t0, 1),
            "wells": len(df), "n_meta": len(meta), "n_feat": len(feats)}


def download(zrow_url: str, members: pd.DataFrame) -> list[dict]:
    PLATE_CACHE.mkdir(parents=True, exist_ok=True)
    recs = []
    with ThreadPoolExecutor(N_THREADS) as ex:
        futs = [ex.submit(fetch_plate, zrow_url, r) for r in members.itertuples()]
        for i, f in enumerate(as_completed(futs), 1):
            recs.append(f.result())  # raises loudly on any failure
            track_disk()
            if i % 10 == 0 or i == len(futs):
                print(f"  {i}/{len(futs)} plate files, {STATE['zenodo_bytes'] / 1e9:.2f} GB fetched", flush=True)
    return recs


# ----------------------------------------------------------------------------- plate identity checks
def plate_label(d: pd.DataFrame, fname: str) -> str:
    """Canonical plate-map name (e.g. B1001_R1). FMP files store Metadata_Plate='B1001' and the
    replicate in Metadata_Batch='R1'; the other sites store 'B1001_R1' in Metadata_Plate."""
    p = d.Metadata_Plate.astype(str).unique()
    b = d.Metadata_Batch.astype(str).unique()
    if len(p) != 1 or len(b) != 1:
        raise RuntimeError(f"{fname}: several plate/batch labels {p} {b}")
    p, b = p[0], b[0]
    if re.fullmatch(r"B\d{4}_R\d", p):
        lab = p
    elif re.fullmatch(r"B\d{4}", p) and re.fullmatch(r"R\d", b):
        lab = f"{p}_{b}"
    else:
        raise RuntimeError(f"{fname}: cannot parse plate label from {p!r}/{b!r}")
    lay, rep = lab.split("_")
    if lay not in fname or rep not in fname:
        raise RuntimeError(f"{fname}: in-file label {lab} disagrees with the file name")
    return lab


def resolve_duplicates(members: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Several files per (site, Metadata_Plate) (IMTM B1001_R4 / B1007_R4 swap, Zenodo README).
    Keep, for each duplicated label, the file whose per-well log object count correlates best with
    the other replicates of the same layout; assert it is also the better fit vs. the alternative layout."""
    rows = []
    for r in members.itertuples():
        d = pd.read_parquet(PLATE_CACHE / (Path(r.name).stem + ".parquet"),
                            columns=["Metadata_Batch", "Metadata_Plate", "Metadata_Well", "Metadata_Object_Count"])
        rows.append(dict(name=r.name, site=r.site, plate=plate_label(d, Path(r.name).name),
                         oc=d.set_index("Metadata_Well").Metadata_Object_Count))
    tab = pd.DataFrame(rows)
    # per-well mean log count over replicates, per (site, layout), using non-duplicated files only
    tab["layout"] = tab.plate.str.split("_").str[0]
    dup = tab.duplicated(["site", "plate"], keep=False)
    ref = {}
    for (site, lay), g in tab[~dup].groupby(["site", "layout"]):
        ref[(site, lay)] = {n: np.log1p(x) for n, x in zip(g.name, g.oc)}

    def fit(oc, site, lay, exclude=None):
        a = np.log1p(oc)
        b = pd.concat([v for n, v in ref[(site, lay)].items() if n != exclude], axis=1).mean(1).reindex(a.index)
        ok = a.notna() & b.notna()
        return float(np.corrcoef(a[ok], b[ok])[0, 1])

    info = {"files_listed": len(tab), "duplicates": {}}
    # general label check on every plate: own layout vs best other layout (same site)
    checks = []
    for r in tab.itertuples():
        own = fit(r.oc, r.site, r.layout, exclude=r.name)
        others = {lay: fit(r.oc, r.site, lay) for (s, lay) in ref if s == r.site and lay != r.layout}
        best_other = max(others, key=others.get)
        checks.append(dict(name=r.name, site=r.site, plate=r.plate, r_own=round(own, 3),
                           best_other=best_other, r_other=round(others[best_other], 3)))
    chk = pd.DataFrame(checks)
    keep = np.ones(len(tab), bool)
    for (site, plate), g in tab[dup].groupby(["site", "plate"]):
        c = chk.loc[g.index]
        winner = c.r_own.idxmax()
        for i in g.index:
            keep[i] = i == winner
        info["duplicates"][f"{site}:{plate}"] = c[["name", "r_own", "best_other", "r_other"]].to_dict("records")
        w = c.loc[winner]
        if not w.r_own > w.r_other:
            raise RuntimeError(f"duplicate {site}:{plate}: chosen file does not fit its layout best: {w.to_dict()}")
    bad = chk[keep & (chk.r_own <= chk.r_other)]
    info["mislabel_suspects"] = bad.to_dict("records")
    info["label_check_r_own_median"] = float(chk.r_own[keep].median())
    info["label_check_r_other_median"] = float(chk.r_other[keep].median())
    info["label_check_min_margin"] = float((chk.r_own - chk.r_other)[keep].min())
    if len(bad):
        raise RuntimeError(f"plates whose object-count pattern fits another layout better: {bad.to_dict('records')}")
    kept = members.loc[tab.index[keep]].copy()
    kept["plate"] = tab.plate[keep].values
    info["dropped_files"] = tab.name[~keep].tolist()
    return kept, info


# ----------------------------------------------------------------------------- loading
def load_wells(kept: pd.DataFrame, pm: pd.DataFrame):
    frames, featsets = [], set()
    for r in kept.itertuples():
        d = pd.read_parquet(PLATE_CACHE / (Path(r.name).stem + ".parquet"))
        assert plate_label(d, Path(r.name).name) == r.plate
        d["Metadata_Plate"] = r.plate
        d = d.rename(columns={"Metadata_Plate": "plate", "Metadata_Well": "well",
                              "Metadata_Object_Count": "object_count"})
        d.insert(0, "site", r.site)
        frames.append(d)
        featsets.add(tuple(c for c in d.columns if not c.startswith("Metadata_")
                           and c not in ("site", "plate", "well", "object_count")))
    if len(featsets) != 1:
        raise RuntimeError(f"feature sets differ between plate files ({len(featsets)} variants)")
    feats = list(next(iter(featsets)))
    prof = pd.concat(frames, ignore_index=True)
    assert not prof.duplicated(["site", "plate", "well"]).any()
    obs = pm.merge(prof[["site", "plate", "well", "object_count"]].assign(_row=np.arange(len(prof))),
                   on=["site", "plate", "well"], how="left", indicator=True)
    extra = prof.merge(pm[["site", "plate", "well"]], how="left", indicator=True)
    info = {"wells_expected": len(pm), "wells_in_files": len(prof),
            "wells_missing": obs[obs._merge == "left_only"][["site", "plate", "well", "eos"]].to_dict("records"),
            "wells_not_in_platemap": int((extra._merge == "left_only").sum())}
    if info["wells_not_in_platemap"]:
        raise RuntimeError(f"{info['wells_not_in_platemap']} profile wells are not in the platemaps")
    obs = obs[obs._merge == "both"].drop(columns="_merge").reset_index(drop=True)
    X = prof[feats].to_numpy(np.float32)[obs._row.to_numpy().astype(np.int64)]
    obs = obs.drop(columns="_row")
    return obs, X, feats, info


# ----------------------------------------------------------------------------- processing
def plate_groups(obs):
    out = []
    dm = obs.is_dmso.to_numpy()
    for _, idx in obs.groupby(["site", "plate"]).indices.items():
        idx = np.sort(idx)
        out.append((idx, idx[dm[idx]]))
    return out


def normalize(obs, X, feats, info):
    f = np.array(feats)
    book = np.array([bool(BOOKKEEPING.search(c)) for c in feats])
    info["drop_bookkeeping"] = int(book.sum())
    info["bookkeeping_features"] = f[book].tolist()
    X = X[:, ~book]
    f = f[~book]
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
    groups = plate_groups(obs)
    meds = np.empty((len(groups), X.shape[1]))
    scales = np.empty_like(meds)
    for g, (_, didx) in enumerate(groups):
        assert len(didx) >= 8, f"too few DMSO wells on a plate ({len(didx)})"
        D = X[didx].astype(np.float64)
        meds[g] = np.nanmedian(D, axis=0)
        scales[g] = MAD_SCALE * np.nanmedian(np.abs(D - meds[g]), axis=0)
    ref = np.nanmedian(scales, axis=0)
    floor = np.maximum(MAD_FLOOR_REL * np.nan_to_num(ref), MAD_FLOOR_ABS)
    floored = scales < floor[None]
    info["plate_feature_scales_floored_frac"] = round(float(np.nanmean(floored)), 4)
    scales = np.where(np.isfinite(scales) & ~floored, scales, floor[None])
    Z = np.empty_like(X)
    n_imp = 0
    for g, (idx, _) in enumerate(groups):
        x = X[idx].astype(np.float64)
        nan = np.isnan(x)
        n_imp += int(nan.sum())
        x = np.where(nan, meds[g][None], x)
        z = (x - meds[g][None]) / scales[g][None]
        z[~np.isfinite(z)] = 0.0
        Z[idx] = z.astype(np.float32)
    info["n_values_imputed"] = n_imp
    maxabs = np.abs(Z).max(0)
    np.clip(Z, -CLIP_Z, CLIP_Z, out=Z)
    return Z, maxabs, f.tolist()


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


def correlation_prune(Z, rng):
    idx = np.sort(rng.choice(len(Z), min(CORR_SUBSAMPLE, len(Z)), replace=False))
    R = np.abs(np.corrcoef(Z[idx].astype(np.float64), rowvar=False))
    R[~np.isfinite(R)] = 0
    order = np.argsort(R.sum(0), kind="stable")
    keep = np.zeros(Z.shape[1], bool)
    kept: list[int] = []
    for j in order:
        if not kept or R[j, kept].max() <= CORR_THRESH:
            keep[j] = True
            kept.append(j)
    return keep


def select_features(Z, maxabs, feats, info, rng):
    f = np.array(feats)
    blk = np.isin(f, BLOCKLIST)
    info["drop_blocklist"] = int(blk.sum())
    out = (~blk) & (maxabs > OUTLIER_Z)
    info["drop_outlier_gt500"] = int(out.sum())
    keep = ~(blk | out)
    vd = np.zeros(len(f), bool)
    vd[keep] = variance_threshold_drop(Z[:, keep])
    info["drop_low_variance"] = int(vd.sum())
    keep &= ~vd
    ck = np.zeros(len(f), bool)
    ck[keep] = correlation_prune(Z[:, keep], rng)
    info["drop_correlated"] = int((keep & ~ck).sum())
    keep &= ck
    info["n_selected"] = int(keep.sum())
    return Z[:, keep], f[keep].tolist()


def fit_sphere(D, info):
    mu = D.mean(0)
    C = np.cov((D - mu).astype(np.float64), rowvar=False)
    lam, U = np.linalg.eigh(C)
    lam, U = lam[::-1].clip(min=0), U[:, ::-1]
    k = int(np.searchsorted(np.cumsum(lam) / lam.sum(), VAR_KEEP) + 1)
    info["sphere_dim"] = k
    info["sphere_eig_first_last"] = [round(float(lam[0]), 3), round(float(lam[k - 1]), 4)]
    info["sphere_eps"] = SPHERE_EPS
    return mu, U[:, :k] / np.sqrt(lam[:k] + SPHERE_EPS)


def unitize(V):
    n = np.linalg.norm(V, axis=-1, keepdims=True)
    return V / np.where(n > 0, n, 1)


# ----------------------------------------------------------------------------- metrics
def ap_rows(S, pos):
    """Average precision for each row. S: (..., m) scores, pos: (..., m) bool (positives)."""
    order = np.argsort(-S, axis=-1, kind="stable")
    rel = np.take_along_axis(pos, order, axis=-1)
    cum = np.cumsum(rel, axis=-1)
    ranks = np.arange(1, S.shape[-1] + 1)
    prec = np.where(rel, cum / ranks, 0.0)
    npos = rel.sum(-1)
    return prec.sum(-1) / np.maximum(npos, 1)


def rep_map(U, rep_idx, neg_idx):
    """Each replicate well queries {other replicates} + negatives, ranked by cosine. mAP over queries."""
    pool = np.concatenate([rep_idx, neg_idx])
    G = U[rep_idx] @ U[pool].T                      # (n, n + m)
    n = len(rep_idx)
    pos = np.zeros_like(G, bool)
    pos[:, :n] = True
    keep = np.ones_like(G, bool)
    keep[np.arange(n), np.arange(n)] = False       # drop self
    S = G[keep].reshape(n, -1)
    P = pos[keep].reshape(n, -1)
    return float(ap_rows(S, P).mean())


def rows_by(obs, mask, key):
    """{group: positional row indices into obs} for the rows selected by mask."""
    pos = np.where(np.asarray(mask))[0]
    return {k: pos[v] for k, v in obs.iloc[pos].groupby(key).indices.items()}


def dmso_by_plate(obs):
    dm = obs.is_dmso.to_numpy()
    return {k: np.sort(v[dm[v]]) for k, v in obs.groupby(["site", "plate"]).indices.items()}


class MapNull:
    """Replicate-mAP null: same plate structure, replicates replaced by DMSO wells (one per replicate
    well, drawn from the replicate's plate without replacement); negatives = all other DMSO wells on
    those plates. Cached by the multiset of plates."""

    def __init__(self, U, dbp, rng, n_draws):
        self.U, self.dbp, self.rng, self.n = U, dbp, rng, n_draws
        self.cache = {}

    def get(self, plates):
        key = tuple(sorted(plates))
        if key in self.cache:
            return self.cache[key]
        uniq, cnt = np.unique(np.array([f"{s}|{p}" for s, p in key]), return_counts=True)
        all_d = np.concatenate([self.dbp[tuple(u.split("|"))] for u in uniq])
        G = self.U[all_d] @ self.U[all_d].T
        loc = {}
        off = 0
        for u in uniq:
            m = len(self.dbp[tuple(u.split("|"))])
            loc[u] = np.arange(off, off + m)
            off += m
        M = len(all_d)
        n = int(cnt.sum())
        draws = np.empty((self.n, n), int)
        for t in range(self.n):
            draws[t] = np.concatenate([self.rng.choice(loc[u], c, replace=False) for u, c in zip(uniq, cnt)])
        vals = np.empty(self.n)
        for t in range(self.n):
            q = draws[t]
            S = G[q]                                   # (n, M)
            P = np.zeros((n, M), bool)
            P[:, q] = True
            keep = np.ones((n, M), bool)
            keep[np.arange(n), q] = False
            vals[t] = ap_rows(S[keep].reshape(n, M - 1), P[keep].reshape(n, M - 1)).mean()
        self.cache[key] = vals
        return vals


def pval(null, x):
    return float((1 + np.sum(null >= x)) / (1 + len(null)))


def consensus(V, sites):
    """median within site, then mean across sites."""
    s = np.asarray(sites)
    return np.mean([np.median(V[s == k], axis=0) for k in sorted(set(s))], axis=0)


def mean_pair_cos(U):
    G = U @ U.T
    n = len(U)
    return float((G.sum() - np.trace(G)) / (n * (n - 1)))


# ----------------------------------------------------------------------------- notes
def write_notes(run: dict, act: pd.DataFrame, ctrl: pd.DataFrame) -> None:
    from importlib.metadata import version
    vers = {p: version(p) for p in ("pandas", "numpy", "pyarrow", "requests")}
    z, f, sm, ps, cs = run["zenodo"], run["features"], run["summary"], run["position"], run["cross_site"]
    zf = next(x for x in z["files"] if x["key"] == ZIP_NAME)
    di = run["plate_identity"]
    miss = pd.DataFrame(run["load"]["wells_missing"])
    unit_eos = set(act.eos_ids.str.split(";").explode())
    dl = run["download_summary"]
    L: list[str] = []
    a = L.append
    a("# EU-OS T2 morphology profiles (HepG2, 24 h, 10 uM): notes (written by step2/morph_euos.py)\n")
    a("## Sources")
    a(f"- Gallery accession (resolved): `{run['gallery']['accession']}`. Per-site workspace prefixes: "
      + "; ".join(f"{s}: {', '.join(x.rstrip('/').split('/')[-1] for x in run['gallery'][s])}" for s in SITES)
      + ". There is **no `workspace/profiles` or `backend` prefix** (only images, load_data_csv, metadata; plus "
        "`all/workspace/metadata_harmonized/` with one 0.88 GB draft metadata CSV, not used).")
    a(f"- Profiles: Zenodo concept {z['concept']} -> latest record **{z['record_id']} (v{z['version']}, "
      f"{z['publication_date']})**, DOI {z['doi']}. File `{ZIP_NAME}` ({zf['size']:,} bytes): {zf['url']}")
    a("- Versions of the record (REST API): " + "; ".join(
        f"{v['id']} v{v['version']} ({v['date']}): {ZIP_NAME} {v['files'].get(ZIP_NAME, 0):,} B"
        for v in run["zenodo_versions"]) + ". Only v3.0.0 carries the IMTM plate-swap fix.")
    a("- The zip is never downloaded whole: its central directory is read with two HTTP range requests and "
      "each HepG2 plate CSV member is fetched with one range request, inflated in memory and stored as "
      "float32 parquet in `step2/euos_cache/plates/` (deleted at the end of the run).")
    a(f"- Zip members used: {run['zip_members_hepg2']} HepG2 CSVs (the 28 `FMP_U2OS/` members are skipped). "
      "Platemaps: gallery `<acc>/<site>/workspace/metadata/platemaps/*HepG2/platemap/*.txt` (cached in "
      "`step2/euos_cache/metadata/platemaps/`).")
    a("\n## File format and quirks")
    a(f"- One CSV per plate, one row per well (CellProfiler features aggregated per well by the authors): "
      f"`Metadata_Batch, Metadata_Plate, Metadata_Well, Metadata_Object_Count` + {f['n_features_raw']} features "
      "(prefixes `Nuc_`, `Cells_`, `Cyto_`). The feature set is **identical at all four sites** (same 2981-column "
      "header), so no intersection was needed.")
    a("- FMP files store the plate as `Metadata_Plate='B1001'` with the replicate in `Metadata_Batch='R1'`; the other "
      "sites store `B1001_R1`. The script builds the canonical `B1001_R1` and checks it against the file name.")
    a(f"- IMTM plate swap (Zenodo README + `2026-02-19_SwappedID_IMTM.pdf`): IMTM has 29 files and `B1001_R4` appears "
      f"twice. Resolved by object-count pattern (per-well log count vs. the leave-one-out mean of the other replicates "
      f"of each layout): {json.dumps(di['duplicates'])}. Dropped: {di['dropped_files']} (the stale pre-fix file named "
      f"in the README). Label check on all 112 kept plates: per-well log-count correlation with own layout median "
      f"{di['label_check_r_own_median']:.3f} vs best other layout {di['label_check_r_other_median']:.3f}; smallest "
      f"margin {di['label_check_min_margin']:.3f}; no suspects.")
    a(f"- Wells: {run['load']['wells_expected']} in platemaps, {run['load']['wells_in_files']} in files, "
      f"{len(miss)} platemap wells absent from the files (cause not recorded; plausibly wells without usable cells "
      f"or images): "
      + (", ".join(f"{k} {v}" for k, v in miss.groupby('site').size().items()) if len(miss) else "none")
      + f"; of these, matched-unit wells: {int(miss.eos.isin(unit_eos).sum()) if len(miss) else 0}. Units with fewer "
        "than 16 wells: " + "; ".join(f"{r.unit} ({r.eos_ids}, {r.n_wells} wells)" for r in
                                     act[(act.n_wells < 16)].itertuples()) + ".")
    a("- Layout: DMSO = 28 wells per plate, all in the edge columns 23/24; compounds in columns 1-22; "
      "nocodazole/tetrandrine 5 uM in columns 23/24 (positions differ by site). Within a site the layout of "
      "B100x_R1..R4 is identical; across sites only the control positions differ. Negative control = "
      "`Metadata_EOS == 'DMSO'` (never concentration 0).")
    a("- The 'concentration 0' control wells (14 nocodazole + 14 tetrandrine at each of MEDINA and USC) behave like "
      "the 5 uM controls (mean cosine 0.71-0.75 to the 5 uM consensus of the same control and site, similar norms; "
      "table below), i.e. a metadata quirk. They are reported, never used.")
    a("\n## Processing (MORPH_SPEC steps 1-7)")
    a(f"- Features: raw {f['n_features_raw']} -> bookkeeping -{f['drop_bookkeeping']} (BoundingBox min/max, "
      f"Center_X/Y, closest-object numbers; deviation, see below) -> NaN/inf >1% -{f['drop_nan_gt1pct']}, constant "
      f"-{f['drop_constant']} = {f['n_after_nan_const']} -> blocklist -{f['drop_blocklist']} (pycytominer list "
      f"translated to `Nuc_`), outliers (pre-clip max|z| > {OUTLIER_Z:g}) -{f['drop_outlier_gt500']}, low variance "
      f"(freq_cut {FREQ_CUT}, unique_cut {UNIQUE_CUT}) -{f['drop_low_variance']}, correlated (greedy |r| > "
      f"{CORR_THRESH}, r estimated on {CORR_SUBSAMPLE} random wells) -{f['drop_correlated']} = "
      f"**{f['n_selected']} selected**.")
    a(f"- Imputed values: {f['n_values_imputed']}; fraction of plate x feature MADs at the floor: "
      f"{f['plate_feature_scales_floored_frac']}.")
    a(f"- Normalization: robust z per plate (site x plate) against its 28 DMSO wells, MAD floor = {MAD_FLOOR_REL} x "
      f"median plate scale, clip +-{CLIP_Z:g}.")
    a(f"- Sphering: PCA-whitening fit on all {sum(v['dmso_wells'] for v in run['per_site'].values())} DMSO wells "
      f"pooled over the four sites, 99% variance, eps {f['sphere_eps']} (same as T3): **{f['sphere_dim']} dimensions** "
      f"(eigenvalues {f['sphere_eig_first_last'][0]} ... {f['sphere_eig_first_last'][1]}).")
    a("- **Site handling:** one feature selection and one sphering for all four sites (the feature sets are "
      "identical, so nothing is lost to an intersection). Consensus = median of wells within each site, then mean "
      "across sites, in that single sphered space. Why: a cross-site mean is only meaningful if the sites share one "
      "space; per-plate robust z against each plate's own DMSO reduces (does not remove) site and plate offsets before "
      f"the joint selection and sphering: replicate-pair cosine is {cs['rep_cos_within_site']} within a site vs "
      f"{cs['rep_cos_across_sites']} across sites. Per-site reproducibility is reported too (within-site mAP in the "
      "shared space). Separate per-site sphering would make the cross-site mean meaningless, so it was not used.")
    a("- Per site: " + "; ".join(f"{s}: {v['plates']} plates, {v['wells']} wells ({v['dmso_wells']} DMSO), "
                                  f"{v['features_raw']} raw features, sphered DMSO mean variance "
                                  f"{v['dmso_sphered_var_mean']}" for s, v in run["per_site"].items()))
    a(f"- Replicate mAP: each replicate well queries its other replicates + all DMSO wells of the replicate plates "
      f"(cosine, sphered). Null: the replicates are replaced by DMSO wells drawn from the same plates (one per "
      f"replicate well), {N_NULL} draws, cached per plate set; reproducible if p < {P_ACTIVE}. Overall = all wells of "
      f"the unit (all sites); per site = that site's 4 wells.")
    a(f"- active_z: robust z of the consensus L2 norm vs all {sm['n_cpds_ref']} treated compounds (median "
      f"{sm['norm_ref_median_all_cpds']}, 1.4826 MAD {sm['norm_ref_mad_all_cpds']}), active if z > {Z_ACTIVE:g}.")
    a(f"- Addendum calls: `active_norm_dmso` (consensus norm vs DMSO pseudo-consensus: one DMSO well per replicate well "
      f"from the same plate, median within site, mean across sites, {N_NULL} draws); `active_vs_cpd` (mean pairwise "
      f"cosine of all replicate wells vs the same score for other compounds' wells drawn from the same plates, distinct "
      f"compounds within a draw); `active_vs_pos` (extra, position-matched: mean within-site replicate-pair cosine vs "
      f"wells at the SAME well position, site and replicate index on plates of other layouts, i.e. different compounds "
      f"at the same position, distinct layouts within a site); `active_cpd_and_pos` = active_vs_cpd AND active_vs_pos. "
      "What each null matches: `active_vs_cpd` matches plates (hence site and batch) but not well position; "
      "`active_vs_pos` matches site, replicate index and well position but not batch where a layout's replicates "
      "were imaged together (see batch sharing below).")
    a("\n## Counts (dose 10 uM, HepG2, 24 h)")
    a(f"- Units: **{sm['n_units']} / 172** covered (EOS ids {run['unit_map']['n_eos_ids']}; one unit has two EOS ids). "
      f"Wells per unit: {sm['n_wells_dist']} (16 = 4 sites x 4 replicates; fewer = wells absent from the files; "
      "32 = two EOS ids).")
    a(f"- All wells: {json.dumps(run['roles'])}; all 2464 treated compounds are used for normalization context, "
      "selection and the active_z reference.")
    a("\n## Reproducibility and activity (172 units)")
    a(f"- **Headline: site effects dominate the replicate structure.** Replicate-pair cosine {cs['rep_cos_within_site']} "
      f"within a site vs {cs['rep_cos_across_sites']} across sites; median replicate mAP {sm['rep_mAP_median']} over all "
      f"16 wells vs " + ", ".join(f"{s} {sm['rep_mAP_median_per_site'][s]}" for s in SITES) + " within site. "
      f"Same unit, consensus at two different sites: cosine {cs['site_consensus_cos_same_unit']} (by site pair "
      f"{cs['site_consensus_cos_same_unit_by_pair']}) vs {cs['site_consensus_cos_diff_unit_diff_site']} for different "
      "units at different sites, so every site pair carries shared compound signal.")
    a(f"- Replicate mAP reproducible (p<0.05), overall: {sm['reproducible_overall']:.1%} (DMSO-null mean mAP "
      f"{sm['mapnull_mean_overall_example']}); saturated, see caveat.")
    a("- Per site (within-site mAP): " + ", ".join(f"{s} {sm['reproducible_per_site'][s]:.1%} (median mAP "
                                                    f"{sm['rep_mAP_median_per_site'][s]})" for s in SITES)
      + f". Number of sites reproducible per unit: {sm['n_sites_reproducible_dist']}.")
    a(f"- active_z {sm['active_z']:.1%}; active_norm_dmso {sm['active_norm_dmso']:.1%}; active_vs_cpd "
      f"{sm['active_vs_cpd']:.1%}; active_vs_pos {sm['active_vs_pos']:.1%} (per site, within-site score: "
      + ", ".join(f"{k} {v:.1%}" for k, v in sm['active_vs_pos_per_site'].items())
      + f"); active_cpd_and_pos {sm['active_cpd_and_pos']:.1%}.")
    a(f"- Batch sharing of within-site replicate pairs (platemap batch folders): "
      f"{ps['frac_within_site_replicate_pairs_same_batch']}. At FMP batch = replicate index, so `active_vs_pos` is "
      "batch-matched there; at IMTM a layout's 4 replicates were imaged in one batch, so replicate pairs share a batch "
      "and null pairs do not. Compare the FMP and IMTM per-site rates above to judge the batch contribution.")
    a("- **Saturation caveat:** `active_rep` and `active_norm_dmso` are saturated and should not be used as activity "
      "calls here. DMSO sits only in edge columns 23/24 while compounds sit in columns 1-22, and all replicates of a "
      f"compound share one well position at every site (replicate pairs at the same well: "
      f"{ps['frac_replicate_pairs_same_well']:.0%}), so any compound set differs from DMSO by an edge/position "
      f"offset. DMSO wells at the same position on different plates: cosine {ps['dmso_cos_same_well_diff_plate']} vs "
      f"{ps['dmso_cos_diff_well_diff_plate']} at different positions; different compounds at the same position: "
      f"{ps['cpd_cos_diff_cpd_same_well']} vs {ps['cpd_cos_diff_cpd_diff_well']} at different positions. Neither "
      "addendum null controls everything (position vs batch, see above); `active_z` is the most conservative call.")
    a("\n## Positive controls (5 uM; excluded from the compound set)")
    a("```\n" + ctrl.to_string(index=False) + "\n```")
    a("\n## Outputs")
    a(f"- `step2/morph_euos.parquet`: 172 units + a `DMSO` reference row (dose_uM 0), f0..f{f['sphere_dim'] - 1} float32.")
    a("- `step2/morph_euos_activity.csv`: spec columns + norm, norm_dmso_p, active_norm_dmso, rep_cos, vs_cpd_p, "
      "active_vs_cpd, rep_cos_within_site, vs_pos_p, active_vs_pos, active_cpd_and_pos, vs_pos_p_<site>, "
      "per-site rep_mAP/rep_p, "
      "n_sites_reproducible, eos_ids, wells.")
    a("- `step2/morph_euos_controls.csv` (positive controls); `step2/morph_euos_wells.parquet` (sphered wells of matched "
      "units and DMSO, with site/batch/plate/layout/rep/well/eos/unit/role/conc/object_count); "
      "`step2/euos_cache/euos_run.json` (full run info incl. selected feature names); `step2/euos_cache/metadata/` "
      "(Zenodo record JSON, zip directory, platemaps).")
    a("\n## Bytes, disk, time, versions")
    a(f"- This run: {dl['bytes']:,} bytes over HTTP (plate members {dl['zenodo_member_bytes_this_run']:,}); peak disk "
      f"{dl['peak_disk_bytes']:,} bytes (filtered float32 plate cache {dl['plate_cache_bytes']:,} bytes, deleted: "
      f"{dl['plate_cache_deleted']}); raw zip/CSV bytes never touch the disk. Wall time {dl['seconds']} s "
      f"(download {run['download']['seconds']} s).")
    a("- Cumulative over development: two full plate fetches (about 2 x 1.08 GB) plus about 13 MB of manual probes "
      "(zip directory listings, one IMTM CSV, the swap PDF, header range reads): about 2.2 GB in total.")
    a(f"- Python 3.13; packages: {vers}.")
    a("\n## Deviations from MORPH_SPEC")
    a("- Profiles come from Zenodo (not the gallery), fetched by range requests from one zip, because the gallery "
      "holds no profiles.")
    a(f"- {f['drop_bookkeeping']} bookkeeping columns (object bounding-box coordinates, centers, closest-object "
      "numbers) are dropped before the NaN/constant step: in these raw per-well aggregates they are not morphology.")
    a("- Feature correlation for the |r| > 0.9 pruning is estimated on a 20,000-well random subsample (memory); all "
      "other steps use all wells.")
    a("- Blocklist names translated from `Nuclei_` to `Nuc_` (this dataset's compartment name).")
    a("- Extra activity calls `active_vs_pos` (position-matched) and `active_cpd_and_pos`, because every replicate shares "
      "a well position and DMSO is confined to the edge columns, which saturates the DMSO-based calls.")
    a("- The two IMTM `B1001_R4` files are resolved from the data (object-count pattern); the result agrees with the "
      "Zenodo README.")
    OUT_NOTES.write_text("\n".join(L) + "\n", encoding="utf-8")


# ----------------------------------------------------------------------------- main
def main():
    t_start = time.time()
    rng = np.random.default_rng(SEED)
    CACHE.mkdir(exist_ok=True)
    META_CACHE.mkdir(parents=True, exist_ok=True)
    run: dict = {}

    # --- sources
    acc = resolve("euos")
    prefixes, _ = s3_list(f"{acc}/", "/")
    gallery = {"accession": acc, "top_prefixes": prefixes}
    for site in SITES:
        p, _ = s3_list(f"{acc}/{site}/workspace/", "/")
        gallery[site] = p
        if any(x.rstrip("/").endswith(("profiles", "backend")) for x in p):
            log(f"NOTE: gallery has a profiles/backend prefix for {site}: {p}")
    has_prof = any(x.rstrip("/").endswith(("profiles", "backend"))
                   for s in SITES for x in gallery[s])
    log(f"gallery {acc}: profiles/backend prefixes present: {has_prof}")
    run["gallery"] = gallery

    zen = zenodo_latest()
    zfile = next((x for x in zen["files"] if x["key"] == ZIP_NAME), None)
    if zfile is None:
        raise RuntimeError(f"{ZIP_NAME} not in Zenodo record {zen['record_id']}")
    log(f"Zenodo record {zen['record_id']} v{zen['version']} ({zen['publication_date']}): {ZIP_NAME} "
        f"{zfile['size'] / 1e9:.3f} GB at {zfile['url']}")
    zdir = zip_directory(zfile["url"], zfile["size"])
    zdir["site"] = zdir.name.map(site_of)
    members = zdir[zdir.site.isin(SITES)].reset_index(drop=True)
    log(f"zip: {len(zdir)} entries; HepG2 CSV members: {members.groupby('site').size().to_dict()} "
        f"({members.csize.sum() / 1e9:.3f} GB compressed)")
    run["zip_members_hepg2"] = members.groupby("site").size().to_dict()

    pm = platemaps(acc)
    track_disk()
    log(f"platemaps: {pm.groupby('site').plate.nunique().to_dict()} plates, {len(pm)} wells")
    units, eos2unit, uinfo = unit_map()
    run["unit_map"] = uinfo
    log(f"units: {uinfo}")

    # --- download (range requests; nothing but filtered parquet ever lands on disk)
    t0 = time.time()
    recs = download(zfile["url"], members)
    t_dl = time.time() - t0
    new = [r for r in recs if r["status"] == "downloaded"]
    log(f"plate files: {len(new)} downloaded ({sum(r['bytes'] for r in new) / 1e9:.3f} GB), "
        f"{len(recs) - len(new)} cached, {t_dl:.0f}s")
    run["download"] = {"n_new": len(new), "bytes_new": sum(r["bytes"] for r in new), "seconds": round(t_dl, 1)}

    kept, dinfo = resolve_duplicates(members)
    run["plate_identity"] = dinfo
    log(f"duplicate resolution: dropped {dinfo['dropped_files']}; label check r_own median "
        f"{dinfo['label_check_r_own_median']:.3f} vs best other layout {dinfo['label_check_r_other_median']:.3f}, "
        f"min margin {dinfo['label_check_min_margin']:.3f}")
    cnt = kept.groupby("site").size()
    if not (cnt.reindex(SITES) == 28).all():
        raise RuntimeError(f"expected 28 plates per site after dedup, got {cnt.to_dict()}")

    # --- load + annotate
    obs, X, feats, linfo = load_wells(kept, pm)
    run["load"] = {k: v for k, v in linfo.items()}
    log(f"wells: {len(obs)} loaded of {linfo['wells_expected']} in platemaps; missing "
        f"{len(linfo['wells_missing'])}: {linfo['wells_missing']}")
    obs["is_dmso"] = obs.eos == DMSO
    obs["is_ctrl"] = obs.eos.isin(CONTROLS)
    obs["is_cpd"] = obs.eos.str.startswith("EOS")
    assert (obs.is_dmso | obs.is_ctrl | obs.is_cpd).all()
    assert (obs.conc[obs.is_cpd] == DOSE_UM).all()
    obs["unit"] = obs.eos.map(eos2unit)
    obs["role"] = np.select([obs.is_dmso, obs.is_ctrl, obs.unit.notna()],
                            ["negative_control", "positive_control", "matched_compound"], "other_compound")
    miss_ids = sorted(set(eos2unit) - set(obs.eos))
    if miss_ids:
        raise RuntimeError(f"unit EOS ids absent from HepG2 platemaps: {miss_ids}")
    run["roles"] = {k: int(v) for k, v in obs.role.value_counts().items()}
    log(f"roles: {run['roles']}; {obs.eos[obs.is_cpd].nunique()} compounds")

    # --- normalize, select, sphere
    t0 = time.time()
    finfo: dict = {"n_features_raw": len(feats)}
    Z, maxabs, f1 = normalize(obs, X, feats, finfo)
    del X
    Zs, fsel = select_features(Z, maxabs, f1, finfo, rng)
    del Z
    finfo["selected_features"] = fsel
    D = Zs[obs.is_dmso.to_numpy()]
    mu, W = fit_sphere(D, finfo)
    V = ((Zs - mu) @ W).astype(np.float32)
    Vs = ((D - mu) @ W)
    finfo["dmso_sphered_var_mean"] = round(float(Vs.var(0).mean()), 3)
    run["features"] = finfo
    log(f"features: raw {finfo['n_features_raw']} -> bookkeeping -{finfo['drop_bookkeeping']}, NaN>1% "
        f"-{finfo['drop_nan_gt1pct']}, constant -{finfo['drop_constant']} = {finfo['n_after_nan_const']} -> "
        f"blocklist -{finfo['drop_blocklist']}, outlier -{finfo['drop_outlier_gt500']}, low variance "
        f"-{finfo['drop_low_variance']}, correlated -{finfo['drop_correlated']} = {finfo['n_selected']} selected; "
        f"sphered dims {finfo['sphere_dim']} ({time.time() - t0:.0f}s)")
    # per-site feature diagnostics (one shared space; report per-site DMSO spread)
    site_arr = obs.site.to_numpy()
    run["per_site"] = {s: {"plates": int(obs.plate[obs.site == s].nunique()),
                           "wells": int((obs.site == s).sum()),
                           "dmso_wells": int((obs.is_dmso & (obs.site == s)).sum()),
                           "features_raw": finfo["n_features_raw"],
                           "dmso_sphered_var_mean": round(float(V[(obs.is_dmso & (obs.site == s)).to_numpy()].var(0).mean()), 3)}
                       for s in SITES}

    U = unitize(V.astype(np.float64))
    dbp = dmso_by_plate(obs)
    mapnull = MapNull(U, dbp, rng, N_NULL)
    plates_key = list(zip(obs.site, obs.plate))

    # --- position diagnostics
    pos = {}
    g = obs[obs.is_cpd].groupby("eos")
    pos["cpd_wells_per_position_within_site"] = float(obs[obs.is_cpd].groupby(["eos", "site"]).well.nunique().groupby("eos").max().mean())
    pos["frac_cpds_same_well_all_sites"] = float((g.well.nunique() == 1).mean())
    pairs_same = pairs_all = 0
    for _, d in g:
        w = d.well.to_numpy()
        n = len(w)
        pairs_all += n * (n - 1) // 2
        pairs_same += int(((w[:, None] == w[None, :]).sum() - n) // 2)
    pos["frac_replicate_pairs_same_well"] = pairs_same / pairs_all
    bshare = {}
    for s in SITES:
        d = obs[obs.is_cpd & (obs.site == s)]
        b = d.groupby("eos").batch.apply(lambda x: (lambda v: ((v[:, None] == v[None, :]).sum() - len(v)) /
                                                     max(len(v) * (len(v) - 1), 1))(x.to_numpy()))
        bshare[s] = round(float(b.mean()), 3)
    pos["frac_within_site_replicate_pairs_same_batch"] = bshare
    ucols = sorted({w[1:] for w in obs.well[obs.is_dmso]})
    pos["dmso_columns"] = ucols
    pos["cpd_columns"] = f"{obs.well[obs.is_cpd].str[1:].astype(int).min()}-{obs.well[obs.is_cpd].str[1:].astype(int).max()}"
    # DMSO cosine: same well on different plates of a site vs different wells
    sa, sd = [], []
    for s in SITES:
        di = np.where((obs.is_dmso & (obs.site == s)).to_numpy())[0]
        G = U[di] @ U[di].T
        w = obs.well.to_numpy()[di]
        p = obs.plate.to_numpy()[di]
        iu = np.triu_indices(len(di), 1)
        diffp = p[iu[0]] != p[iu[1]]
        same = (w[iu[0]] == w[iu[1]]) & diffp
        sa.append(G[iu][same])
        sd.append(G[iu][diffp & ~same])
    pos["dmso_cos_same_well_diff_plate"] = round(float(np.concatenate(sa).mean()), 4)
    pos["dmso_cos_diff_well_diff_plate"] = round(float(np.concatenate(sd).mean()), 4)
    # compounds: cosine of DIFFERENT compounds at the same well (other layout, same site, same replicate
    # index) vs at different wells; and replicate pairs within site vs across sites
    cpd = obs[obs.is_cpd]
    key = cpd.set_index(["site", "rep", "well"])
    same_c, diff_c = [], []
    rs = np.random.default_rng(SEED + 1)
    ci = cpd.index.to_numpy()
    for s in SITES:
        for rp in ["R1", "R2", "R3", "R4"]:
            sub = cpd[(cpd.site == s) & (cpd.rep == rp)]
            ii = sub.index.to_numpy()
            G = U[ii] @ U[ii].T
            w = sub.well.to_numpy()
            lay = sub.layout.to_numpy()
            iu = np.triu_indices(len(ii), 1)
            difflay = lay[iu[0]] != lay[iu[1]]
            sw = (w[iu[0]] == w[iu[1]]) & difflay
            same_c.append(G[iu][sw])
            dw = np.where(difflay & ~sw)[0]
            diff_c.append(G[iu][rs.choice(dw, min(len(dw), 50000), replace=False)])
    pos["cpd_cos_diff_cpd_same_well"] = round(float(np.concatenate(same_c).mean()), 4)
    pos["cpd_cos_diff_cpd_diff_well"] = round(float(np.concatenate(diff_c).mean()), 4)
    run["position"] = pos
    log(f"position: {pos}")

    # --- consensus for ALL compounds (needed for active_z reference distribution)
    t0 = time.time()
    cons_all = {}
    for e, idx in rows_by(obs, obs.is_cpd, "eos").items():
        cons_all[e] = consensus(V[idx], site_arr[idx])
    norms_all = np.array([np.linalg.norm(v) for v in cons_all.values()])
    med_n = np.median(norms_all)
    mad_n = MAD_SCALE * np.median(np.abs(norms_all - med_n))

    # --- per unit
    cpd_rows_by_plate = {k: np.sort(v[obs.is_cpd.to_numpy()[v]]) for k, v in obs.groupby(["site", "plate"]).indices.items()}
    eos_arr = obs.eos.to_numpy()
    rows, cons_rows = [], []
    norm_null_cache = {}

    def norm_null(idx):
        key = tuple(sorted(plates_key[i] for i in idx))
        if key in norm_null_cache:
            return norm_null_cache[key]
        pk = [plates_key[i] for i in idx]
        sites_i = site_arr[idx]
        vals = np.empty(N_NULL)
        # draw one DMSO per replicate well from its plate (without replacement within a plate)
        by_p = {}
        for j, p in enumerate(pk):
            by_p.setdefault(p, []).append(j)
        for t in range(N_NULL):
            pick = np.empty(len(idx), int)
            for p, js in by_p.items():
                pick[js] = rng.choice(dbp[p], len(js), replace=False)
            vals[t] = np.linalg.norm(consensus(V[pick], sites_i))
        norm_null_cache[key] = vals
        return vals

    def vs_cpd(idx, own_eos):
        """mean pairwise cosine of replicates vs a null of other compounds' wells from the same plates
        (one per replicate well; distinct compounds within a draw)."""
        obsv = mean_pair_cos(U[idx])
        pk = [plates_key[i] for i in idx]
        null = np.empty(N_NULL)
        cand = {p: cpd_rows_by_plate[p][~np.isin(eos_arr[cpd_rows_by_plate[p]], own_eos)] for p in set(pk)}
        for t in range(N_NULL):
            used = set()
            pick = np.empty(len(idx), int)
            order = rng.permutation(len(idx))
            for j in order:
                c = cand[pk[j]]
                while True:
                    w = c[rng.integers(len(c))]
                    if eos_arr[w] not in used:
                        break
                used.add(eos_arr[w])
                pick[j] = w
            null[t] = mean_pair_cos(U[pick])
        return obsv, pval(null, obsv), float(null.mean())

    row_of = {k: i for i, k in enumerate(zip(obs.site, obs.layout, obs.rep, obs.well))}
    layouts = sorted(obs.layout.unique())
    lay_arr, rep_arr, well_arr = obs.layout.to_numpy(), obs.rep.to_numpy(), obs.well.to_numpy()

    def within_site_cos(ix, groups, gsites):
        """mean within-group pair cosine, overall and per site."""
        tot, cnt = {}, {}
        for gidx, gs in zip(groups, gsites):
            if len(gidx) >= 2:
                G = U[ix[gidx]] @ U[ix[gidx]].T
                tot[gs] = tot.get(gs, 0.0) + G.sum() - np.trace(G)
                cnt[gs] = cnt.get(gs, 0) + len(gidx) * (len(gidx) - 1)
        out = {k: tot[k] / cnt[k] for k in tot}
        out["all"] = sum(tot.values()) / sum(cnt.values()) if cnt else np.nan
        return out

    def vs_pos(idx):
        """Position-matched test. Score = mean cosine over WITHIN-SITE replicate pairs (same compound id,
        same well, different replicate plates). Null: each replicate well is replaced by the well at the
        SAME position, same site and same replicate index on a plate of another layout (i.e. a different
        compound), distinct layouts within a site, so null pairs share position but never compound."""
        groups, gl = [], {}
        for j, i in enumerate(idx):
            gl.setdefault((site_arr[i], eos_arr[i]), []).append(j)
        groups = [np.array(v) for v in gl.values()]
        gsites = [k[0] for k in gl]
        obsv = within_site_cos(np.asarray(idx), groups, gsites)
        null = []
        for t in range(N_NULL):
            pick = np.empty(len(idx), int)
            ok = True
            for (s, e), js in gl.items():
                own_lay = lay_arr[idx[js[0]]]
                perm = [l for l in rng.permutation(layouts) if l != own_lay]
                k = 0
                for j in js:
                    i = idx[j]
                    while k < len(perm) and (s, perm[k], rep_arr[i], well_arr[i]) not in row_of:
                        k += 1
                    if k == len(perm):
                        ok = False
                        break
                    pick[j] = row_of[(s, perm[k], rep_arr[i], well_arr[i])]
                    k += 1
            if ok:
                null.append(within_site_cos(pick, groups, gsites))
        res = {}
        for k, v in obsv.items():
            nv = np.array([d[k] for d in null])
            res[k] = (v, pval(nv, v), float(nv.mean()))
        return res, len(null)

    for un, idx in rows_by(obs, obs.role == "matched_compound", "unit").items():
        idx = np.sort(idx)
        own = sorted(set(eos_arr[idx]))
        pl = [plates_key[i] for i in idx]
        neg = np.concatenate([dbp[p] for p in sorted(set(pl))])
        m = rep_map(U, idx, neg)
        p = pval(mapnull.get(pl), m)
        row = dict(unit=un, dose_uM=DOSE_UM, cell_line=CELL_LINE, time_h=TIME_H, n_wells=len(idx),
                   rep_mAP=round(m, 4), rep_p=p, active_rep=p < P_ACTIVE)
        c = consensus(V[idx], site_arr[idx])
        nrm = float(np.linalg.norm(c))
        row["norm"] = round(nrm, 4)
        row["norm_z"] = round((nrm - med_n) / mad_n, 3)
        row["active_z"] = row["norm_z"] > Z_ACTIVE
        nn = norm_null(idx)
        row["norm_dmso_p"] = pval(nn, nrm)
        row["active_norm_dmso"] = row["norm_dmso_p"] < P_ACTIVE
        sc, pc, nm = vs_cpd(idx, own)
        row.update(rep_cos=round(sc, 4), rep_cos_null_cpd_mean=round(nm, 4), vs_cpd_p=pc, active_vs_cpd=pc < P_ACTIVE)
        vp, nn_ok = vs_pos(idx)
        sp, pp, npm = vp["all"]
        row.update(rep_cos_within_site=round(sp, 4), rep_cos_null_samepos_mean=round(npm, 4), vs_pos_p=pp,
                   active_vs_pos=pp < P_ACTIVE, n_null_pos=nn_ok)
        for s_ in SITES:
            row[f"vs_pos_p_{s_}"] = vp[s_][1] if s_ in vp else np.nan
        for s in SITES:
            si = idx[site_arr[idx] == s]
            if len(si) >= 2:
                ps = [plates_key[i] for i in si]
                ms = rep_map(U, si, np.concatenate([dbp[q] for q in sorted(set(ps))]))
                row[f"rep_mAP_{s}"] = round(ms, 4)
                row[f"rep_p_{s}"] = pval(mapnull.get(ps), ms)
            else:
                row[f"rep_mAP_{s}"], row[f"rep_p_{s}"] = np.nan, np.nan
        row["n_sites_reproducible"] = int(sum((row[f"rep_p_{s}"] < P_ACTIVE) for s in SITES))
        row["eos_ids"] = ";".join(own)
        row["wells"] = ";".join(sorted(set(obs.well.to_numpy()[idx])))
        rows.append(row)
        cons_rows.append(dict(unit=un, dose_uM=DOSE_UM, cell_line=CELL_LINE, time_h=TIME_H, n_wells=len(idx),
                              n_plates=len(set(pl)), plates=";".join(f"{s}:{q}" for s, q in sorted(set(pl))),
                              sites=";".join(sorted(set(site_arr[idx]))), vec=c))
    act = pd.DataFrame(rows).sort_values("unit").reset_index(drop=True)
    act.insert(act.columns.get_loc("active_vs_pos") + 1, "active_cpd_and_pos",
               act.active_vs_cpd & act.active_vs_pos)

    # --- cross-site agreement (matched units): replicate cosine within vs across sites, and agreement of
    # per-site consensus profiles (same unit vs different units)
    wi, cr, site_cons = [], [], {}
    for un, idx in rows_by(obs, obs.role == "matched_compound", "unit").items():
        G = U[idx] @ U[idx].T
        ss = site_arr[idx]
        iu = np.triu_indices(len(idx), 1)
        same = ss[iu[0]] == ss[iu[1]]
        wi.append(G[iu][same])
        cr.append(G[iu][~same])
        site_cons[un] = {k: unitize(np.median(V[idx[ss == k]], axis=0)[None])[0] for k in SITES if (ss == k).any()}
    same_u = [float(a @ b) for d in site_cons.values() for (k1, a), (k2, b) in combinations(sorted(d.items()), 2)]
    pair_cos = {f"{k1}-{k2}": round(float(np.mean([d[k1] @ d[k2] for d in site_cons.values() if k1 in d and k2 in d])), 4)
                for k1, k2 in combinations(SITES, 2)}
    names = sorted(site_cons)
    diff_u = []
    rs = np.random.default_rng(SEED + 2)
    for _ in range(5000):
        a, b = rs.choice(len(names), 2, replace=False)
        k1, k2 = rs.choice(SITES, 2, replace=False)
        if k1 in site_cons[names[a]] and k2 in site_cons[names[b]]:
            diff_u.append(float(site_cons[names[a]][k1] @ site_cons[names[b]][k2]))
    run["cross_site"] = {"rep_cos_within_site": round(float(np.concatenate(wi).mean()), 4),
                         "rep_cos_across_sites": round(float(np.concatenate(cr).mean()), 4),
                         "site_consensus_cos_same_unit": round(float(np.mean(same_u)), 4),
                         "site_consensus_cos_diff_unit_diff_site": round(float(np.mean(diff_u)), 4),
                         "site_consensus_cos_same_unit_by_pair": pair_cos}
    log(f"cross-site: {run['cross_site']}")
    log(f"unit metrics: {len(act)} units ({time.time() - t0:.0f}s)")

    # DMSO reference consensus
    di = np.where(obs.is_dmso.to_numpy())[0]
    cons_rows.append(dict(unit="DMSO", dose_uM=0.0, cell_line=CELL_LINE, time_h=TIME_H, n_wells=len(di),
                          n_plates=len(set(plates_key[i] for i in di)), plates="all",
                          sites=";".join(SITES), vec=consensus(V[di], site_arr[di])))

    # --- positive controls (excluded from the compound set, reported)
    crow = []
    for name in CONTROLS:
        for s in SITES + ["all"]:
            sel = (obs.eos == name) & (obs.conc == 5)
            if s != "all":
                sel &= obs.site == s
            idx = np.where(sel.to_numpy())[0]
            pl = [plates_key[i] for i in idx]
            neg = np.concatenate([dbp[q] for q in sorted(set(pl))])
            m = rep_map(U, idx, neg)
            r = dict(control=name, dose_uM=5.0, site=s, n_wells=len(idx), n_plates=len(set(pl)), rep_mAP=round(m, 4))
            if s != "all":
                null = MapNull(U, dbp, rng, N_NULL_CTRL).get(pl)
                r["rep_p"] = pval(null, m)
                r["null_mAP_mean"] = round(float(null.mean()), 4)
            c = consensus(V[idx], site_arr[idx])
            r["norm"] = round(float(np.linalg.norm(c)), 3)
            r["norm_z_vs_cpds"] = round((r["norm"] - med_n) / mad_n, 2)
            z0 = np.where(((obs.eos == name) & (obs.conc == 0) & (obs.site == s if s != "all" else True)).to_numpy())[0]
            r["n_wells_conc0"] = len(z0)
            if len(z0):
                r["conc0_cos_to_5uM_consensus"] = round(float((U[z0] @ unitize(c[None])[0]).mean()), 3)
                r["conc0_norm"] = round(float(np.linalg.norm(consensus(V[z0], site_arr[z0]))), 3)
            crow.append(r)
    ctrl = pd.DataFrame(crow)
    ctrl.to_csv(OUT_CONTROLS, index=False)
    log("controls:\n" + ctrl.to_string(index=False))

    # --- outputs
    k = V.shape[1]
    fcols = [f"f{i}" for i in range(k)]
    prof = pd.DataFrame([{kk: vv for kk, vv in r.items() if kk != "vec"} for r in cons_rows])
    prof = pd.concat([prof, pd.DataFrame(np.vstack([r["vec"] for r in cons_rows]).astype(np.float32),
                                         columns=fcols)], axis=1)
    prof.to_parquet(OUT_PROFILES, index=False)
    spec_cols = ["unit", "dose_uM", "cell_line", "time_h", "n_wells", "rep_mAP", "rep_p", "active_rep",
                 "norm_z", "active_z"]
    act = act[spec_cols + [c for c in act.columns if c not in spec_cols]]
    act.to_csv(OUT_ACTIVITY, index=False)
    wsel = obs.role.isin(["matched_compound", "negative_control"]).to_numpy()
    wells = obs.loc[wsel, ["site", "batch", "plate", "layout", "rep", "well", "eos", "unit", "role", "conc",
                           "object_count"]].reset_index(drop=True)
    wells = pd.concat([wells, pd.DataFrame(V[wsel], columns=fcols)], axis=1)
    wells.to_parquet(OUT_WELLS, index=False)

    # --- summaries
    def frac(col, d=act):
        return round(float(d[col].mean()), 4)

    summ = {"n_units": len(act), "reproducible_overall": frac("active_rep"),
            "reproducible_per_site": {s: round(float((act[f"rep_p_{s}"] < P_ACTIVE).mean()), 4) for s in SITES},
            "n_sites_reproducible_dist": act.n_sites_reproducible.value_counts().sort_index().to_dict(),
            "active_z": frac("active_z"), "active_norm_dmso": frac("active_norm_dmso"),
            "active_vs_cpd": frac("active_vs_cpd"), "active_vs_pos": frac("active_vs_pos"),
            "active_cpd_and_pos": frac("active_cpd_and_pos"),
            "active_vs_pos_per_site": {s: round(float((act[f"vs_pos_p_{s}"] < P_ACTIVE).mean()), 4) for s in SITES},
            "rep_mAP_median": round(float(act.rep_mAP.median()), 4),
            "rep_mAP_median_per_site": {s: round(float(act[f"rep_mAP_{s}"].median()), 4) for s in SITES},
            "n_wells_dist": act.n_wells.value_counts().sort_index().to_dict(),
            "norm_ref_median_all_cpds": round(float(med_n), 4), "norm_ref_mad_all_cpds": round(float(mad_n), 4),
            "n_cpds_ref": len(norms_all)}
    # compound wells null-mean of mAP for orientation
    summ["mapnull_mean_overall_example"] = round(float(np.mean([v.mean() for kk, v in mapnull.cache.items()
                                                                if len(kk) >= 12])), 4)
    run["summary"] = summ
    log(f"summary: {json.dumps(summ)}")

    track_disk()
    peak = STATE["peak_disk"]
    cache_plates = dir_bytes(PLATE_CACHE)
    if not KEEP_CACHE:
        shutil.rmtree(PLATE_CACHE, ignore_errors=True)
    t_total = time.time() - t_start
    total_bytes = Bytes.n
    dl = {"bytes": int(total_bytes), "peak_disk_bytes": int(peak), "seconds": round(t_total, 1),
          "zenodo_member_bytes_this_run": int(STATE["zenodo_bytes"]),
          "plate_cache_bytes": int(cache_plates), "plate_cache_deleted": not KEEP_CACHE,
          "note": "bytes = all HTTP bytes this run (gallery listings, platemaps, Zenodo API, zip directory, "
                  "plate members). A run from a warm plate cache downloads ~0 plate bytes."}
    if len(new) == len(members):
        OUT_SUMMARY.write_text(json.dumps(dl, indent=1), encoding="utf-8")
    run["download_summary"] = dl
    run["zenodo"] = {k2: v for k2, v in zen.items() if k2 != "all_versions"}
    run["zenodo_versions"] = zen["all_versions"]
    run["log"] = LOG
    OUT_RUN.write_text(json.dumps(run, indent=1, default=str), encoding="utf-8")
    write_notes(run, act, ctrl)
    log(f"done in {t_total:.0f}s; bytes {total_bytes / 1e9:.3f} GB; peak disk {peak / 1e9:.3f} GB")


if __name__ == "__main__":
    main()
