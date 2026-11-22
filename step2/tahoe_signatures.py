#!/usr/bin/env python
"""Tahoe-100M transcriptomic signatures for the ShapeShift morphology <-> transcriptome comparison.

Source: tahoebio/Tahoe-100M  metadata/pseudobulk_differential_expression/train-*-of-01026.parquet
(pseudobulk DE: one log2FoldChange per gene x cell line x drug x dose x plate, vs same-plate DMSO).
Everything is streamed with HTTP range requests; nothing large is written to disk.

Stages (`--stage all` = pilot -> fetch2 -> index -> build):
  pilot  read a few footers, fetch ~20 full combinations, verify layout + fixed gene order,
         pick the candidate gene set (tahoe_cache/pilot_genes.parquet, candidate_genes.txt)
  fetch  (v1, historical: 186 drugs x {0.5, 5} uM) read every needed footer -> cache a compact
         row-group index (files/NNNN_rg.parquet) + fragments (files/NNNN_frag.npz). Its outputs are
         reused as a cache by fetch2; it is not required for a from-scratch run.
  fetch2 ALL combinations (every drug, every dose; needed for the 324-drug x 3-dose scope and the
         all-drug activity baseline). Uses the cached row-group index (reads the footer only when
         the index is missing), fetches each needed log2FoldChange column chunk with one range
         request and decodes it with tahoe_pqdecode.py (validated against pyarrow) -> files/NNNN_frag2.npz
  index  consolidate the row-group index + combination labels (tahoe_cache/de_index_*.parquet), QC
  build  well-level combos -> per-line unit signatures -> consensus, consistency, null, activity, outputs
  perplate  per-plate (not plate-averaged) signatures for unit x dose on >=2 plates (noise ceiling input)

Layout facts this relies on (verified by the pilot and re-checked in build):
  * every combination has exactly 62,710 rows, genes in one fixed (Ensembl-ID) order
  * rows are contiguous across row groups (1000 rows) and across files; files 0-99 have 3,986,181
    rows and files 100-1025 have 3,986,180, so combo index = global_row // 62710
  * files are sorted by cell line, then plate (numeric), then drug; NOT by plate.
"""
import argparse, bisect, hashlib, io, json, os, random, re, struct, sys, threading, time
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))   # -I drops the script dir; needed for the decoder
from tahoe_pqdecode import decode_float_chunk  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
STEP2 = ROOT / "step2"
CACHE = STEP2 / "tahoe_cache"
FILES_DIR = CACHE / "files"
URL_T = ("https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/metadata/"
         "pseudobulk_differential_expression/train-{:05d}-of-01026.parquet")
N_FILES = 1026
ROWS_BASE, N_EXTRA = 3_986_180, 100          # files < 100 carry one extra row
G = 62_710                                   # genes (rows) per combination
DOSES = {"high": 5.0, "low": 0.5}
SKIP_PLATES = {"1", "4", "7", "10"}         # plates with only 0.05 uM (+DMSO)
COLS = dict(gene_name=0, baseMean=1, lfc=2, plate=7, n_trt=8, n_ctrl=9, cellosaur=10,
            depmap=11, drug=12, conc=13, cellname=15)
KEYS = ["plate", "cellosaur", "depmap", "drug", "conc", "cellname"]
N_CAND_GENES = 5000
MIN_CELLS_TRT = 20         # v1 only (superseded by MIN_CELLS = 100 in build)
MAX_NAN_FRAC = 0.20
PILOT_FILES = [3, 118, 233, 348, 463, 578, 693, 808, 923, 1025]
SEED = 0


def file_row_start(i):
    return i * ROWS_BASE + min(i, N_EXTRA)


def file_nrows(i):
    return ROWS_BASE + (1 if i < N_EXTRA else 0)


# --------------------------------------------------------------------------------------- HTTP
_tls = threading.local()


def _session():
    s = getattr(_tls, "s", None)
    if s is None:
        s = requests.Session()
        s.mount("https://", requests.adapters.HTTPAdapter(pool_connections=4, pool_maxsize=64))
        s.headers["User-Agent"] = "python-requests (range reader)"
        _tls.s = s
    return s


class Remote:
    """Range reader for one DE file; caches the CDN redirect, retries 429/5xx with backoff."""

    def __init__(self, i):
        self.i, self.url, self.size = i, None, None
        self.nbytes = self.nreq = self.retries = 0
        self.codes = {}
        self.lock = threading.Lock()

    def _get(self, rng):
        for k in range(12):
            url = self.url or URL_T.format(self.i)
            code, wait = None, None
            try:
                r = _session().get(url, headers={"Range": rng}, timeout=90)
                code = r.status_code
                if code == 206:
                    data = r.content
                    cr = r.headers.get("Content-Range", "")
                    a, b = cr.split(" ")[1].split("/")[0].split("-")
                    if len(data) == int(b) - int(a) + 1:
                        with self.lock:
                            if self.url is None:
                                self.url = r.url
                            if self.size is None:
                                self.size = int(cr.split("/")[1])
                            self.nbytes += len(data)
                            self.nreq += 1
                        return data
                elif code in (401, 403, 404, 410) and self.url is not None:
                    self.url = None                       # signed URL expired -> re-resolve
                elif code not in (408, 425, 429, 500, 502, 503, 504):
                    raise RuntimeError(f"file {self.i} {rng}: HTTP {code}")
                ra = r.headers.get("Retry-After")
                if ra and ra.strip().isdigit():
                    wait = min(120.0, float(ra))
            except requests.RequestException as e:
                code = type(e).__name__
            with self.lock:
                self.retries += 1
                self.codes[str(code)] = self.codes.get(str(code), 0) + 1
            time.sleep(wait if wait is not None else min(60.0, 0.5 * 2 ** k) + random.random())
        raise RuntimeError(f"file {self.i} {rng}: giving up (last HTTP {code})")

    def get(self, start, end_excl):
        end_excl = min(end_excl, self.size) if self.size else end_excl
        return self._get(f"bytes={start}-{end_excl - 1}")

    def tail(self, n):
        return self._get(f"bytes=-{n}")


class Sparse(io.RawIOBase):
    """Read-only file object over a set of fetched byte ranges (enough for pyarrow)."""

    def __init__(self, size):
        self.size, self.pos, self.starts, self.chunks = size, 0, [], []

    def add(self, start, data):
        k = bisect.bisect_left(self.starts, start)
        self.starts.insert(k, start)
        self.chunks.insert(k, data)

    def readable(self):
        return True

    def seekable(self):
        return True

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else self.size + off)
        return self.pos

    def tell(self):
        return self.pos

    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        k = bisect.bisect_right(self.starts, self.pos) - 1
        # several chunks may start at or before pos; search backwards for one that covers it
        while k >= 0:
            st, data = self.starts[k], self.chunks[k]
            if self.pos + n <= st + len(data):
                b[:n] = data[self.pos - st:self.pos - st + n]
                self.pos += n
                return n
            k -= 1
        raise IOError(f"range not fetched: {self.pos}+{n}")


def chunk_range(c):
    st = c.dictionary_page_offset if c.has_dictionary_page else c.data_page_offset
    return st, st + c.total_compressed_size


def open_footer(i):
    """Fetch + parse the footer of file i. Returns (remote, sparse, ParquetFile, rg DataFrame)."""
    rem = Remote(i)
    tail = rem.tail(8)
    assert tail[4:] == b"PAR1", tail
    flen = struct.unpack("<i", tail[:4])[0]
    fstart = rem.size - 8 - flen
    sp = Sparse(rem.size)
    sp.add(fstart, rem.get(fstart, rem.size))
    pf = pq.ParquetFile(sp, pre_buffer=False)
    md = pf.metadata
    assert md.num_rows == file_nrows(i), (i, md.num_rows)
    names = md.schema.names
    assert names[COLS["lfc"]] == "log2FoldChange" and names[COLS["drug"]] == "drug", names
    rec = {k: [] for k in ["rg", "num_rows", "row_start", "c0_off", "c2_end", "lfc_off", "lfc_size",
                           "rg_end"] + [f"{k}_{m}" for k in KEYS + ["n_trt", "n_ctrl"] for m in ("min", "max")]}
    row = file_row_start(i)
    for g in range(md.num_row_groups):
        rg = md.row_group(g)
        rec["rg"].append(g)
        rec["num_rows"].append(rg.num_rows)
        rec["row_start"].append(row)
        row += rg.num_rows
        c0, c2 = rg.column(0), rg.column(COLS["lfc"])
        rec["c0_off"].append(chunk_range(c0)[0])
        a, b = chunk_range(c2)
        rec["lfc_off"].append(a)
        rec["lfc_size"].append(b - a)
        rec["c2_end"].append(b)
        rec["rg_end"].append(chunk_range(rg.column(rg.num_columns - 1))[1])
        for k in KEYS + ["n_trt", "n_ctrl"]:
            st = rg.column(COLS[k]).statistics
            rec[f"{k}_min"].append(st.min if st is not None and st.has_min_max else None)
            rec[f"{k}_max"].append(st.max if st is not None and st.has_min_max else None)
    df = pd.DataFrame(rec)
    df.insert(0, "file", i)
    return rem, sp, pf, df


def label_combos(df):
    """Assign combo ids to row groups by arithmetic; verify against stats; labels from interior rgs."""
    c_first = df.row_start // G
    c_last = (df.row_start + df.num_rows - 1) // G
    df["combo_first"], df["combo_last"] = c_first, c_last
    same = np.ones(len(df), bool)
    for k in KEYS + ["n_trt", "n_ctrl"]:
        same &= (df[f"{k}_min"].values == df[f"{k}_max"].values)
    single = (c_first == c_last).values
    bad = int((single != same).sum())       # interior rg must have min==max; spanning rg must not
    inner = df[single]
    lab = inner.groupby("combo_first").agg(
        **{k: (f"{k}_min", "first") for k in KEYS + ["n_trt", "n_ctrl"]},
        n_rg=("rg", "size"),
        **{f"nu_{k}": (f"{k}_min", "nunique") for k in KEYS + ["n_trt", "n_ctrl"]})
    nu = lab[[c for c in lab.columns if c.startswith("nu_")]]
    bad += int((nu > 1).any(axis=1).sum())
    lab = lab.drop(columns=nu.columns).reset_index().rename(columns={"combo_first": "combo"})
    return lab, bad


def load_targets():
    comp = pd.read_csv(STEP2 / "compounds.csv")
    names = sorted({n.strip() for x in comp.tahoe_drugs for n in str(x).split(";") if n.strip()})
    return comp, names


def is_needed(drug, conc, names_set):
    if drug is None or conc is None:
        return False
    return drug.strip() in names_set and any(abs(conc - v) <= 1e-3 * v for v in DOSES.values())


# --------------------------------------------------------------------------------------- pilot
def pilot(args):
    _, names = load_targets()
    names_set = set(names)
    rng = random.Random(SEED)
    t0 = time.time()
    total_bytes = 0
    ref_genes = None
    recs, bms, lfcs, meta = [], [], [], []
    n_boundary_checked = 0
    for i in PILOT_FILES:
        rem, sp, pf, df = open_footer(i)
        lab, bad = label_combos(df)
        assert bad == 0, f"pilot file {i}: {bad} stat/arithmetic mismatches"
        fs, fe = file_row_start(i), file_row_start(i) + file_nrows(i)
        lab["inside"] = (lab.combo * G >= fs) & ((lab.combo + 1) * G <= fe)
        cand = lab[lab.inside.values & np.array([is_needed(d, c, names_set) for d, c in zip(lab.drug, lab.conc)], bool)]
        picks = []
        for dose_v in DOSES.values():
            sub = cand[np.isclose(cand.conc, dose_v, rtol=1e-3)]
            if len(sub):
                picks.append(sub.iloc[rng.randrange(len(sub))])
        if len(cand) > len(picks):
            rest = cand[~cand.combo.isin([p.combo for p in picks])]
            picks.append(rest.iloc[rng.randrange(len(rest))])
        for p in picks:
            rgs = df[(df.combo_last >= p.combo) & (df.combo_first <= p.combo)]
            a, b = int(rgs.c0_off.min()), int(rgs.rg_end.max())
            sp.add(a, rem.get(a, b + 128))
            t = pf.read_row_groups(rgs.rg.tolist(), columns=["gene_name", "baseMean", "log2FoldChange",
                                                               "plate", "drug", "concentration", "Cell_Name_Vevo"]).to_pandas()
            grow = np.concatenate([np.arange(s, s + n) for s, n in zip(rgs.row_start, rgs.num_rows)])
            in_combo = (grow // G) == p.combo
            # exact boundary check: labels on every row agree with the arithmetic assignment
            match = ((t.drug == p.drug) & (t.plate == p.plate) & (t.Cell_Name_Vevo == p.cellname)
                     & np.isclose(t.concentration, p.conc)).values
            assert (match[in_combo]).all() and not match[~in_combo].any(), f"boundary mismatch file {i}"
            n_boundary_checked += 1
            sub = t[in_combo].reset_index(drop=True)
            assert len(sub) == G
            if ref_genes is None:
                ref_genes = sub.gene_name.values.copy()
            assert (sub.gene_name.values == ref_genes).all(), f"gene order differs in file {i}"
            bms.append(sub.baseMean.values.astype(np.float64))
            lfcs.append(sub.log2FoldChange.values.astype(np.float64))
            meta.append(dict(file=i, combo=int(p.combo), drug=p.drug, conc=float(p.conc), plate=p.plate,
                             cell_line=p.cellname, n_trt=int(p.n_trt), n_ctrl=int(p.n_ctrl)))
        total_bytes += rem.nbytes
        print(f"pilot file {i}: {len(picks)} combos, {rem.nbytes / 1e6:.1f} MB", flush=True)
    bm, lf = np.vstack(bms), np.vstack(lfcs)
    pg = pd.DataFrame({"pos": np.arange(G), "gene_name": ref_genes,
                       "median_baseMean": np.median(bm, 0), "min_baseMean": bm.min(0),
                       "frac_lfc_ok": np.isfinite(lf).mean(0)})
    excl = pg.gene_name.str.match(r"^(ENSG\d|MT-|Y_RNA|hsa-mir|snoZ|7SK|5S_rRNA|5_8S_rRNA)") | \
        pg.gene_name.str.contains("rRNA")
    pg["eligible"] = (~excl) & (pg.frac_lfc_ok >= 0.9) & pg.gene_name.notna()
    pg["rank"] = pg.median_baseMean.where(pg.eligible).rank(ascending=False, method="first")
    pg["candidate"] = pg["rank"] <= N_CAND_GENES
    CACHE.mkdir(parents=True, exist_ok=True)
    pg.to_parquet(CACHE / "pilot_genes.parquet", index=False)
    pd.DataFrame(meta).to_csv(CACHE / "pilot_combos.csv", index=False)
    cand = pg[pg.candidate].sort_values("pos")
    (CACHE / "candidate_genes.txt").write_text("\n".join(f"{p}\t{g}" for p, g in zip(cand.pos, cand.gene_name)) + "\n")
    info = dict(n_pilot_combos=len(meta), n_pilot_files=len(PILOT_FILES), bytes=total_bytes,
                seconds=round(time.time() - t0, 1), n_boundary_checks=n_boundary_checked,
                gene_order_fixed=True, n_candidates=int(pg.candidate.sum()),
                candidate_min_median_baseMean=float(cand.median_baseMean.min()),
                n_excluded_by_name=int(excl.sum()))
    (CACHE / "pilot_info.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))


def load_candidates():
    lines = (CACHE / "candidate_genes.txt").read_text().splitlines()
    pos = np.array([int(l.split("\t")[0]) for l in lines])
    names = [l.split("\t")[1] for l in lines]
    return pos, names


def gene_hash(pos):
    return hashlib.md5(np.asarray(pos, np.int64).tobytes()).hexdigest()[:12]


# --------------------------------------------------------------------------------------- fetch
def process_file(i, threads):
    t0 = time.time()
    pos, _ = load_candidates()
    gh = gene_hash(pos)
    out_rg, out_fr = FILES_DIR / f"{i:04d}_rg.parquet", FILES_DIR / f"{i:04d}_frag.npz"
    if out_rg.exists() and out_fr.exists():
        z = np.load(out_fr)
        if str(z["gene_hash"]) == gh:
            return dict(file=i, status="cached")
    _, names = load_targets()
    names_set = set(names)
    rem, sp, pf, df = open_footer(i)
    footer_bytes = rem.nbytes
    lab, bad = label_combos(df)
    fs, fe = file_row_start(i), file_row_start(i) + file_nrows(i)
    all_combos = np.arange(fs // G, (fe - 1) // G + 1)
    labelled = dict(zip(lab.combo, zip(lab.drug, lab.conc)))
    # needed = labelled & needed, or not labelled in this file (edge fragment; decided at build)
    need = [c for c in all_combos if (c not in labelled) or is_needed(*labelled[c], names_set)]
    need = np.array(need, np.int64)
    rgs = df[[bool(np.isin(np.arange(a, b + 1), need).any()) for a, b in zip(df.combo_first, df.combo_last)]]
    data_bytes0 = rem.nbytes

    def fetch(row):
        a, n = int(row[0]), int(row[1])
        return a, rem.get(a, a + n + 128)

    if len(rgs):
        with ThreadPoolExecutor(threads) as ex:
            for a, d in ex.map(fetch, zip(rgs.lfc_off, rgs.lfc_size)):
                sp.add(a, d)
        vals = pf.read_row_groups(rgs.rg.tolist(), columns=["log2FoldChange"]).column(0).to_numpy(zero_copy_only=False)
        grow = np.concatenate([np.arange(s, s + n) for s, n in zip(rgs.row_start, rgs.num_rows)])
        assert len(vals) == len(grow)
    else:
        vals, grow = np.zeros(0, np.float32), np.zeros(0, np.int64)
    rank = np.full(G, -1, np.int64)
    rank[pos] = np.arange(len(pos))
    cidx = {c: k for k, c in enumerate(need)}
    frag = np.full((len(need), len(pos)), np.nan, np.float32)
    cov = np.zeros((len(need), len(pos)), bool)
    combo, p = grow // G, grow % G
    m = np.isin(combo, need) & (rank[p] >= 0)
    ci = np.array([cidx[c] for c in combo[m]], np.int64) if m.any() else np.zeros(0, np.int64)
    frag[ci, rank[p[m]]] = vals[m]
    cov[ci, rank[p[m]]] = True
    FILES_DIR.mkdir(parents=True, exist_ok=True)
    keep = ["file", "rg", "num_rows", "row_start", "lfc_off", "lfc_size", "combo_first", "combo_last"] + \
        [f"{k}_{mm}" for k in KEYS + ["n_trt", "n_ctrl"] for mm in ("min", "max")]
    tmp = out_rg.with_suffix(".tmp")
    df[keep].to_parquet(tmp, index=False)
    np.savez(FILES_DIR / f"{i:04d}_frag.tmp.npz", combos=need, frag=frag, cov=cov, gene_hash=gh)
    os.replace(FILES_DIR / f"{i:04d}_frag.tmp.npz", out_fr)
    os.replace(tmp, out_rg)
    return dict(file=i, status="done", stat_mismatch=bad, n_needed_combos=len(need), n_rg_fetched=len(rgs),
                footer_bytes=footer_bytes, data_bytes=rem.nbytes - data_bytes0, n_requests=rem.nreq,
                retries=rem.retries, seconds=round(time.time() - t0, 2))


def edge_labels(i):
    """(first rg labels, last rg labels) of a fetched file from its cached rg index."""
    df = pd.read_parquet(FILES_DIR / f"{i:04d}_rg.parquet", columns=["plate_min", "plate_max",
                                                                     "cellname_min", "cellname_max"])
    f, l = df.iloc[0], df.iloc[-1]
    return (f.plate_min, f.plate_max, f.cellname_min, f.cellname_max), (l.plate_min, l.plate_max, l.cellname_min, l.cellname_max)


def fetch(args):
    log_path = CACHE / "fetch_log.csv"
    old = pd.read_csv(log_path) if log_path.exists() else pd.DataFrame()
    logs = [old] if len(old) else []
    t0 = time.time()

    def run(files, phase):
        res = []
        with ProcessPoolExecutor(args.workers) as ex:
            futs = {ex.submit(process_file, i, args.threads): i for i in files}
            for k, f in enumerate(as_completed(futs)):
                r = f.result()
                r["phase"] = phase
                res.append(r)
                if r["status"] == "done" and (k % 20 == 0 or k == len(files) - 1):
                    done = pd.DataFrame([x for x in res if x["status"] == "done"])
                    print(f"[{phase}] {k + 1}/{len(files)} files, {time.time() - t0:.0f}s, "
                          f"{(done.footer_bytes.sum() + done.data_bytes.sum()) / 1e9:.2f} GB this run, "
                          f"retries {done.retries.sum()}", flush=True)
        return res

    even = [i for i in range(N_FILES) if i % 2 == 0 or i == N_FILES - 1]
    res = run(even, "even")
    odd, skipped = [], []
    for i in range(1, N_FILES - 1, 2):
        if (FILES_DIR / f"{i:04d}_frag.npz").exists() and (FILES_DIR / f"{i:04d}_rg.parquet").exists():
            odd.append(i)
            continue
        _, last_prev = edge_labels(i - 1)
        first_next, _ = edge_labels(i + 1)
        # sorted by (cell line, plate, drug): if the previous file ends and the next file starts in
        # the same (line, 0.05-only plate) block, file i lies wholly inside it -> nothing needed
        if (last_prev[0] == last_prev[1] == first_next[0] == first_next[1] and last_prev[0] in SKIP_PLATES
                and last_prev[2] == last_prev[3] == first_next[2] == first_next[3]):
            skipped.append(dict(file=i, status="skipped", phase="odd", block=f"{last_prev[2]}|plate{last_prev[0]}"))
        else:
            odd.append(i)
    res += run(odd, "odd") + skipped
    new = pd.DataFrame(res)
    new["run_t0"] = t0
    new["run_wall_s"] = round(time.time() - t0, 1)
    logs.append(new)
    log = pd.concat(logs, ignore_index=True)
    log.to_csv(log_path, index=False)
    done = new[new.status == "done"]
    print(f"fetch: {len(done)} files fetched, {len(skipped)} skipped, "
          f"{(done.footer_bytes.sum() + done.data_bytes.sum()) / 1e9:.2f} GB, {time.time() - t0:.0f}s")


# --------------------------------------------------------------------------------------- fetch2
DOSES3 = (0.05, 0.5, 5.0)
MIN_CELLS = 100            # per-condition treated-cell threshold (same as step2/gate_t2.py)
MIN_DRUGS_REF = 5          # min distinct drugs on a plate x dose (x line) for null / robust z
LINE_A549, LINE_HEPG2 = "CVCL_0023", "CVCL_1098"
RG_KEEP = ["file", "rg", "num_rows", "row_start", "lfc_off", "lfc_size", "combo_first", "combo_last"] + \
    [f"{k}_{mm}" for k in KEYS + ["n_trt", "n_ctrl"] for mm in ("min", "max")]


def process_file2(i, threads):
    """Fetch + decode log2FoldChange for every combination of file i not already in the v1 fragment."""
    t0 = time.time()
    pos, _ = load_candidates()
    gh = gene_hash(pos)
    out = FILES_DIR / f"{i:04d}_frag2.npz"
    if out.exists():
        try:
            if str(np.load(out)["gene_hash"]) == gh:
                return dict(file=i, status="cached")
        except Exception:
            pass
    FILES_DIR.mkdir(parents=True, exist_ok=True)
    rg_path = FILES_DIR / f"{i:04d}_rg.parquet"
    if rg_path.exists():
        df = pd.read_parquet(rg_path)
        rem = Remote(i)
        rem.tail(8)
        footer_bytes = rem.nbytes
    else:
        rem, _, _, df = open_footer(i)
        label_combos(df)
        footer_bytes = rem.nbytes
        df[RG_KEEP].to_parquet(rg_path.with_suffix(".tmp"), index=False)
        os.replace(rg_path.with_suffix(".tmp"), rg_path)
    oldz, old = None, np.zeros(0, np.int64)
    if (FILES_DIR / f"{i:04d}_frag.npz").exists():
        oldz = np.load(FILES_DIR / f"{i:04d}_frag.npz")
        if str(oldz["gene_hash"]) == gh:
            old = oldz["combos"].astype(np.int64)
            oldfrag = oldz["frag"]
        else:
            oldz = None
    fs, fe = file_row_start(i), file_row_start(i) + file_nrows(i)
    allc = np.arange(fs // G, (fe - 1) // G + 1)
    new = allc[~np.isin(allc, old)]
    rgs = df[np.isin(df.combo_first, new) | np.isin(df.combo_last, new)]
    rank = np.full(G, -1, np.int64)
    rank[pos] = np.arange(len(pos))
    frag = np.full((len(new), len(pos)), np.nan, np.float32)
    cov = np.zeros((len(new), len(pos)), bool)
    chk = dict(n=0, bad=0)
    b0 = rem.nbytes

    def work(r):
        a = int(r.lfc_off)
        return int(r.row_start), decode_float_chunk(rem.get(a, a + int(r.lfc_size) + 128), int(r.num_rows))

    with ThreadPoolExecutor(threads) as ex:
        for rs, v in ex.map(work, rgs.itertuples(index=False)):
            g = np.arange(rs, rs + len(v))
            c, rk = g // G, rank[g % G]
            m = np.isin(c, new) & (rk >= 0)
            if m.any():
                ci = np.searchsorted(new, c[m])
                frag[ci, rk[m]] = v[m]
                cov[ci, rk[m]] = True
            if oldz is not None:                      # free decoder cross-check vs pyarrow (v1) values
                mo = np.isin(c, old) & (rk >= 0)
                if mo.any():
                    a, e = v[mo], oldfrag[np.searchsorted(old, c[mo]), rk[mo]]
                    chk["n"] += int(mo.sum())
                    chk["bad"] += int(((np.isnan(a) != np.isnan(e)) | (~np.isnan(a) & (a != e))).sum())
    tmp = FILES_DIR / f"{i:04d}_frag2.tmp.npz"
    np.savez_compressed(tmp, combos=new, frag=frag, cov=cov, gene_hash=gh)
    os.replace(tmp, out)
    return dict(file=i, status="done", n_new_combos=len(new), n_rg_fetched=len(rgs), footer_bytes=footer_bytes,
                data_bytes=rem.nbytes - b0, n_requests=rem.nreq, retries=rem.retries, codes=json.dumps(rem.codes),
                decode_checked=chk["n"], decode_mismatch=chk["bad"], seconds=round(time.time() - t0, 2))


def fetch2(args):
    t0 = time.time()
    res = []
    with ProcessPoolExecutor(args.workers) as ex:
        futs = {ex.submit(process_file2, i, args.threads): i for i in range(N_FILES)}
        for k, f in enumerate(as_completed(futs)):
            try:
                r = f.result()
            except Exception as e:                    # never let one file kill the run; re-run resumes
                r = dict(file=futs[f], status="failed", error=repr(e)[:300])
            r["phase"] = "fetch2"
            res.append(r)
            if k % 25 == 0 or k == N_FILES - 1:
                d = pd.DataFrame(res)
                dd = d[d.status == "done"]
                print(f"[fetch2] {k + 1}/{N_FILES} {time.time() - t0:.0f}s "
                      f"{(dd.get('footer_bytes', 0).sum() + dd.get('data_bytes', 0).sum()) / 1e9 if len(dd) else 0:.2f} GB "
                      f"retries {dd.retries.sum() if len(dd) else 0} failed {(d.status == 'failed').sum()} "
                      f"decode_mismatch {dd.decode_mismatch.sum() if len(dd) else 0}", flush=True)
    new = pd.DataFrame(res)
    new["run_t0"], new["run_wall_s"] = t0, round(time.time() - t0, 1)
    lp = CACHE / "fetch_log.csv"
    log = pd.concat([pd.read_csv(lp), new], ignore_index=True) if lp.exists() else new
    log.to_csv(lp, index=False)
    disk = sum(p.stat().st_size for p in CACHE.rglob("*") if p.is_file())
    print(f"fetch2 done: {(new.status == 'done').sum()} done, {(new.status == 'cached').sum()} cached, "
          f"{(new.status == 'failed').sum()} failed, {time.time() - t0:.0f}s; cache {disk / 1e9:.2f} GB")
    if (new.status == "failed").any():
        print(new[new.status == "failed"][["file", "error"]].to_string())


# --------------------------------------------------------------------------------------- index
def dose_of(conc):
    for d in DOSES3:
        if abs(conc - d) <= 1e-3 * d:
            return d
    return np.nan


def index(args):
    rg_files = sorted(FILES_DIR.glob("*_rg.parquet"))
    have = {int(p.name[:4]) for p in rg_files}
    miss = sorted(set(range(N_FILES)) - have)
    assert not miss, f"row-group index missing for {len(miss)} files (run fetch2): {miss[:10]}"
    rgi = pd.concat([pd.read_parquet(p) for p in rg_files], ignore_index=True)
    n_bad, labs = 0, []
    for f, df in rgi.groupby("file", sort=True):
        assert df.num_rows.sum() == file_nrows(f) and df.row_start.iloc[0] == file_row_start(f), f
        lab, bad = label_combos(df.copy())
        n_bad += bad
        lab["file"] = f
        labs.append(lab)
    labs = pd.concat(labs, ignore_index=True)
    kk = KEYS + ["n_trt", "n_ctrl"]
    agg = labs.groupby("combo").agg(**{k: (k, "first") for k in kk}, n_rg_interior=("n_rg", "sum"),
                                    files=("file", lambda x: ";".join(map(str, sorted(set(x))))),
                                    **{f"nu_{k}": (k, "nunique") for k in kk})
    nu = [c for c in agg.columns if c.startswith("nu_")]
    n_conflict = int((agg[nu] > 1).any(axis=1).sum())
    combos = agg.drop(columns=nu).reset_index()
    combos["drug_stripped"] = combos.drug.str.strip()
    combos["dose_uM"] = combos.conc.map(dose_of)
    n_total = (file_row_start(N_FILES - 1) + file_nrows(N_FILES - 1)) // G
    cs = combos.sort_values("combo")
    blk = cs.cellname + "|" + cs.plate
    runs, n_blocks = int((blk != blk.shift()).sum()), int(blk.nunique())
    ctrl_nu = combos.groupby(["cellname", "plate"]).n_ctrl.nunique()
    dup = combos[combos.duplicated(["drug", "conc", "plate", "cellname"], keep=False)]
    dup_strip = combos[combos.duplicated(["drug_stripped", "dose_uM", "plate", "cellname"], keep=False)]
    rgi[["file", "rg", "row_start", "num_rows", "lfc_off", "lfc_size"]].astype(
        {"file": "int16", "rg": "int16", "num_rows": "int16", "lfc_size": "int32"}).to_parquet(
        CACHE / "de_index_rowgroups.parquet", index=False)
    combos.to_parquet(CACHE / "de_index_combos.parquet", index=False)
    summ = dict(n_rowgroups=len(rgi), n_files=len(have), n_combos_total=int(n_total), n_combos_labelled=len(combos),
                stat_arith_mismatches=n_bad, label_conflicts=n_conflict, n_line_plate_blocks=n_blocks,
                n_block_runs=runs, n_ctrl_nonconstant_blocks=int((ctrl_nu > 1).sum()),
                dup_key_combos=len(dup), dup_key_combos_after_strip=len(dup_strip),
                combos_per_dose=combos.dose_uM.value_counts(dropna=False).to_dict(),
                n_drugs=int(combos.drug_stripped.nunique()), n_lines=int(combos.cellosaur.nunique()),
                n_plates=int(combos.plate.nunique()))
    (CACHE / "index_summary.json").write_text(json.dumps(summ, indent=1, default=str))
    print(json.dumps(summ, indent=1, default=str))
    if len(dup_strip):
        print(dup_strip.sort_values(["drug_stripped", "plate", "cellname"]).head(20).to_string())
    return combos


# --------------------------------------------------------------------------------------- build
def _unit_rows(X):
    Z = X - X.mean(1, keepdims=True)
    n = np.linalg.norm(Z, axis=1, keepdims=True)
    n[n == 0] = 1
    return Z / n


def corr_matrix(X):
    Z = _unit_rows(np.asarray(X, np.float64))
    return Z @ Z.T


def mean_offdiag(R):
    k = len(R)
    return float((R.sum() - np.trace(R)) / (k * (k - 1))) if k > 1 else np.nan


def boot_ci(R, rng, B=1000):
    k = len(R)
    if k < 3:
        return np.nan, np.nan
    idx = rng.integers(k, size=(B, k))
    sub = R[idx[:, :, None], idx[:, None, :]]
    m = idx[:, :, None] != idx[:, None, :]
    v = (sub * m).sum((1, 2)) / m.sum((1, 2))
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def robust_z(x, ref):
    ref = np.asarray(ref, np.float64)
    med = np.median(ref)
    mad = 1.4826 * np.median(np.abs(ref - med))
    return (np.asarray(x, np.float64) - med) / mad if mad > 0 else np.full(np.shape(x), np.nan)


def load_units():
    u = pd.read_csv(STEP2 / "units.csv")
    rows = []
    for unit, g in u.groupby("unit"):
        sets = [frozenset(n.strip() for n in str(x).split(";") if n.strip()) for x in g.tahoe_drugs]
        rows.append(dict(unit=unit, drugs=sorted(set().union(*sets)), consistent=len(set(sets)) == 1,
                         tiers=";".join(sorted(set(g.tier))), moa=";".join(sorted({str(m) for m in g.moa if pd.notna(m)}))))
    U = pd.DataFrame(rows)
    m = pd.DataFrame([(d, r.unit) for r in U.itertuples() for d in r.drugs], columns=["drug_stripped", "unit"])
    names = sorted({l.strip() for l in (STEP2 / "tahoe_drugs_needed.txt").read_text(encoding="utf8").splitlines() if l.strip()})
    return U, m, names


def write_pq(df, path):
    import pyarrow as pa
    t = pa.Table.from_pandas(df, preserve_index=False)
    fcols = [c for c in df.columns if df[c].dtype == np.float32]
    pq.write_table(t, path, compression="zstd", use_dictionary=[c for c in df.columns if c not in fcols],
                   use_byte_stream_split=fcols)
    return path.stat().st_size


def build(args):
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    combos = pd.read_parquet(CACHE / "de_index_combos.parquet").sort_values("combo").reset_index(drop=True)
    pos, gnames = load_candidates()
    gh = gene_hash(pos)
    U, name2unit, names = load_units()
    assert not name2unit.drug_stripped.duplicated().any(), "a Tahoe name maps to >1 unit"
    assert set(names) == set(name2unit.drug_stripped), "tahoe_drugs_needed.txt != units.csv names"
    # ---- combo x gene matrix (all labelled combos)
    row_of = pd.Series(np.arange(len(combos)), index=combos.combo.values)
    X = np.full((len(combos), len(pos)), np.nan, np.float32)
    C = np.zeros((len(combos), len(pos)), np.int8)
    for p in sorted(FILES_DIR.glob("*_frag*.npz")):
        if ".tmp" in p.name:
            continue
        z = np.load(p)
        if str(z["gene_hash"]) != gh:
            continue
        cc = z["combos"].astype(np.int64)
        ok = np.isin(cc, row_of.index.values)
        fr, cv = z["frag"], z["cov"]
        for k in np.where(ok)[0]:
            r = row_of[cc[k]]
            X[r, cv[k]] = fr[k, cv[k]]
            C[r] += cv[k]
    cov_ok = (C == 1).all(1)
    del C
    assert cov_ok.all(), f"{(~cov_ok).sum()} combos not covered exactly once"
    combos["target"] = combos.drug_stripped.isin(names)
    combos = combos.merge(name2unit, on="drug_stripped", how="left")
    combos["ident"] = combos.unit.fillna("drug:" + combos.drug_stripped)      # unit for targets, name otherwise
    # ---- threshold
    combos["pass"] = combos.n_trt >= MIN_CELLS
    tg = combos[combos.target]
    line_tab = combos.groupby(["cellosaur", "cellname"]).agg(median_n_trt=("n_trt", "median"),
                                                             frac_dropped=("pass", lambda x: 1 - x.mean()),
                                                             n_conditions=("pass", "size")).reset_index()
    line_tab_t = tg.groupby("cellosaur").agg(target_conditions=("pass", "size"), target_dropped=("pass", lambda x: (~x).sum()))
    line_tab = line_tab.merge(line_tab_t, on="cellosaur", how="left").sort_values("median_n_trt")
    line_tab.to_csv(CACHE / "lines_cell_counts.csv", index=False)
    thr = dict(threshold=MIN_CELLS, conditions_all=len(combos), dropped_all=int((~combos["pass"]).sum()),
               conditions_target=len(tg), dropped_target=int((~tg["pass"]).sum()),
               dropped_target_by_dose={str(d): int((~g["pass"]).sum()) for d, g in tg.groupby("dose_uM")},
               lines_all_target_dropped=line_tab[line_tab.target_dropped == line_tab.target_conditions].cellosaur.tolist(),
               lines_gt50pct_dropped=line_tab[line_tab.frac_dropped > 0.5][["cellosaur", "cellname", "median_n_trt", "frac_dropped"]].to_dict("records"),
               n_lines_with_any_target_kept=int(tg[tg["pass"]].cellosaur.nunique()),
               n_trt_quantiles_target={str(q): float(tg.n_trt.quantile(q)) for q in (0, .01, .05, .1, .25, .5, .75, .9, 1)},
               n_lt20_target=int((tg.n_trt < 20).sum()))
    keep = combos["pass"].values
    nanfrac = np.isnan(X[keep]).mean(0)
    gk = nanfrac <= MAX_NAN_FRAC
    genes = [g for g, k in zip(gnames, gk) if k]
    K = combos[keep].reset_index(drop=True)
    XK = np.nan_to_num(X[keep][:, gk], nan=0.0)
    del X
    print(f"kept {len(K)}/{len(combos)} combos (n_trt>={MIN_CELLS}); genes {len(genes)}/{len(gnames)}", flush=True)
    norms = np.linalg.norm(XK.astype(np.float64), axis=1)
    K["norm"] = norms
    # plate x dose x line mean over ALL drugs on the plate: the component every drug shares through the
    # common DMSO control (+ generic treated-vs-DMSO effects). Used for diagnostics, not removed from outputs.
    XKc = XK.copy()
    pm_rows, pm_X = [], []
    for (plate, dose, cl), g in K.groupby(["plate", "dose_uM", "cellosaur"]):
        mu = XK[g.index.values].mean(0)
        XKc[g.index.values] -= mu
        pm_rows.append(dict(plate=plate, dose_uM=dose, cell_line=cl, n_drugs=int(g.ident.nunique())))
        pm_X.append(mu)
    sizes = {}
    sizes["step2/tahoe_plate_line_mean.parquet"] = write_pq(pd.concat([pd.DataFrame(pm_rows), pd.DataFrame(np.vstack(pm_X).astype(np.float32), columns=genes)], axis=1),
                                                            STEP2 / "tahoe_plate_line_mean.parquet")

    # ---- null + well-level activity z, per plate x dose (all drugs present)
    null_rows = []
    K["z_well"] = np.nan
    for (plate, dose), g in K.groupby(["plate", "dose_uM"]):
        idx = g.index.values
        ident, cl = g.ident.values, g.cellosaur.values
        nd = len(set(ident))
        rec = dict(plate=plate, dose_uM=dose, n_drugs=nd, n_combos=len(idx))
        for l, gl in g.groupby("cellosaur"):          # robust z per line against all drugs on plate x dose
            if gl.ident.nunique() >= MIN_DRUGS_REF:
                K.loc[gl.index, "z_well"] = robust_z(gl.norm.values, gl.norm.values)
        if nd >= MIN_DRUGS_REF:
            R = corr_matrix(XK[idx])
            dc, dl = ident[:, None] != ident[None, :], cl[:, None] != cl[None, :]
            rec.update(null_r_diffdrug_diffline=float(R[dc & dl].mean()), null_r_diffdrug_sameline=float(R[dc & ~dl].mean()),
                       samedrug_r_diffline=float(R[~dc & dl].mean()))
            Rc = corr_matrix(XKc[idx])
            rec.update(centred_null_r_diffdrug_diffline=float(Rc[dc & dl].mean()), centred_samedrug_r_diffline=float(Rc[~dc & dl].mean()))
            del Rc
            by_line = {l: np.where(cl == l)[0] for l in np.unique(cl)}
            fake = []
            for _ in range(200):
                pick, used = [], set()
                for l in rng.permutation(list(by_line)):
                    ok = [j for j in by_line[l] if ident[j] not in used]
                    if ok:
                        j = ok[rng.integers(len(ok))]
                        pick.append(j)
                        used.add(ident[j])
                fake.append(mean_offdiag(R[np.ix_(pick, pick)]))
            rec.update(fake_consistency_mean=float(np.mean(fake)), fake_consistency_p95=float(np.percentile(fake, 95)),
                       fake_consistency_p99=float(np.percentile(fake, 99)))
            # consensus-level activity reference: per drug, cross-line mean on this plate
            ids = sorted(set(ident))
            cn = {d: float(np.linalg.norm(XK[idx[ident == d]].astype(np.float64).mean(0))) for d in ids}
            zc = dict(zip(ids, robust_z([cn[d] for d in ids], list(cn.values()))))
            rec["_cons"] = {d: (cn[d], zc[d]) for d in ids}
        null_rows.append(rec)
    cons_ref = {(r["plate"], r["dose_uM"], d): v for r in null_rows if "_cons" in r for d, v in r["_cons"].items()}
    null = pd.DataFrame([{k: v for k, v in r.items() if k != "_cons"} for r in null_rows])
    null = null.assign(pn=null.plate.astype(int)).sort_values(["dose_uM", "pn"]).drop(columns="pn")
    null.to_csv(STEP2 / "tahoe_null_by_plate.csv", index=False)

    # ---- targets: plate replicates, drug-level per line
    T = K[K.target].copy()
    rep = []
    for _, g in T.groupby(["drug_stripped", "dose_uM", "cellosaur"]):
        if len(g) > 1:
            R = corr_matrix(XK[g.index.values])
            pl = g.plate.values
            rep += [(float(R[a, b]), pl[a] == pl[b]) for a in range(len(g)) for b in range(a + 1, len(g))]
    dl_rows, dl_idx = [], []
    for (d, dose, cl), g in T.groupby(["drug_stripped", "dose_uM", "cellosaur"], sort=True):
        dl_rows.append(dict(drug_stripped=d, dose_uM=dose, cellosaur=cl, cellname=g.cellname.iloc[0], unit=g.unit.iloc[0],
                            plates=sorted(set(g.plate), key=int), n_wells=len(g), n_trt=int(g.n_trt.sum())))
        dl_idx.append(g.index.values)
    DL = pd.DataFrame(dl_rows)
    DLX = np.vstack([XK[ix].mean(0) for ix in dl_idx]).astype(np.float32)
    DLXc = np.vstack([XKc[ix].mean(0) for ix in dl_idx]).astype(np.float32)
    # ---- salt duplicates -> per-line unit signatures
    salt, salt_cd, salt_c, salt_cd_c = [], {}, [], {}
    pl_rows, pl_X, pl_Xc = [], [], []
    for (unit, dose, cl), g in DL.groupby(["unit", "dose_uM", "cellosaur"], sort=True):
        ix = g.index.values
        if len(ix) > 1:
            R = corr_matrix(DLX[ix])
            rs = [float(R[a, b]) for a in range(len(ix)) for b in range(a + 1, len(ix))]
            salt += rs
            salt_cd.setdefault((unit, dose), []).extend(rs)
            Rc_ = corr_matrix(DLXc[ix])
            rc = [float(Rc_[a, b]) for a in range(len(ix)) for b in range(a + 1, len(ix))]
            salt_c += rc
            salt_cd_c.setdefault((unit, dose), []).extend(rc)
        pl_rows.append(dict(unit=unit, dose_uM=dose, cell_line=cl, cell_name=g.cellname.iloc[0],
                            tahoe_drugs=";".join(sorted(g.drug_stripped)),
                            plates=";".join(sorted({p for x in g.plates for p in x}, key=int)),
                            n_wells=int(g.n_wells.sum()), n_cells_trt=int(g.n_trt.sum())))
        pl_X.append(DLX[ix].mean(0))
        pl_Xc.append(DLXc[ix].mean(0))
    PL = pd.DataFrame(pl_rows)
    PLX = np.vstack(pl_X).astype(np.float32)
    PLXc = np.vstack(pl_Xc).astype(np.float32)
    sizes["step2/tahoe_perline.parquet"] = write_pq(pd.concat([PL, pd.DataFrame(PLX, columns=genes)], axis=1), STEP2 / "tahoe_perline.parquet")
    for line, fn in ((LINE_A549, "tahoe_A549.parquet"), (LINE_HEPG2, "tahoe_HepG2C3A.parquet")):
        m = (PL.cell_line == line).values
        sizes[f"step2/{fn}"] = write_pq(pd.concat([PL[m].reset_index(drop=True), pd.DataFrame(PLX[m], columns=genes)], axis=1), STEP2 / fn)

    # ---- consensus, consistency (+bootstrap CI), null per unit
    nl = null.set_index(["plate", "dose_uM"])
    pnorm = np.linalg.norm(PLX.astype(np.float64), axis=1)
    pcnorm = np.linalg.norm(PLXc.astype(np.float64), axis=1)
    tg_by = {k: v.index.values for k, v in T.groupby(["unit", "dose_uM"])}

    def crosswell(cidx, M):
        """mean r over (line i, plate A) x (line j, plate B) pairs, i != j and A != B (no shared well/control)."""
        sub = T.loc[cidx]
        keys, rows = [], []
        for (cl, pl), gg in sub.groupby(["cellosaur", "plate"]):
            keys.append((cl, pl))
            rows.append(M[gg.index.values].mean(0))
        if len({k[1] for k in keys}) < 2:
            return np.nan, 0
        R = corr_matrix(np.vstack(rows))
        cl = np.array([k[0] for k in keys])
        pl = np.array([k[1] for k in keys])
        m = (cl[:, None] != cl[None, :]) & (pl[:, None] != pl[None, :])
        return float(R[m].mean()), len(set(pl))

    cs_rows, cs_X = [], []
    for (unit, dose), g in PL.groupby(["unit", "dose_uM"], sort=True):
        ix = g.index.values
        R = corr_matrix(PLX[ix])
        cw, npl = crosswell(tg_by[(unit, dose)], XK)
        cwc, _ = crosswell(tg_by[(unit, dose)], XKc)
        lo, hi = boot_ci(R, rng)
        plates = sorted({p for x in g.plates for p in x.split(";")}, key=int)

        def nlv(col, f=np.nanmean):
            v = [nl.loc[(p, dose), col] for p in plates if (p, dose) in nl.index and col in nl.columns]
            v = [x for x in v if pd.notna(x)]
            return float(f(v)) if v else np.nan
        cs_rows.append(dict(unit=unit, dose_uM=dose, n_lines=len(ix), plates=";".join(plates), magnitude=float(pnorm[ix].mean()),
                            consistency=mean_offdiag(R), consistency_ci_lo=lo, consistency_ci_hi=hi,
                            consistency_centred=mean_offdiag(corr_matrix(PLXc[ix])), magnitude_centred=float(pcnorm[ix].mean()),
                            consistency_crosswell=cw, consistency_crosswell_centred=cwc, n_plates=npl if npl else len(plates),
                            null_r_sameplate_diffdrug=nlv("null_r_diffdrug_diffline"),
                            null_r_sameplate_diffdrug_sameline=nlv("null_r_diffdrug_sameline"),
                            null_fake_consistency_p95=nlv("fake_consistency_p95", np.nanmax),
                            salt_dup_r=float(np.median(salt_cd[(unit, dose)])) if (unit, dose) in salt_cd else np.nan,
                            salt_dup_r_centred=float(np.median(salt_cd_c[(unit, dose)])) if (unit, dose) in salt_cd_c else np.nan,
                            n_wells=int(g.n_wells.sum()), tahoe_drugs=";".join(sorted({d for x in g.tahoe_drugs for d in x.split(";")}))))
        cs_X.append(PLX[ix].mean(0))
    CS = pd.DataFrame(cs_rows)
    CSX = np.vstack(cs_X).astype(np.float32)
    sizes["tahoe_signatures.parquet"] = write_pq(pd.concat([CS[["unit", "dose_uM", "n_lines", "plates", "magnitude"]],
                                                            pd.DataFrame(CSX, columns=genes)], axis=1), ROOT / "tahoe_signatures.parquet")
    CS["above_null_p95"] = CS.consistency > CS.null_fake_consistency_p95
    CS = CS.merge(U[["unit", "tiers", "moa"]], on="unit", how="left")
    CS.to_csv(STEP2 / "tahoe_consistency.csv", index=False)

    # ---- activity: per well robust z (plate x dose x line, all drugs) -> unit; consensus per plate -> unit
    act = []
    for (unit, dose), g in T.groupby(["unit", "dose_uM"]):
        for ctx, line in (("A549", LINE_A549), ("HepG2C3A", LINE_HEPG2)):
            w = g[g.cellosaur == line]
            m = ((PL.unit == unit) & (PL.dose_uM == dose) & (PL.cell_line == line)).values
            if not m.any():
                continue
            z = w.z_well.dropna()
            act.append(dict(unit=unit, dose_uM=dose, context=ctx, cell_line=line, n_wells=len(w), n_wells_scored=len(z),
                            n_cells_trt=int(w.n_trt.sum()),
                            z_norm=float(z.mean()) if len(z) else np.nan, n_wells_active=int((z > 3).sum()),
                            norm=float(np.linalg.norm(PLX[m][0].astype(np.float64))),
                            frac_abs_lfc_gt1=float((np.abs(PLX[m][0]) > 1).mean()), plates=";".join(sorted(set(w.plate), key=int))))
        keys = sorted({(p, dose, unit) for p in g.plate})
        zs = [cons_ref[k][1] for k in keys if k in cons_ref and pd.notna(cons_ref[k][1])]
        m = ((CS.unit == unit) & (CS.dose_uM == dose)).values
        act.append(dict(unit=unit, dose_uM=dose, context="consensus", cell_line="", n_wells=len(keys), n_wells_scored=len(zs),
                        n_cells_trt=int(g.n_trt.sum()),
                        z_norm=float(np.mean(zs)) if zs else np.nan, n_wells_active=int(sum(z > 3 for z in zs)),
                        norm=float(np.linalg.norm(CSX[m][0].astype(np.float64))),
                        frac_abs_lfc_gt1=float((np.abs(CSX[m][0]) > 1).mean()), plates=";".join(sorted({k[0] for k in keys}, key=int))))
    A = pd.DataFrame(act)
    A["active"] = A.z_norm > 3
    # does the norm track the treated-cell count? (sampling noise grows as cells fall)
    from scipy.stats import spearmanr
    rho = [spearmanr(g.norm, g.n_trt).statistic for _, g in K.groupby(["plate", "dose_uM", "cellosaur"]) if len(g) >= 10]
    cellcheck = dict(median_within_group_spearman_norm_vs_ntrt=float(np.nanmedian(rho)), n_groups=len(rho),
                     median_n_trt_per_well={f"{c}|{'active' if act_ else 'inactive'}": float(gg.n_cells_trt.median() / max(1, gg.n_wells.median()))
                                            for (c, act_), gg in A.groupby(["context", "active"])},
                     median_n_trt_unit_sum={f"{c}|{'active' if act_ else 'inactive'}": float(gg.n_cells_trt.median()) for (c, act_), gg in A.groupby(["context", "active"])})
    A.to_csv(STEP2 / "tahoe_activity.csv", index=False)

    # ---- summary
    piv = CS.pivot(index="unit", columns="dose_uM", values="magnitude")
    found = {str(d): sorted(set(T[T.dose_uM == d].drug_stripped)) for d in DOSES3}
    strong = CS[CS.moa.fillna("").str.contains("Proteasome|HDAC|MTOR|mTOR|MEK|CDK|Topoisomerase|HSP90|tubulin|Tubulin", regex=True)]
    summ = dict(
        threshold=thr, n_genes_candidates=len(gnames), n_genes=len(genes),
        gene_nanfrac_quantiles={str(q): float(np.quantile(nanfrac, q)) for q in (.5, .9, .99, 1)},
        units_inconsistent_drug_lists=U[~U.consistent].unit.tolist(), n_units=len(U), n_target_names=len(names),
        names_found_per_dose={d: len(v) for d, v in found.items()},
        names_missing_per_dose={d: sorted(set(names) - set(v)) for d, v in found.items()},
        plate_rep_r_median=float(np.median([r for r, s in rep if not s])) if rep else None,
        plate_rep_n=int(sum(not s for _, s in rep)), same_plate_rep_n=int(sum(s for _, s in rep)),
        salt_dup_r_median=float(np.median(salt)) if salt else None,
        salt_dup_r_centred_median=float(np.median(salt_c)) if salt_c else None, activity_cellcount_check=cellcheck, salt_dup_n_pairs=len(salt), salt_dup_n_unit_dose=len(salt_cd),
        n_perline_rows=len(PL), perline_rows_per_dose={str(k): int(v) for k, v in PL.dose_uM.value_counts().sort_index().items()},
        n_consensus_rows=len(CS), units_per_dose={str(k): int(v) for k, v in CS.dose_uM.value_counts().sort_index().items()},
        n_lines_per_unit_dose={str(d): CS[CS.dose_uM == d].n_lines.describe().round(2).to_dict() for d in DOSES3},
        A549_rows_per_dose={str(k): int(v) for k, v in PL[PL.cell_line == LINE_A549].dose_uM.value_counts().sort_index().items()},
        HepG2_rows_per_dose={str(k): int(v) for k, v in PL[PL.cell_line == LINE_HEPG2].dose_uM.value_counts().sort_index().items()},
        magnitude_median_per_dose={str(d): float(CS[CS.dose_uM == d].magnitude.median()) for d in DOSES3},
        frac_5_gt_0p5=float((piv[5.0] > piv[0.5]).mean()), frac_0p5_gt_0p05=float((piv[0.5] > piv[0.05]).mean()),
        consistency_quantiles={str(d): CS[CS.dose_uM == d].consistency.quantile([0, .1, .25, .5, .75, .9, 1]).round(4).to_dict() for d in DOSES3},
        frac_above_null_p95={str(d): float(CS[CS.dose_uM == d].above_null_p95.mean()) for d in DOSES3},
        consistency_centred_median={str(d): float(CS[CS.dose_uM == d].consistency_centred.median()) for d in DOSES3},
        consistency_crosswell_centred_median={str(d): float(CS[CS.dose_uM == d].consistency_crosswell_centred.median()) for d in DOSES3},
        n_crosswell={str(d): int(CS[CS.dose_uM == d].consistency_crosswell.notna().sum()) for d in DOSES3},
        magnitude_centred_median_per_dose={str(d): float(CS[CS.dose_uM == d].magnitude_centred.median()) for d in DOSES3},
        frac_5_gt_0p5_centred=float((CS.pivot(index="unit", columns="dose_uM", values="magnitude_centred")[5.0] > CS.pivot(index="unit", columns="dose_uM", values="magnitude_centred")[0.5]).mean()),
        z_norm_median={f"{c}|{d}": float(g.z_norm.median()) for (c, d), g in A.groupby(["context", "dose_uM"])},
        dose_trend_active_units=(lambda pa: dict(n_units=len(pa), frac_5_gt_0p5=float((pa[5.0] > pa[0.5]).mean()),
                                                  frac_0p5_gt_0p05=float((pa[0.5] > pa[0.05]).mean()),
                                                  median_magnitude={str(d): float(pa[d].median()) for d in DOSES3}))(
            piv.dropna().loc[lambda x: x.index.isin(set(A[(A.context == "consensus") & A.active].unit))]),
        frac_any_context_active_per_dose={str(d): float(A[A.dose_uM == d].groupby("unit").active.any().mean()) for d in DOSES3},
        centred_null_r_by_dose={str(d): float(null[null.dose_uM == d].centred_null_r_diffdrug_diffline.mean()) for d in DOSES3},
        centred_samedrug_r_by_dose={str(d): float(null[null.dose_uM == d].centred_samedrug_r_diffline.mean()) for d in DOSES3},
        null_r_diffdrug_diffline_by_dose={str(d): float(null[null.dose_uM == d].null_r_diffdrug_diffline.mean()) for d in DOSES3},
        null_r_diffdrug_sameline_by_dose={str(d): float(null[null.dose_uM == d].null_r_diffdrug_sameline.mean()) for d in DOSES3},
        strong_moa_consistency_median={str(d): float(strong[strong.dose_uM == d].consistency.median()) for d in DOSES3},
        all_consistency_median={str(d): float(CS[CS.dose_uM == d].consistency.median()) for d in DOSES3},
        activity_frac_active={f"{c}|{d}": float(g.active.mean()) for (c, d), g in A.groupby(["context", "dose_uM"])},
        activity_n={f"{c}|{d}": int(len(g)) for (c, d), g in A.groupby(["context", "dose_uM"])},
        activity_frac_gt1_median={f"{c}|{d}": float(g.frac_abs_lfc_gt1.median()) for (c, d), g in A.groupby(["context", "dose_uM"])},
        output_bytes=sizes, build_s=round(time.time() - t0, 1))
    (CACHE / "build_summary.json").write_text(json.dumps(summ, indent=1, default=str))
    # ---- download summary over the whole pipeline
    log = pd.read_csv(CACHE / "fetch_log.csv")
    d = log[log.status == "done"]
    pil = json.loads((CACHE / "pilot_info.json").read_text())
    byph = {ph: dict(footer_bytes=int(g.footer_bytes.fillna(0).sum()), data_bytes=int(g.data_bytes.fillna(0).sum()),
                     n_files=int(g.file.nunique()), n_requests=int(g.n_requests.fillna(0).sum()), retries=int(g.retries.fillna(0).sum()),
                     wall_s=float(g.drop_duplicates("run_t0").run_wall_s.sum()))
            for ph, g in d.assign(ph=np.where(d.phase == "fetch2", "fetch2", "fetch_v1")).groupby("ph")}
    total = pil["bytes"] + sum(v["footer_bytes"] + v["data_bytes"] for v in byph.values())
    disk = sum(p.stat().st_size for p in CACHE.rglob("*") if p.is_file())
    ds = dict(total_bytes_streamed=int(total), total_GB=round(total / 1e9, 3), pilot_bytes=pil["bytes"], phases=byph,
              note="excludes ~0.1 GB of interactive layout probes before the pipeline existed; no DE file was downloaded in full",
              peak_disk_tahoe_cache_bytes=int(disk), peak_disk_note="cache only grows (nothing deleted), so its final size is the peak",
              output_bytes=sizes)
    (STEP2 / "tahoe_download_summary.json").write_text(json.dumps(ds, indent=1))
    print(json.dumps(summ, indent=1, default=str))
    print(json.dumps(ds, indent=1))
    print(line_tab.head(8).to_string())
    print(CS.sort_values("consistency", ascending=False)[["unit", "dose_uM", "n_lines", "consistency", "magnitude", "moa", "tahoe_drugs"]].head(25).to_string())

def perplate(args):
    """Per-plate (not plate-averaged) signatures for unit x dose observed on >=2 plates; cache only, no download."""
    t0 = time.time()
    combos = pd.read_parquet(CACHE / "de_index_combos.parquet").sort_values("combo").reset_index(drop=True)
    pos, gnames = load_candidates()
    gh = gene_hash(pos)
    genes = [c for c in pq.read_schema(STEP2 / "tahoe_perline.parquet").names if c in set(gnames)]
    gsel = np.array([gnames.index(g) for g in genes])
    U, name2unit, names = load_units()
    combos = combos.merge(name2unit, on="drug_stripped", how="inner")
    combos = combos[combos.n_trt >= MIN_CELLS]
    npl = combos.groupby(["unit", "dose_uM"]).plate.nunique()
    multi = npl[npl >= 2].index
    sel = combos.set_index(["unit", "dose_uM"]).loc[multi].reset_index().sort_values("combo").reset_index(drop=True)
    row_of = pd.Series(np.arange(len(sel)), index=sel.combo.values)
    X = np.full((len(sel), len(pos)), np.nan, np.float32)
    C = np.zeros((len(sel), len(pos)), np.int8)
    files = sorted({int(f) for fs in sel.files for f in fs.split(";")})
    files = sorted(set(files) | {f - 1 for f in files if f > 0} | {f + 1 for f in files if f < N_FILES - 1})
    for f in files:
        for suf in ("frag", "frag2"):
            pth = FILES_DIR / f"{f:04d}_{suf}.npz"
            if not pth.exists():
                continue
            z = np.load(pth)
            if str(z["gene_hash"]) != gh:
                continue
            for k, cc in enumerate(z["combos"]):
                r = row_of.get(int(cc))
                if r is not None:
                    cv = z["cov"][k]
                    X[r, cv] = z["frag"][k, cv]
                    C[r] += cv
    assert (C == 1).all(), "cache does not cover every needed combination exactly once"
    X = np.nan_to_num(X[:, gsel], nan=0.0)
    rows, PX = [], []
    for (unit, dose, plate, cl), g in sel.groupby(["unit", "dose_uM", "plate", "cellosaur"], sort=True):
        rows.append(dict(unit=unit, dose_uM=dose, plate=plate, cell_line=cl, n_cells_trt=int(g.n_trt.sum()),
                         n_drugs=len(g), tahoe_drugs=";".join(sorted(g.drug_stripped))))
        PX.append(X[g.index.values].mean(0))       # salt variants on the same plate are averaged
    P = pd.DataFrame(rows)
    PX = np.vstack(PX).astype(np.float32)
    P = P.assign(pn=P.plate.astype(int)).sort_values(["unit", "dose_uM", "pn", "cell_line"]).drop(columns="pn")
    PX = PX[P.index.values]
    P = P.reset_index(drop=True)
    sz1 = write_pq(pd.concat([P[["unit", "dose_uM", "plate", "cell_line", "n_cells_trt"]], pd.DataFrame(PX, columns=genes)], axis=1),
                   STEP2 / "tahoe_perplate.parquet")
    crow, CX = [], []
    for (unit, dose, plate), g in P.groupby(["unit", "dose_uM", "plate"], sort=False):
        crow.append(dict(unit=unit, dose_uM=dose, plate=plate, n_lines=len(g), n_cells_trt=int(g.n_cells_trt.sum())))
        CX.append(PX[g.index.values].mean(0))
    CP = pd.DataFrame(crow)
    sz2 = write_pq(pd.concat([CP, pd.DataFrame(np.vstack(CX).astype(np.float32), columns=genes)], axis=1),
                   STEP2 / "tahoe_perplate_consensus.parquet")
    up = CP.groupby("dose_uM").unit.nunique()
    summ = dict(n_rows_perplate=len(P), n_rows_consensus=len(CP), units_per_dose={str(k): int(v) for k, v in up.items()},
                unit_dose_by_n_plates={str(k): int(v) for k, v in CP.groupby(["unit", "dose_uM"]).size().value_counts().sort_index().items()},
                rows_per_dose={str(k): int(v) for k, v in P.dose_uM.value_counts().sort_index().items()},
                n_salt_merged_rows=int((P.n_drugs > 1).sum()) if "n_drugs" in P else 0,
                lines_per_unit_dose_plate=CP.n_lines.describe().round(1).to_dict(), n_genes=len(genes),
                bytes={"step2/tahoe_perplate.parquet": sz1, "step2/tahoe_perplate_consensus.parquet": sz2},
                min_cells=MIN_CELLS, centred=False, seconds=round(time.time() - t0, 1))
    (CACHE / "perplate_summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all", choices=["pilot", "fetch", "fetch2", "index", "build", "perplate", "all"])
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--files", type=str, default=None, help="debug: run fetch2 on a comma list of files")
    args = ap.parse_args()
    if args.files:
        for i in map(int, args.files.split(",")):
            print(process_file2(i, args.threads))
        sys.exit()
    if args.stage in ("pilot", "all") and not (CACHE / "candidate_genes.txt").exists():
        pilot(args)
    if args.stage == "fetch":
        fetch(args)
    if args.stage in ("fetch2", "all"):
        fetch2(args)
    if args.stage in ("index", "all"):
        index(args)
    if args.stage in ("build", "all"):
        build(args)
    if args.stage in ("perplate", "all"):
        perplate(args)
