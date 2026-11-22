"""Minimal decoder for one Parquet column chunk of a nullable FLOAT column (dictionary or plain encoded,
SNAPPY/uncompressed, data page v1 or v2). Lets tahoe_signatures.py decode log2FoldChange chunks from
byte offsets cached in its row-group index, without re-reading the 5.75 MB file footers.
Validated against pyarrow on fetched fragments (see tahoe_signatures.py --stage validate_decoder)."""
import numpy as np
import pyarrow as pa

_SNAPPY = pa.Codec("snappy")


def _varint(b, i):
    r = s = 0
    while True:
        x = b[i]
        i += 1
        r |= (x & 0x7F) << s
        s += 7
        if not x & 0x80:
            return r, i


def _zz(n):
    return (n >> 1) ^ -(n & 1)


def _val(b, i, t):
    if t == 1:
        return True, i
    if t == 2:
        return False, i
    if t == 3:
        return b[i], i + 1
    if t in (4, 5, 6):
        v, i = _varint(b, i)
        return _zz(v), i
    if t == 7:
        return b[i:i + 8], i + 8
    if t == 8:
        n, i = _varint(b, i)
        return b[i:i + n], i + n
    if t in (9, 10):
        h = b[i]
        i += 1
        n, et = h >> 4, h & 0x0F
        if n == 15:
            n, i = _varint(b, i)
        out = []
        for _ in range(n):
            if et in (1, 2):          # bool list elements are one byte each
                out.append(b[i] == 1)
                i += 1
            else:
                v, i = _val(b, i, et)
                out.append(v)
        return out, i
    if t == 12:
        return _struct(b, i)
    raise ValueError(f"thrift type {t} unsupported")


def _struct(b, i):
    out, fid = {}, 0
    while True:
        h = b[i]
        i += 1
        if h == 0:
            return out, i
        t, d = h & 0x0F, h >> 4
        if d:
            fid += d
        else:
            v, i = _varint(b, i)
            fid = _zz(v)
        out[fid], i = _val(b, i, t)


def _rle_hybrid(buf, i, w, n):
    out = np.zeros(n, np.int64)
    if w == 0:
        return out, i
    k, bw = 0, (w + 7) // 8
    weights = (1 << np.arange(w, dtype=np.int64))
    while k < n:
        h, i = _varint(buf, i)
        if h & 1:
            groups = h >> 1
            nb = groups * w
            bits = np.unpackbits(np.frombuffer(buf, np.uint8, nb, i), bitorder="little")
            vals = bits.reshape(-1, w).astype(np.int64) @ weights
            take = min(groups * 8, n - k)
            out[k:k + take] = vals[:take]
            i += nb
        else:
            run = h >> 1
            v = int.from_bytes(buf[i:i + bw], "little")
            i += bw
            take = min(run, n - k)
            out[k:k + take] = v
        k += take
    return out, i


def _decompress(codec, data, usize):
    if codec == "snappy":
        return _SNAPPY.decompress(data, decompressed_size=usize, asbytes=True)
    return bytes(data)


def decode_float_chunk(buf, num_rows, codec="snappy"):
    """buf: bytes of the whole column chunk (may have trailing slack). Returns float32[num_rows], NaN for nulls."""
    mv = memoryview(buf)
    i, dic, parts, got = 0, None, [], 0
    while got < num_rows:
        hdr, i = _struct(mv, i)
        ptype, usize, csize = hdr[1], hdr[2], hdr[3]
        body = mv[i:i + csize]
        i += csize
        if ptype == 2:                                   # dictionary page
            dec = _decompress(codec, body, usize)
            dic = np.frombuffer(dec, "<f4", hdr[7][1])
            continue
        if ptype == 0:                                   # data page v1
            dph = hdr[5]
            nv, enc = dph[1], dph[2]
            dec = _decompress(codec, body, usize)
            j = 0
            dlen = int.from_bytes(dec[0:4], "little")
            defs, _ = _rle_hybrid(dec, 4, 1, nv)
            j = 4 + dlen
        elif ptype == 3:                                 # data page v2
            dph = hdr[8]
            nv, enc, dl, rl = dph[1], dph[4], dph[5], dph.get(6, 0)
            lev = bytes(body[:dl + rl])
            defs, _ = _rle_hybrid(lev, rl, 1, nv)
            rest = body[dl + rl:]
            comp = dph.get(7, True)
            dec = _decompress(codec, rest, usize - dl - rl) if comp else bytes(rest)
            j = 0
        else:
            raise ValueError(f"page type {ptype}")
        nn = int(defs.sum())
        if enc in (8, 2):                                # RLE_DICTIONARY / PLAIN_DICTIONARY
            w = dec[j]
            idx, _ = _rle_hybrid(dec, j + 1, w, nn)
            vals = dic[idx]
        elif enc == 0:                                   # PLAIN
            vals = np.frombuffer(dec, "<f4", nn, j)
        else:
            raise ValueError(f"encoding {enc}")
        out = np.full(nv, np.nan, np.float32)
        out[defs.astype(bool)] = vals
        parts.append(out)
        got += nv
    return np.concatenate(parts)[:num_rows]
