"""Resolve Cell Painting Gallery accessions from the gallery's own S3 index (no hardcoded accession ids).

    from gallery_index import resolve, s3_list, s3_url, Bytes
    acc = resolve("lincs")          # -> e.g. 'cpg0004-lincs'
    prefixes, keys = s3_list(f"{acc}/broad/workspace/", delimiter="/")

Fails loudly (RuntimeError) if the index is unreachable or a dataset can't be resolved uniquely.
Also provides a process-wide byte counter (Bytes) so every script can report bytes downloaded.
"""
import re
import threading
import xml.etree.ElementTree as ET

import requests

BUCKET = "https://cellpainting-gallery.s3.amazonaws.com"
NS = "{http://s3.amazonaws.com/doc/2006-03-01/}"
# keyword patterns matched against the gallery's top-level accession prefixes
DATASETS = {
    "lincs": r"^cpg\d{4}-lincs$",                       # LINCS Cell Painting (A549)
    "cpjump1": r"^cpg\d{4}-jump-pilot$",                # CPJUMP1 pilot (A549 + U2OS)
    "euos": r"^cpg\d{4}-eu-os-bioactives$",            # EU-OPENSCREEN bioactives (HepG2 + U2OS)
    "jump": r"^cpg\d{4}-jump$",                         # main JUMP production set (U2OS)
}


class _Bytes:
    def __init__(self):
        self.n = 0
        self._l = threading.Lock()

    def add(self, k):
        with self._l:
            self.n += int(k)


Bytes = _Bytes()


def s3_list(prefix="", delimiter="/", max_pages=1000):
    """List an S3 prefix anonymously. Returns (common_prefixes, [(key, size)])."""
    prefixes, keys, token = [], [], None
    for _ in range(max_pages):
        params = {"list-type": "2", "prefix": prefix, "delimiter": delimiter}
        if token:
            params["continuation-token"] = token
        r = requests.get(BUCKET + "/", params=params, timeout=60)
        if r.status_code != 200:
            raise RuntimeError(f"Cell Painting Gallery index unreachable ({r.status_code}) for prefix {prefix!r}")
        Bytes.add(len(r.content))
        root = ET.fromstring(r.content)
        prefixes += [p.find(NS + "Prefix").text for p in root.findall(NS + "CommonPrefixes")]
        keys += [(c.find(NS + "Key").text, int(c.find(NS + "Size").text)) for c in root.findall(NS + "Contents")]
        if root.find(NS + "IsTruncated").text != "true":
            return prefixes, keys
        token = root.find(NS + "NextContinuationToken").text
    raise RuntimeError("S3 listing did not terminate")


def accessions():
    prefixes, _ = s3_list("", "/")
    return [p.rstrip("/") for p in prefixes]


def resolve(name):
    pat = re.compile(DATASETS[name], re.I)
    hits = [a for a in accessions() if pat.match(a)]
    if len(hits) != 1:
        raise RuntimeError(f"could not resolve gallery accession for {name!r}: {hits}")
    return hits[0]


def s3_url(key):
    return f"{BUCKET}/{key}"


if __name__ == "__main__":
    for k in DATASETS:
        print(k, "->", resolve(k))
