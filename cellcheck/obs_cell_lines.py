"""List the screened Tahoe-100M cell lines from obs_metadata.parquet without downloading it:
reads the pandas schema metadata (footer) and the dictionary-encoded `cell_line`/`cell_name` columns
of single row groups via HTTP range requests (~0.8 MB per column chunk)."""
import json, sys, pathlib, fsspec, pyarrow.parquet as pq
URL = "https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/metadata/obs_metadata.parquet"
fs = fsspec.filesystem("https")
rgs = [int(x) for x in sys.argv[1:]] or [0]
with fs.open(URL, "rb", block_size=2**20) as f:
    pf = pq.ParquetFile(f)
    pdmeta = json.loads(pf.schema_arrow.metadata[b"pandas"])
    for c in pdmeta["columns"]:
        if c["name"] in ("cell_line", "cell_name"):
            print(c["name"], c.get("metadata"))
    seen = {}
    for rg in rgs:
        t = pf.read_row_group(rg, columns=["cell_line", "cell_name"])
        for col in ("cell_line", "cell_name"):
            ch = t.column(col).chunk(0)
            print(f"rg{rg} {col}: dictionary size {len(ch.dictionary)}")
        pairs = t.to_pandas().astype(str).drop_duplicates()
        print(f"rg{rg}: {len(pairs)} distinct (cell_line, cell_name) pairs present")
        for a, b in pairs.itertuples(index=False):
            seen[a] = b
    d1 = t.column("cell_line").chunk(0).dictionary.to_pylist()
    d2 = t.column("cell_name").chunk(0).dictionary.to_pylist()
    json.dump({"cell_line_dictionary": d1, "cell_name_dictionary": d2, "pairs_seen": seen},
              open(pathlib.Path(__file__).parent / "data/obs_cell_lines.json", "w"), indent=1)
    print(sorted(seen.items(), key=lambda x: x[1].lower()))
