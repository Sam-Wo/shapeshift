"""Read only the parquet footer of Tahoe-100M obs_metadata.parquet (2.3 GB) over HTTP range requests."""
import fsspec, pyarrow.parquet as pq
URL = "https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/metadata/obs_metadata.parquet"
fs = fsspec.filesystem("https")
with fs.open(URL, "rb", block_size=2**20) as f:
    pf = pq.ParquetFile(f)
    md = pf.metadata
    print(md.num_rows, md.num_row_groups, pf.schema_arrow)
    rg = md.row_group(0)
    for j in range(rg.num_columns):
        c = rg.column(j)
        print(c.path_in_schema, c.encodings, c.statistics.has_min_max if c.statistics else None,
              c.statistics.min if c.statistics and c.statistics.has_min_max else '', c.statistics.max if c.statistics and c.statistics.has_min_max else '',
              c.total_compressed_size, c.dictionary_page_offset)
