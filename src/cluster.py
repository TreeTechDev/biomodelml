import numpy
import pandas
import hdbscan
import multiprocessing
import os
import sys
import pyarrow
import pyarrow.parquet as pq

directory = sys.argv[1]
processes = int(multiprocessing.cpu_count()-1)
manager = multiprocessing.Manager()
output_dfs = manager.list()
parquet_file = os.path.join(directory, "cluster_matrix.parquet")

def clusterize(path: str, filename: str) -> pandas.DataFrame:
    if not filename.endswith(".csv"):
        return
    csv_file = os.path.join(path, filename)
    basename = filename.split(".csv")[0]
    df = pandas.read_csv(csv_file).set_index("sequence")
    clusterer = hdbscan.HDBSCAN(
        cluster_selection_method="leaf",
        min_samples=2,
        min_cluster_size=3,
        prediction_data=True,
        approx_min_span_tree=False,
        metric="manhattan")
    clusterer.fit(df)
    soft_clusters = hdbscan.all_points_membership_vectors(clusterer)
    columns = soft_clusters.shape[1]
    score_df = pandas.DataFrame(
        soft_clusters.T,
        index=[f"{basename}_{c}" for c in range(columns)],
        columns=[c.replace("(", "").replace(")", "") for c in df.index]
    )
    score_df.index.name = "index"
    output_dfs.append(score_df)
    print(f"using file {csv_file} with {clusterer.labels_.max()} clusters")

args = [(directory, filename) for filename in os.listdir(directory)]

print(f"running clusteres with {len(args)} files...")

with multiprocessing.Pool(processes) as pool:
    pool.starmap(
        clusterize,
        args
    )

print("writing parquet with clusteres...")
parquet_schema = pyarrow.Table.from_pandas(df=output_dfs[0], preserve_index=True).schema
parquet_writer = pq.ParquetWriter(parquet_file, parquet_schema, compression='snappy')

for result_df in output_dfs:
    table = pyarrow.Table.from_pandas(result_df, schema=parquet_schema, preserve_index=True)
    parquet_writer.write_table(table)
parquet_writer.close()