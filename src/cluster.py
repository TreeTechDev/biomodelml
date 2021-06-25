import pandas
import hdbscan
import multiprocessing
import os
import sys

directory = sys.argv[1]
processes = int(multiprocessing.cpu_count()/2)
manager = multiprocessing.Manager()
output_dfs = manager.list()

def clusterize(path: str, filename: str) -> pandas.DataFrame:
    print(filename)
    csv_file = os.path.join(path, filename)
    basename = csv_file.split(".csv")[0]
    df = pandas.read_csv(csv_file).set_index("sequence")
    clusterer = hdbscan.HDBSCAN(
        cluster_selection_method="leaf",
        min_samples=8,
        min_cluster_size=5,
        prediction_data=True,
        approx_min_span_tree=False,
        metric="manhattan")
    clusterer.fit(df)
    soft_clusters = hdbscan.all_points_membership_vectors(clusterer)
    columns = soft_clusters.shape[1]
    score_df = pandas.DataFrame(
        soft_clusters.T,
        index=[f"{basename}_{c}" for c in range(columns)],
        columns=df.index
    )
    output_dfs.append(score_df)

args = [(directory, filename) for filename in os.listdir(directory)]
with multiprocessing.Pool(processes) as pool:
    pool.starmap(
        clusterize,
        args
    )

result_df = pandas.concat(output_dfs, copy=False)

result_df.to_csv(os.path.join(directory, "cluster_matrix.csv"))