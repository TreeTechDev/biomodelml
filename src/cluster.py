import numpy
import ray
ray.init()
import modin.pandas as pandas
import hdbscan
import multiprocessing
import os
import sys

directory = sys.argv[1]
processes = int(multiprocessing.cpu_count()-1)
output_path = sys.argv[2]

def clusterize(path: str, filename: str) -> pandas.DataFrame:
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
    lines, columns = soft_clusters.shape
    header = ["index"] + list(df.index)
    output_name = os.path.join(output_path, filename)
    
    matrix = numpy.zeros((columns,lines+1)).astype("object")
    matrix[:, 1:] = soft_clusters.T
    matrix[:, 0] = [f"{basename}_{c}" for c in range(columns)]

    numpy.savetxt(
        output_name, matrix, delimiter=",", header=",".join(header), fmt="%s", comments="")
    print(
        f"using file {csv_file} to {output_name} with {clusterer.labels_.max()} clusters")

args = []
for filename in os.listdir(directory):
    if filename.endswith(".csv"):
        args.append((directory, filename))

print(f"running clusteres with {len(args)} files...")

with multiprocessing.Pool(processes) as pool:
    pool.starmap(
        clusterize,
        args
    )