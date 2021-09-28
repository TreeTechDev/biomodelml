import modin.pandas as pandas
import sys
import os
from multiprocessing.dummy import Pool

directory = sys.argv[1]
output_path = os.path.join(directory, "variance.txt")
output_all_path = os.path.join(directory, "gene_cluster.parquet")
n_threads = 5
to_load = []
df = pandas.DataFrame()

def read_csv(filename):
    df = pandas.read_csv(os.path.join(directory, filename), index_col="cluster")
    df.rename(
        columns={"prob": filename.split(".csv")[0].replace("NC_000913.3:", "")},
        inplace=True,
        copy=False
    )
    return df

for i, filename in enumerate(os.listdir(directory)):
    if filename.endswith(".csv") and not output_path.endswith(filename):
        print(f"adding file {i+1} to df...")
        to_load.append(filename)
        if not (i+1) % n_threads:
            with Pool(len(to_load)) as pool:
                dfs = list(pool.map(read_csv, to_load))
            to_load = []
            df = pandas.concat(df + dfs, axis=1)
            print(f"new shape: {df.shape}")

df = df.T.fillna(0.0)
df.var().to_csv(output_path)
print("Saving compressed...")
df.to_parquet(output_all_path)