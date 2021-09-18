import pandas
import sys
import os

directory = sys.argv[1]
output_path = os.path.join(directory, "variance.txt")
output_all_path = os.path.join(directory, "gene_cluster.parquet")

df = pandas.DataFrame()

for i, filename in enumerate(os.listdir(directory)):
    if filename.endswith(".csv") and not output_path.endswith(filename):
        print(f"adding file {i+1} to list...")
        df = pandas.concat(
            [df, pandas.read_csv(os.path.join(directory, filename), index_col="cluster")],
        copy=False, index=1)

df.T.var().to_csv(output_path)
print("Saving compressed...")
df.to_parquet(output_all_path)