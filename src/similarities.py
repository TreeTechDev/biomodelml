import sys
import os
import vaex
import multiprocessing
from scipy.spatial.distance import jensenshannon

directory = sys.argv[1]
output_path = os.path.join(directory, "similarity.csv")
chunksize = 30_000
processes = 1  # need more memory


def similarity(df, column, compare):
    print(f"doing similarities between {column} and {compare}")
    similarity = df.mutual_information(column, compare)[0]
    output = f"\n{column},{compare},{similarity}"
    if column != compare:
       output += f"\n{compare},{column},{similarity}"
    with open(output_path, "w") as sim:
        sim.write(output)
        sim.flush()
    print(f"wrote {column} and {compare} similarity {similarity} to file")

print(f"openning all hdf5 files...")

df = vaex.open(os.path.join(directory, "*.csv.hdf5"))

with open(output_path, "w") as sim:
    sim.write("column,compare,similarity")

columns = list(df.columns)[1:]
args = [(df, column, compare) for i, column in enumerate(columns) for compare in columns[i:]]

print(f"starting similarity check with {len(args)} combinations...")
with multiprocessing.Pool(processes) as pool:
    pool.starmap(
        similarity,
        args
    )
