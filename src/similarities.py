import sys
import pandas
import vaex
import multiprocessing
from scipy.spatial.distance import jensenshannon

file_path = sys.argv[1]
chunksize = 30_000
processes = int(1)
manager = multiprocessing.Manager()
output_values = manager.list()


def similarity(df, column, compare):
    print(f"doing similarities between {column} and {compare}")
    similarity = 1.0 - jensenshannon(df[column].values, df[compare].values)
    output = f"\n{column},{compare},{similarity}"
    if column != compare:
       output += f"\n{compare},{column},{similarity}"
    output_values.append(output)

df = vaex.open(file_path)

with open(f"{file_path}.similarity.csv", "w") as sim:
    sim.write("column,compare,similarity")

columns = list(df.columns)[1:]
args = [(df, column, compare) for i, column in enumerate(columns) for compare in columns[i:]]

print(f"starting similarity check with {len(args)} combinations...")
with multiprocessing.Pool(processes) as pool:
    pool.starmap(
        similarity,
        args
    )

print("writing similarities to file...")
with open(f"{file_path}.similarity.csv", "a") as sim:
    sim.writelines(output_values)
