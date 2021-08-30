import sys
import os
import pandas
import multiprocessing
from scipy.spatial.distance import jensenshannon

directory = sys.argv[1]
output_path = os.path.join(directory, "similarity.csv")
processes = int(multiprocessing.cpu_count()-1)
args = []

def similarity(df, filename, filename_compare):
    basename = filename.split(".csv")[0]
    basename_compare = filename_compare.split("csv")[0]
    filepath = os.path.join(directory, filename_compare)
    df_compare = pandas.read_csv(
        filepath, index_col=0, header=None).T
    df_join = pandas.concat(
        [df, df_compare], join="inner", ignore_index=True, copy=False).drop("cluster", axis=1)
    print(f"doing similarities between {basename} and {basename_compare} with {len(df_join.columns)} features")
    try:
        similarity = 1 - jensenshannon(*df_join.values)
    except Exception as e:
        print(f"fail to compute similarity between {basename} and {basename_compare}")
        print(e)
        return
    output = f"\n{basename},{basename_compare},{similarity}"
    if basename != basename_compare:
       output += f"\n{basename_compare},{basename},{similarity}"
    with open(output_path, "a") as sim:
        sim.write(output)
        sim.flush()
    print(f"wrote {basename} and {basename_compare} similarity {similarity} to file")

print(f"openning all csv files...")

for filename in os.listdir(directory):
    if filename.endswith(".csv") and not output_path.endswith(filename):
        args.append((directory, filename))

with open(output_path, "w") as sim:
    sim.write("gene,compare,similarity")


print(f"starting similarity check with {len(args)} combinations...")

for i, path in enumerate(args):
    filepath = os.path.join(*path)
    df = pandas.read_csv(filepath, index_col=0, header=None).T
    compare_args = [(df, path[1], f) for d, f in args[i:]]
    with multiprocessing.Pool(processes) as pool:
        pool.starmap(
            similarity,
            compare_args
        )
