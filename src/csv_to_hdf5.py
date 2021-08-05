import vaex
import sys
import os
import multiprocessing

processes = int(multiprocessing.cpu_count()-1)
directory = sys.argv[1]
args = [(directory, filename) for filename in os.listdir(directory)]

def convert(filename):
    vaex.from_csv(filename, convert=True, copy_index=True)

print(f"converting {len(args)} files to hdf5...")

with multiprocessing.Pool(processes) as pool:
    pool.starmap(
        convert,
        args
    )

print("grouping to one hdf5...")
df = vaex.open(os.path.join(directory, "*.csv.hdf5"))
df.export(os.path.join(directory, "clusters_combined.hdf5"))