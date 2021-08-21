import vaex
import sys
import os
import multiprocessing
from pandas.errors import EmptyDataError

processes = int(multiprocessing.cpu_count()-1)
directory = sys.argv[1]
args = []
for filename in os.listdir(directory):
    if filename.endswith(".csv") and not os.path.exists(f"{os.path.join(directory, filename)}.hdf5"):
        args.append((directory, filename))

def convert(directory, filename):
    filepath = os.path.join(directory, filename)
    try:
        vaex.from_csv(filepath, convert=True, copy_index=True)
    except EmptyDataError:
        print(f"no data on file {filename}")

print(f"converting {len(args)} files to hdf5...")

with multiprocessing.Pool(processes) as pool:
    pool.starmap(
        convert,
        args
    )