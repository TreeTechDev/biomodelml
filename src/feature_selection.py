import sys
import os
import pickle
import pandas
import multiprocessing

directory = sys.argv[1]
output_path = os.path.join(directory, "features.pkl")
processes = int(multiprocessing.cpu_count()-1)*2
args = []
manager = multiprocessing.Manager()
Global = manager.Namespace()
Global.min_features = set()
Global.max_features = set()


def feature_compare(filename):
    basename = filename.split(".csv")[0]
    filepath = os.path.join(directory, filename)
    df = pandas.read_csv(
        filepath, index_col=0, header=None, skiprows=1).T
    cols = set(df.dropna(axis=1).columns)
    Global.max_features = cols.union(Global.max_features)
    if Global.min_features:
        cols = cols.intersection(Global.min_features)
    Global.min_features = cols

    print(f"After {basename} with {len(Global.min_features)} min features and {len(Global.max_features)} max features returned")

print(f"openning all csv files...")

for filename in os.listdir(directory):
    if filename.endswith(".csv") and not output_path.endswith(filename):
        args.append((filename,))

print(f"starting feature check with {len(args)} files...")


with multiprocessing.Pool(processes) as pool:
    pool.starmap(
        feature_compare,
        args
    )
features = {"min": Global.min_features, "max": Global.max_features}

print(f"process finished with {len(Global.min_features)} min features and {len(Global.max_features)} max features")
with open(output_path, "wb") as f:
    pickle.dump(features, f)

