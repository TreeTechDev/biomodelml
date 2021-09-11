import sys
import os
import pandas
import pickle
import asyncio
import multiprocessing
from scipy.spatial.distance import jensenshannon

directory = sys.argv[1]
output_path = os.path.join(directory, "similarity.txt")
processes = int(multiprocessing.cpu_count()-1)
checkpoint = os.path.join(directory, "checkpoint.pkl")
with open(os.path.join(directory, "features.pkl"), "rb") as f:
    features = pickle.load(f)["max"]

print(f"loaded {len(features)} features")

async def write_checkpoint(compare_args):
    print("dumping checkpoint to pickle...")
    with open(checkpoint, "wb") as f:
        pickle.dump(compare_args, f)
    print("checkpoint did")

def similarity(df, filename, filename_compare):
    basename = filename.split(".csv")[0]
    basename_compare = filename_compare.split(".csv")[0]
    filepath = os.path.join(directory, filename_compare)
    df_compare = pandas.read_csv(
        filepath, index_col=0, header=None, skiprows=1).T
    df_join = pandas.concat(
        [df, df_compare], join="outer", ignore_index=True, copy=False).reindex(
            columns=features, fill_value=0.0, copy=False).fillna(0.0)
    try:
        similarity = 1 - jensenshannon(*df_join.values)
    except Exception as e:
        print(f"fail to compute similarity between {basename} and {basename_compare}")
        raise e
    output = f"\n{basename},{basename_compare},{similarity}"
    if basename != basename_compare:
       output += f"\n{basename_compare},{basename},{similarity}"
    with open(output_path, "a") as sim:
        sim.write(output)
        sim.flush()
    print(f"wrote {basename} and {basename_compare} similarity {similarity} to file")

async def main():
    if os.path.exists(checkpoint):
        print("getting files remaining from checkpoint...")
        with open(checkpoint, "rb") as f:
            args = pickle.load(f)
    else:
        print("openning all csv files...")
        args = []
        for filename in os.listdir(directory):
            if filename.endswith(".csv") and not output_path.endswith(filename):
                args.append((directory, filename))

    with open(output_path, "w") as sim:
        sim.write("gene,compare,similarity")

    print(f"starting similarity check with {len(args)} combinations...")

    # combinations_with_replacement 2 x 2
    for i, path in enumerate(args):
        filepath = os.path.join(*path)
        df = pandas.read_csv(filepath, index_col=0, header=None, skiprows=1).T
        compare_args = [(df, path[1], f) for d, f in args[i:]]
        print(f"comparing {i+1} from {len(args)} sequences...")
        asyncio.ensure_future(write_checkpoint(compare_args))
        with multiprocessing.Pool(processes) as pool:
            pool.starmap(
                similarity,
                compare_args
            )
    os.remove(checkpoint)

if __name__ == '__main__':
    loop = asyncio.get_event_loop()
    loop.run_until_complete(main())
    
    # Let's also finish all running tasks:
    pending = asyncio.Task.all_tasks()
    loop.run_until_complete(asyncio.gather(*pending))