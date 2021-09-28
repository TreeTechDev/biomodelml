import modin.pandas as pandas
import sys
import os
import multiprocessing

processes = int(multiprocessing.cpu_count()-1)
directory = sys.argv[1]
output = sys.argv[2]
args, out_args = [], []

for filename in os.listdir(directory):
    if filename.endswith(".csv"):
        args.append((directory, filename))
        out_args.append((output, filename))

def create_file(directory, filename):
    filepath = os.path.join(directory, filename)
    with open(filepath, "w") as csv:
        csv.write(f"cluster,prob\n")
        csv.flush()

def convert(directory, df, basename):
    filepath = os.path.join(directory, f"{basename}.csv")
    if not os.path.exists(filepath):
        print(f"file {filepath} do not exists")
        return 
    with open(filepath, "a") as f:
        f.write("\n".join(
            (df["index"].str.replace("NC_000913.3:", "") + "," + df[basename].astype(str)).values
            ) + "\n")
        f.flush()

print(f"Creating {len(out_args)} vector files...")

with multiprocessing.Pool(processes) as pool:
    pool.starmap(
        create_file,
        out_args
    )

print(f"Building {len(args)} col args to vector convertion...")

for i, path in enumerate(args):
    filepath = os.path.join(*path)
    df = pandas.read_csv(filepath)
    columns = list(df.columns)[1:]

    col_args = []
    for column in columns:
        col_args.append((output, df, column))

    print(f"converting {len(col_args)} columns {i}/{len(args)} files to vector for file {path[1]}...")
    with multiprocessing.Pool(processes) as pool:
        pool.starmap(
            convert,
            col_args
        )

