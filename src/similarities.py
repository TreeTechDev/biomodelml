import sys
import pandas
import vaex
import multiprocessing
from scipy.spatial.distance import cosine
from scipy.stats import spearmanr, kendalltau

file_path = sys.argv[1]
chunksize = 30_000
processes = int(multiprocessing.cpu_count()-1)
manager = multiprocessing.Manager()
output_values = manager.list()


def similarity(df, column, compare):
    similarity = str(1 - cosine(df[column].values, df[compare].values))
    rho, p_rho = spearmanr(df[column].values, df[compare].values)
    tau, p_tau = kendalltau(df[column].values, df[compare].values)
    correlations = f"{similarity},{rho},{p_rho},{tau},{p_tau}"
    output = f"\n{column},{compare},{correlations}"
    if column != compare:
        output += f"\n{compare},{column},{correlations}"
    output_values.append(output)

df = vaex.open(file_path)

with open(f"{file_path}.similarity.csv", "w") as sim:
    sim.write("column,compare,similarity,rho,p_rho,tau,p_tau")

columns = list(df.columns)[1:]

args = [(df, column, compare) for i, column in enumerate(columns) for compare in columns[i:]]
with multiprocessing.Pool(processes) as pool:
    pool.starmap(
        similarity,
        args
    )

with open(f"{file_path}.similarity.csv", "a") as sim:
    sim.writelines(output_values)