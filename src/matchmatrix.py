import sys
import itertools
import multiprocessing
from typing import List, Iterable

fasta = sys.argv[1]
outpath = sys.argv[2]
processes = int(multiprocessing.cpu_count()/2)
output_csvs = {}

def build_matrix(seq1: str, seq2: str) -> List[str]:
    matrix = [0]*len(seq1)
    for idx, ref_base in enumerate(seq1):
        for base in seq2:
            if ref_base == base:
                matrix[idx]+=5
            else:
                matrix[idx]-=4
    return matrix

def get_file(title: str, sequence: str):
    if output_csvs.get(title):
        return output_csvs[title]
    else:
        output_csvs[title] = open(f"{outpath}/{title}.csv", "a")
        output_csvs[title].write("sequence,"+",".join(sequence[:-1])+"\n")
        return output_csvs[title]

def parallel_iterate(ref_fasta: Iterable[str], fasta: Iterable[str]):
    ref_title, ref_sequence = ref_fasta
    title, sequence = fasta
    matrix = build_matrix(ref_sequence[:-1], sequence[:-1])
    str_matrix = ",".join([str(m) for m in matrix])
    output_csv = get_file(ref_title[1:-1], ref_sequence)
    output_csv.write(f"{title[1:-1]},{str_matrix}\n")


with open(fasta, "r") as sequences:
    with multiprocessing.Pool(processes) as pool:
        pool.starmap(
            parallel_iterate,
            itertools.product(
                itertools.zip_longest(sequences, sequences), repeat=2))

for csv in output_csvs:
    csv.close()