import sys
import itertools
import multiprocessing
from numpy import linalg
from typing import List, Iterable, TextIO

fasta = sys.argv[1]
outpath = sys.argv[2]
processes = int(multiprocessing.cpu_count()-1)
opened_csvs = dict()

def build_matrix(seq1: str, seq2: str) -> List[str]:
    matrix = []
    for ref_base in seq1:
        row = []
        for base in seq2:
            if ref_base == base:
                row.append(5)
            else:
                row.append(-4)
        matrix.append(row)
    
    return linalg.eigvalsh(matrix)

def create_file(title: str, sequence: str):
    with open(f"{outpath}/{title[1:-1]}.csv", "a") as csv:
        csv.write("sequence,"+",".join(sequence[:-1])+"\n")
        csv.flush()

def get_file(title: str) -> TextIO:
    if not opened_csvs.get(title[1:-1]):
        opened_csvs[title[1:-1]] = open(f"{outpath}/{title[1:-1]}.csv", "a")
    
    return opened_csvs[title[1:-1]]

def parallel_iterate(ref_fasta: Iterable[str], fasta: Iterable[str]):
    ref_title, ref_sequence = ref_fasta
    title, sequence = fasta
    matrix = build_matrix(ref_sequence[:-1], sequence[:-1])
    str_matrix = ",".join([str(m) for m in matrix])
    opened_csvs = get_file(ref_title)
    opened_csvs.write(f"{title[1:-1]},{str_matrix}\n")
    opened_csvs.flush()

with open(fasta, "r") as sequences:
    print(f"starting to build files for {len(sequences.readlines())} sequences")
with open(fasta, "r") as sequences:
    with multiprocessing.Pool(processes) as pool:
        pool.starmap(
            create_file,
            itertools.zip_longest(sequences, sequences)
        )

print("starting to build matrix for sequences")
with open(fasta, "r") as sequences:
    with multiprocessing.Pool(processes) as pool:
        pool.starmap(
            parallel_iterate,
            itertools.product(
                itertools.zip_longest(sequences, sequences), repeat=2))

for csv in opened_csvs:
    csv.close()