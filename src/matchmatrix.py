import sys
import itertools
import multiprocessing
from typing import List, Iterable

fasta = sys.argv[1]
processes = int(multiprocessing.cpu_count()/2)
output_csvs = [open(f"{fasta}.{p}.csv", "a") for p in range(processes)]
iter_csvs = itertools.cycle(output_csvs)

def build_matrix(seq1: str, seq2: str) -> List[str]:
    matrix = []
    for ref_base in seq1:
        for base in seq2:
            if ref_base == base:
                matrix.append("1")
            else:
                matrix.append("0")
    return matrix

def parallel_iterate(ref_fasta: Iterable[str], fasta: Iterable[str]):
    ref_title, ref_sequence = ref_fasta
    title, sequence = fasta
    matrix = build_matrix(ref_sequence, sequence)
    str_matrix = ",".join(matrix)
    next(iter_csvs).write(f"{ref_title},{title},{str_matrix}\n")


with open(fasta, "r") as sequences:
    with multiprocessing.Pool(processes) as pool:
        pool.starmap(
            parallel_iterate,
            itertools.product(
                itertools.zip_longest(sequences, sequences), repeat=2))

for csv in output_csvs:
    csv.close()