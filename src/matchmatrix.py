import os
import sys
import itertools
from multiprocessing.dummy import Pool
import numpy
import joblib
from Bio.Seq import Seq
from typing import Iterable

fasta = sys.argv[1]
outpath = sys.argv[2]
threads = 30
Global = dict()

WINDOW=10

def build_matrix(seq1: Seq, seq2: Seq, max_window: int = WINDOW, min_window: int = 1):
    """
    Primeira sequência é a coluna e segunda é a linha.
    Retorna a soma de todas as janelas em 2 dimensões, na normal e na reversa.

    :rtype: numpy array (len segunda, len primeira, 2)
    """
    len2 = len(seq2)
    len1 = len(seq1)
    seq1 = str(seq1.upper())
    seq2_reverse = str(seq2.upper().reverse_complement())
    seq2_invert = str(seq2.upper())[::-1]
    seq2 = str(seq2.upper())
    rows = numpy.zeros((len2, len1, 3))

    for w in range(min_window, max_window+1):
        for b in range(0, len2-w+1):
            for r in range(0, len1-w+1):
                if seq1[r:r+w] == seq2[b:b+w]:
                    numpy.fill_diagonal(rows[b:b+w, r:r+w, 0], w)
                if seq1[r:r+w] == seq2_reverse[b:b+w]:
                    numpy.fill_diagonal(rows[:, ::-1][b:b+w, r:r+w, 1], w)
                if seq1[r:r+w] == seq2_invert[b:b+w]:
                    numpy.fill_diagonal(rows[:, ::-1][b:b+w, r:r+w, 2], w)

    return rows

def create_dict(ref_fasta: Iterable[str], fasta: Iterable[str]):
    ref_title, ref_sequence = ref_fasta
    Global[ref_title[1:-1]] = {}

def parallel_iterate(ref_fasta: Iterable[str], fasta: Iterable[str]):
    ref_title, ref_sequence = ref_fasta
    title, sequence = fasta
    matrix = build_matrix(ref_sequence[:-1], sequence[:-1])
    Global[ref_title[1:-1]][title[1:-1]] = matrix
    print(ref_title[1:-1], len(Global[ref_title[1:-1]]))

def main():
    print("Avoiding race condition creating dict structure first...")
    with open(fasta, "r") as sequences:
        with Pool(threads) as pool:
            pool.starmap(
                create_dict,
                itertools.product(
                    itertools.zip_longest(sequences, sequences), repeat=2))

    with open(fasta, "r") as sequences:
        print(f"starting to build matrix for {len(sequences.readlines())} sequences")

    with open(fasta, "r") as sequences:
        with Pool(threads) as pool:
            pool.starmap(
                parallel_iterate,
                itertools.product(
                    itertools.zip_longest(sequences, sequences), repeat=2))
    print("writing files...")
    joblib.dump(Global, os.path.join(outpath, "matrices.job"))
    print("done!")

if __name__ == "__main__":
    main()