import sys
import itertools
import multiprocessing
import numpy
from typing import List, Iterable, TextIO

fasta = sys.argv[1]
outpath = sys.argv[2]
processes = int(multiprocessing.cpu_count()-1)
manager = multiprocessing.Manager()
Global = manager.dict()

WINDOW=10

def build_matrix(seq1: str, seq2: str):
    """
    Primeira sequência é a coluna e segunda é a linha.
    Retorna a soma de todas as janelas em 2 dimensões, na normal e na reversa.

    :rtype: numpy array (len segunda, len primeira, 2)
    """
    len2 = len(seq2)
    len1 = len(seq1)
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    window = WINDOW
    rows = numpy.zeros((len2, len1, window, 2))

    for w in range(window):
        for b in range(0, len2-w):
            for r in range(0, len1-w):
                if seq1[r:r+w+1] == seq2[b:b+w+1]:
                    for i in range(w+1):
                        rows[b+i, r+i, w, 0] = 1
                if seq1[r:r+w+1] == seq2[::-1][b:b+w+1]:
                    for i in range(w+1):
                        rows[len2-b-i-1, r+i, w, 1] = 1
    return numpy.sum(rows, 2)


def parallel_iterate(ref_fasta: Iterable[str], fasta: Iterable[str]):
    ref_title, ref_sequence = ref_fasta
    title, sequence = fasta
    matrix = build_matrix(ref_sequence[:-1], sequence[:-1])
    Global.setdefault(ref_title, {})[title] = matrix

with open(fasta, "r") as sequences:
    print(f"starting to build files for {len(sequences.readlines())} sequences")

print("starting to build matrix for sequences")
with open(fasta, "r") as sequences:
    with multiprocessing.Pool(processes) as pool:
        pool.starmap(
            parallel_iterate,
            itertools.product(
                itertools.zip_longest(sequences, sequences), repeat=2))
