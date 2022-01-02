import os
import sys
from multiprocessing import Pool, cpu_count
from matplotlib import pyplot
import numpy
from Bio.Seq import Seq
from Bio import SeqIO


def build_matrix(seq1: Seq, seq2: Seq, max_window: int, min_window: int):
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


def save_image_by_matrices(
        name1: str, name2: str, seq1: Seq, seq2: Seq,
        max_window: int, min_window: int, output_path: str):
    matrix = build_matrix(seq1, seq2, max_window, min_window)
    max_rgb = 255
    filename = f"{name1}x{name2}.png" if name1 != name2 else f"{name1}.png"
    pyplot.imsave(
        os.path.join(output_path, filename),
        (matrix*max_rgb/max_window).astype(numpy.uint8)
    )


def main(fasta_file: str, output_path: str):

    max_window = 20
    min_window = 1
    procs = cpu_count()

    with open(fasta_file, "r") as handle:
        sequences = SeqIO.parse(handle, "fasta")
        print(f"File with {len(list(sequences))} sequences")

    with open(fasta_file) as handle:
        sequences = SeqIO.parse(handle, "fasta")
        to_run = []
        for s in sequences:
            if not os.path.exists(os.path.join(output_path, f"{s.description}.png")):
                to_run.append(
                    (s.description, s.description, s.seq, s.seq, max_window,
                  min_window, output_path)
                )
        print(f"starting to build image matrix for {len(to_run)} sequences")
            
        with Pool(procs) as pool:
            pool.starmap(
                save_image_by_matrices,
                to_run
            )


if __name__ == "__main__":
    fasta_file = sys.argv[1]
    outpath = sys.argv[2]
    outpath = os.path.join(outpath, fasta_file.split(".")[0].split("/")[-1])
    os.makedirs(outpath, exist_ok=True)
    main(fasta_file, outpath)
