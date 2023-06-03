import os
import sys
import numpy
from multiprocessing import Pool, cpu_count
from matplotlib import pyplot
from Bio.Seq import Seq
from Bio import SeqIO


def _weight_seqs(seq1: Seq, seq2: Seq, rows: numpy.ndarray, max_window: int):
    indexes = dict()
    for line, letter in enumerate(seq2):
        if letter == "N":
            continue
        elif letter not in indexes:
            indexes[letter] = numpy.where(numpy.array(list(seq1)) == seq2[line])[0]
        idx = indexes[letter]
        rows[line, idx] = max_window
    return rows


def build_matrix(seq1: Seq, seq2: Seq, max_window: int):
    """
    Primeira sequência é a coluna e segunda é a linha.
    Retorna a soma de todas as janelas em 2 dimensões, na normal e na reversa.

    :rtype: numpy array (len segunda, len primeira, 2)
    """
    len2 = len(seq2)
    len1 = len(seq1)
    seq1 = str(seq1)
    seq2_complement = str(seq2.complement())
    seq2 = str(seq2)
    rows = numpy.zeros((len2, len1, 3), numpy.int8)

    #  red
    rows[:, :, 0] = _weight_seqs(seq1, seq2, rows[:, :, 0], max_window)
    #  green
    rows[:, :, 1] = _weight_seqs(seq1, seq2_complement, rows[:, :, 1], max_window)
    #  blue    
    all_lines, all_columns = numpy.where((rows[:, :, 0] == 0) & (rows[:, :, 1] == 0))
    rows[all_lines, all_columns, 2] = max_window
    return rows


def _produce_results_images(
    matrix: numpy.ndarray, output_path: str,
    filename: str, max_window: int
):
    new_name = f"{max_window}_{filename}"
    matrix = numpy.invert(matrix)

    os.makedirs(os.path.join(output_path, "result"), exist_ok=True)
    pyplot.imsave(
        os.path.join(output_path, "result", f"red_{new_name}"),
        matrix[:, :, 0], cmap=pyplot.cm.gray
    )
    pyplot.imsave(
        os.path.join(output_path, "result", f"green_{new_name}"),
        matrix[:, :, 1], cmap=pyplot.cm.gray
    )
    pyplot.imsave(
        os.path.join(output_path, "result", f"blue_{new_name}"),
        matrix[:, :, 2], cmap=pyplot.cm.gray
    )
    # matrix = numpy.invert(matrix)

    # black_pixels = numpy.where(
    #     (matrix[:, :, 0] < 25) & 
    #     (matrix[:, :, 1] < 25) & 
    #     (matrix[:, :, 2] < 25)
    # )
    # matrix[black_pixels] = [255, 255, 255]

    # pyplot.imsave(
    #     os.path.join(output_path, f"white_{new_name}"),
    #     matrix
    # )

def _produce_grayscale_groups(
    matrix: numpy.ndarray,
    output_path: str,
    filename: str):

    groups = {
        "gray_max": numpy.max,
        "gray_mean": numpy.mean
    }

    for name, group in groups.items():
        os.makedirs(os.path.join(output_path, name), exist_ok=True)
        pyplot.imsave(
            os.path.join(output_path, name, filename),
            group(matrix, axis=2), cmap=pyplot.cm.gray
        )

def _produce_grayscale_by_channel(
    channel_name: str,
    matrix: numpy.ndarray,
    output_path: str,
    filename: str):

    channels = {
        "gray_r": 0,
        "gray_g": 1,
        "gray_b": 2
    }


    os.makedirs(os.path.join(output_path, channel_name), exist_ok=True)
    pyplot.imsave(
        os.path.join(output_path, channel_name, filename),
        matrix[:, :, channels[channel_name]], cmap=pyplot.cm.gray
    )


def _produce_by_channel(
    channel_name: str,
    matrix: numpy.ndarray,
    output_path: str,
    filename: str
):
    channels_not = {
        "red": (1, 2),
        "green": (0, 2),
        "blue": (0, 1),
        "red_blue": (1,),
        "red_green": (2,),
        "green_blue": (0,)
    }
    new_matrix = matrix.copy()
    os.makedirs(os.path.join(output_path, channel_name), exist_ok=True)
    for channel in channels_not.get(channel_name, []):
        new_matrix[:, :, channel] = 0
    pyplot.imsave(
        os.path.join(output_path, channel_name, filename),
        new_matrix
    )

def _produce_channel_images(**kwargs):
    _produce_by_channel("red", **kwargs)
    _produce_by_channel("green", **kwargs)
    _produce_by_channel("blue", **kwargs)
    _produce_by_channel("red_blue", **kwargs)
    _produce_by_channel("red_green", **kwargs)
    _produce_by_channel("green_blue", **kwargs)
    _produce_by_channel("full", **kwargs)
    _produce_grayscale_by_channel("gray_r", **kwargs)
    _produce_grayscale_by_channel("gray_g", **kwargs)
    _produce_grayscale_by_channel("gray_b", **kwargs)
    _produce_grayscale_groups(**kwargs)

def save_image_by_matrices(
        name1: str, name2: str, seq1: Seq, seq2: Seq,
        max_window: int, output_path: str):
    matrix = build_matrix(seq1, seq2, max_window)
    filename = f"{name1}x{name2}.png" if name1 != name2 else f"{name1}.png"
    _produce_channel_images(
        matrix=matrix, output_path=output_path, filename=filename)
    _produce_results_images(
        matrix, output_path, filename, max_window)


def main(fasta_file: str, output_path: str):

    max_window = 255
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
                    (s.description, s.description, s.seq, s.seq, max_window, output_path)
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
