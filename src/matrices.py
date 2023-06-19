import os
import numpy
from matplotlib import pyplot
from Bio.Seq import Seq
from biotite.sequence.align import SubstitutionMatrix


def _weight_ptns(seq1: Seq, seq2: Seq, max_window: int):
    subst_matrix = SubstitutionMatrix.dict_from_db("BLOSUM62")
    min_subst = min(subst_matrix.values())
    max_subst = max(subst_matrix.values()) + abs(min_subst)

    rows = numpy.zeros((len(seq2), len(seq1)), numpy.uint8)
    for line, letter1 in enumerate(seq2):
        for col, letter2 in enumerate(seq1):
            rows[line, col] = round((subst_matrix[letter1, letter2]+abs(min_subst))*max_window/max_subst)
    return rows

def _weight_seqs(seq1: Seq, seq2: Seq, rows: numpy.ndarray, max_window: int):
    indexes = dict()
    for line, letter in enumerate(seq2):
        if letter not in indexes:
            indexes[letter] = numpy.where(numpy.array(list(seq1)) == seq2[line])[0]
        idx = indexes[letter]
        rows[line, idx] = max_window
    return rows


def build_matrix(seq1: Seq, seq2: Seq, max_window: int, seq_type: str):
    """
    Primeira sequência é a coluna e segunda é a linha.
    Retorna a soma de todas as janelas em 2 dimensões, na normal e na reversa.

    :rtype: numpy array (len segunda, len primeira, 2)
    """
    if seq_type == "N":
        len2 = len(seq2)
        len1 = len(seq1)
        seq1 = str(seq1)
        seq2_complement = str(seq2.complement())
        seq2 = str(seq2)
        rows = numpy.zeros((len2, len1, 3), numpy.uint8)

        #  red
        rows[:, :, 0] = _weight_seqs(seq1, seq2, rows[:, :, 0], max_window)
        #  green
        rows[:, :, 1] = _weight_seqs(seq1, seq2_complement, rows[:, :, 1], max_window)
        #  blue    
        all_lines, all_columns = numpy.where((rows[:, :, 0] == 0) & (rows[:, :, 1] == 0))
        rows[all_lines, all_columns, 2] = max_window
    elif seq_type == "P":
        rows = _weight_ptns(seq1, seq2, max_window)
    return rows

def _save_grayscale(
        matrix: numpy.ndarray,
        output_path: str,
        filename: str
    ):
    os.makedirs(output_path, exist_ok=True)
    pyplot.imsave(
            os.path.join(output_path, filename), matrix, cmap=pyplot.cm.gray
        )

def _produce_grayscale_groups(
        matrix: numpy.ndarray,
        output_path: str,
        filename: str):

    groups = {
        "gray_max": numpy.max,
        "gray_mean": numpy.mean
    }

    for name, group in groups.items():
        _save_grayscale(group(matrix, axis=2), os.path.join(output_path, name), filename)


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

    _save_grayscale(matrix[:, :, channels[channel_name]], os.path.join(output_path, channel_name), filename)


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
    if len(kwargs["matrix"].shape) == 3:
        _produce_by_channel("full", **kwargs)
        _produce_by_channel("red", **kwargs)
        _produce_by_channel("green", **kwargs)
        _produce_by_channel("blue", **kwargs)
        _produce_by_channel("red_blue", **kwargs)
        _produce_by_channel("red_green", **kwargs)
        _produce_by_channel("green_blue", **kwargs)
        _produce_grayscale_by_channel("gray_r", **kwargs)
        _produce_grayscale_by_channel("gray_g", **kwargs)
        _produce_grayscale_by_channel("gray_b", **kwargs)
        _produce_grayscale_groups(**kwargs)
    else:
        kwargs["output_path"] = os.path.join(kwargs["output_path"], "full")
        _save_grayscale(**kwargs)


def save_image_by_matrices(
        name1: str, name2: str, seq1: Seq, seq2: Seq,
        max_window: int, output_path: str, seq_type: str):
    matrix = build_matrix(seq1, seq2, max_window, seq_type)
    filename = f"{name1}x{name2}.png" if name1 != name2 else f"{name1}.png"
    _produce_channel_images(
        matrix=matrix, output_path=output_path, filename=filename)
