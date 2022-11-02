import os
import numpy
from matplotlib import pyplot
from Bio.Seq import Seq


def _weight_seqs(seq1: Seq, seq2: Seq, rows: numpy.ndarray):
    indexes = dict()
    for line, letter in enumerate(seq2):
        if letter not in indexes:
            indexes[letter] = numpy.where(numpy.array(list(seq1)) == seq2[line])[0]
        idx = indexes[letter]
        rows[line, idx] = 255
    # all_lines, all_columns = numpy.where(rows == 1)
    # lines, columns = all_lines.copy(), all_columns.copy()
    # line = lines[0]
    # col = columns[0]
    # total_r = 1
    # while(lines.size):
    #     if numpy.where(all_columns[numpy.where(all_lines == line + total_r)[0]] == col + total_r)[0].size == 1:
    #         total_r += 1
    #         for k in range(total_r):
    #             rows[line+k, col+k] = max(rows[line+k, col+k], total_r)
    #     else:
    #         line = lines[0]
    #         lines = lines[1:]
    #         col = columns[0]
    #         columns = columns[1:]
    #         total_r = 1
    # lines, columns = all_lines.copy(), all_columns.copy()
    # line = lines[0]
    # col = columns[0]
    # total_l = 1
    # while(lines.size):
    #     if numpy.where(all_columns[numpy.where(all_lines == line + total_l)[0]] == col - total_l)[0].size == 1:
    #         total_l += 1
    #         for k in range(total_l):
    #             rows[line+k, col-k] = max(rows[line+k, col-k], total_l)

    #     else:
    #         line = lines[0]
    #         lines = lines[1:]
    #         col = columns[0]
    #         columns = columns[1:]
    #         total_l = 1
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
    seq2_reverse = str(seq2.reverse_complement())
    seq2 = str(seq2)
    rows = numpy.zeros((len2, len1, 3))
    norm = numpy.zeros((1, 1, 3))

    #  red
    rows[:, :, 0] = _weight_seqs(seq1, seq2, rows[:, :, 0])
    # norm[:, :, 0] = max(rows[:, :, 0].max(), numpy.e)
    #  green
    rows[:, :, 1] = _weight_seqs(seq1, seq2_reverse, rows[:, :, 1])
    # norm[:, :, 1] = max(rows[:, :, 1].max(), numpy.e)
    #  blue    
    normalizer = numpy.zeros((len2, len1))
    all_lines, all_columns = numpy.where(rows[:, :, 0] == 0)
    for i in range(all_lines.size):
        line = all_lines[i]
        col = all_columns[i]
        diags = ((line-1, col-1), (line-1, col+1), (line+1, col+1), (line+1, col-1))
        parents = ((line-2, col-2), (line-2, col+2), (line+2, col+2), (line+2, col-2))
        args = {}
        for j, (x, y) in enumerate(diags):
            if x < rows.shape[0] and y < rows.shape[1]:
                val = rows[x, y, 2]
                px, py = parents[j]
                if px < rows.shape[0] and py < rows.shape[1] and val > rows[px, py, 2]:
                    continue
                args[j] = val                

        if not args or (args.get(0) or args.get(1)) and (args.get(2) or args.get(3)): continue
        normalizer[line, col] = 255 #max(list(args.values()))
    rows[:, :, 2] = normalizer
    # norm[:, :, 2] = max(rows[:, :, 2].max(), numpy.e)
    #  norm
    # rows = numpy.ma.log(rows).filled(0) * max_window / numpy.ma.log(norm).filled(1)
    # rows = rows * max_window / norm
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
    _produce_by_channel("full", **kwargs)
    _produce_by_channel("red_blue", **kwargs)
    _produce_by_channel("red_green", **kwargs)
    _produce_by_channel("green_blue", **kwargs)
    _produce_grayscale_by_channel("gray_r", **kwargs)
    _produce_grayscale_by_channel("gray_g", **kwargs)
    _produce_grayscale_by_channel("gray_b", **kwargs)
    _produce_grayscale_groups(**kwargs)


def save_image_by_matrices(
        name1: str, name2: str, seq1: Seq, seq2: Seq,
        max_window: int, output_path: str):
    matrix = build_matrix(seq1, seq2, max_window)
    filename = f"{name1}x{name2}.png" if name1 != name2 else f"{name1}.png"
    color_matrix = matrix.astype(numpy.uint8)
    _produce_channel_images(
        matrix=color_matrix, output_path=output_path, filename=filename)
    _produce_results_images(
        color_matrix, output_path, filename, max_window)
