#!/usr/bin/env python
# coding: utf-8

import os
import sys
from typing import Optional
from pathlib import Path
from src.experiment import Experiment
from src.variants.control import ControlVariant
from src.variants.sw import SmithWatermanVariant
from src.variants.nw import NeedlemanWunschVariant
from src.variants.uqi import UQIVariant
from src.variants.resized_ssim import ResizedSSIMVariant
from src.variants.resized_ssim_multiscale import ResizedSSIMMultiScaleVariant
from src.variants.windowed_ssim_multiscale import WindowedSSIMMultiScaleVariant
from src.variants.greedy_ssim import GreedySSIMVariant
from src.variants.unrestricted_ssim import UnrestrictedSSIMVariant
from src.variants.deep_search.variant import DeepSearchVariant


def main(fasta_file: str, output_path: str, sequence_type: str, image_path: Optional[str] = None):
    Experiment(
        Path(output_path),
        # ControlVariant(fasta_file, sequence_type),
        # SmithWatermanVariant(fasta_file, sequence_type),
        # NeedlemanWunschVariant(fasta_file, sequence_type),
        # ResizedSSIMVariant(fasta_file, sequence_type, image_path),
        # ResizedSSIMMultiScaleVariant(fasta_file, sequence_type, image_path),
        # WindowedSSIMMultiScaleVariant(fasta_file, sequence_type, image_path),
        # GreedySSIMVariant(fasta_file, sequence_type, image_path),
        # UnrestrictedSSIMVariant(fasta_file, sequence_type, image_path),
        # UQIVariant(fasta_file, sequence_type, image_path),
        DeepSearchVariant(fasta_file, sequence_type, image_path)
    ).run_and_save()

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    output_path = sys.argv[2]
    sequence_type = sys.argv[3] if sys.argv[3] in ["P", "N"] else None
    output_path = os.path.join(output_path, fasta_file.split(".")[0].split("/")[-1])
    try:
        image_path = sys.argv[4]
        if not os.path.exists(image_path): raise IOError("Image path not exists")
    except IndexError:
        image_path = None
    os.makedirs(output_path, exist_ok=True)
    main(fasta_file, output_path, sequence_type, image_path)
