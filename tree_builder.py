import os
import sys
from typing import Optional
from src.experiment import Experiment
from src.variants.control import ControlVariant
from src.variants.sw import SmithWatermanVariant
from src.variants.nw import NeedlemanWunschVariant
from src.variants.uqi import UQIVariant
from src.variants.ssim import SSIMVariant
from src.variants.ssim_multiscale import SSIMMultiScaleVariant
from src.variants.deep_search.variant import DeepSearchVariant


def main(fasta_file: str, output_path: str, sequence_type: str, image_path: Optional[str] = None):
    Experiment(
        output_path,
        ControlVariant(fasta_file, sequence_type),
        SmithWatermanVariant(fasta_file, sequence_type),
        NeedlemanWunschVariant(fasta_file, sequence_type),
        SSIMVariant(fasta_file, sequence_type, image_path),
        SSIMMultiScaleVariant(fasta_file, sequence_type, image_path),
        UQIVariant(fasta_file, sequence_type, image_path),
        DeepSearchVariant(fasta_file, sequence_type, image_path)
    ).run().save()

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    output_path = sys.argv[2]
    sequence_type = sys.argv[3] if sys.argv[3] in ["P", "N"] else None
    output_path = os.path.join(output_path, fasta_file.split(".")[0].split("/")[-1])
    os.makedirs(output_path, exist_ok=True)
    try:
        image_path = sys.argv[4]
    except IndexError:
        image_path = None
    main(fasta_file, output_path, sequence_type, image_path)
