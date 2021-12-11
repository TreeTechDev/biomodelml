import os
import sys
from typing import Optional
from src.experiment import Experiment
from src.variants.control import ControlVariant
from src.variants.uqi import UQIVariant
from src.variants.ssim_multiscale import SSIMMultiScaleVariant
from src.variants.deep_search.variant import DeepSearchVariant


def main(fasta_file: str, output_path: str, image_path: Optional[str] = None):
    Experiment(
        output_path,
        # ControlVariant(fasta_file),
        SSIMMultiScaleVariant(fasta_file, image_path)
        # UQIVariant(fasta_file, image_path),
        # DeepSearchVariant(fasta_file, image_path)
    ).run().save()

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    output_path = sys.argv[2]
    output_path = os.path.join(output_path, fasta_file.split(".")[0].split("/")[-1])
    os.makedirs(output_path, exist_ok=True)
    try:
        image_path = sys.argv[3]
    except IndexError:
        image_path = None
    main(fasta_file, output_path, image_path)