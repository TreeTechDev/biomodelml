import tensorflow
from tensorflow import Tensor
from src.variants.ssim_base import SSIMVariant
from concurrent.futures import ThreadPoolExecutor


class SSIMMultiScaleVariant(SSIMVariant):

    name = "MultiScale Structural Similarity Index Measure"

    def __init__(self, fasta_file: str, sequence_type: str, image_folder: str, **alg_params):
        super().__init__(fasta_file, sequence_type, image_folder, **alg_params)
        self._executor = ThreadPoolExecutor(max_workers=1)

    def _call_alg(self, image: Tensor, other: Tensor) -> float:

        return tensorflow.image.ssim_multiscale(
            image,
            other,
            **self._alg_params)[0].numpy()
