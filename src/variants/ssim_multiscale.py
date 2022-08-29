import numpy
import tensorflow
from tensorflow import Tensor
from src.variants.ssim import SSIMVariant

_MSSSIM_WEIGHTS = (0.0448, 0.2856, 0.3001, 0.2363, 0.1333)  # original paper
RESIZE_FACTOR = 2**(len(_MSSSIM_WEIGHTS)-1)  # https://github.com/tensorflow/tensorflow/issues/33840#issuecomment-633715778


class SSIMMultiScaleVariant(SSIMVariant):

    name = "MultiScale Structural Similarity Index Measure"
    filter_resized = 11
    filter_size = filter_resized * RESIZE_FACTOR

    def _call_alg(self, image: Tensor, other: Tensor) -> float:
        filter_size = min(self.filter_resized, max(
            1, min(image.shape[2]//RESIZE_FACTOR, other.shape[2]//RESIZE_FACTOR)
        ))
        # if filter_size != self.filter_size:
        #     print(f"Filter size redefined from {self.filter_resized} to {filter_size}")
        return tensorflow.image.ssim_multiscale(
                                image,
                                other,
                                max_val=255, filter_size=filter_size,
                                power_factors=_MSSSIM_WEIGHTS,
                                filter_sigma=1.5, k1=0.01, k2=0.03)[0].numpy()
