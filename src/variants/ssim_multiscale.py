import numpy
import tensorflow
from tensorflow import Tensor
from src.variants.ssim import SSIMVariant

_MSSSIM_WEIGHTS = (0.0448, 0.2856, 0.3001, 0.2363, 0.1333)  # original paper
RESIZE_FACTOR = 2**(len(_MSSSIM_WEIGHTS)-1)  # https://github.com/tensorflow/tensorflow/issues/33840#issuecomment-633715778


class SSIMMultiScaleVariant(SSIMVariant):

    name = "MultiScale Structural Similarity Index Measure"

    def _call_alg(self, image: Tensor, other: Tensor) -> float:
        filter_resized = self.filter_size / RESIZE_FACTOR
        filter_size = min(filter_resized, max(
            1, min(image.shape[2]//RESIZE_FACTOR, other.shape[2]//RESIZE_FACTOR)
        ))
        self._alg_params["filter_size"] = filter_size
        # if filter_size != self.filter_size:
        #     print(f"Filter size redefined from {self.filter_resized} to {filter_size}")
        return tensorflow.image.ssim_multiscale(
                                image,
                                other,
                                power_factors=_MSSSIM_WEIGHTS,
                                **self._alg_params)[0].numpy()
