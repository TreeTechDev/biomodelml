import tensorflow
from tensorflow import Tensor
from src.variants.ssim import SSIMVariant



class SSIMMultiScaleVariant(SSIMVariant):

    name = "MultiScale Structural Similarity Index Measure"

    def _call_alg(self, image: Tensor, other: Tensor) -> float:
        return tensorflow.image.ssim_multiscale(
                                image,
                                other,
                                **self._alg_params)[0].numpy()
