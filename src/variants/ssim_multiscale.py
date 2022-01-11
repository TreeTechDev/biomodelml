import numpy
import tensorflow
from tensorflow import Tensor
from src.variants.ssim import SSIMVariant



class SSIMMultiScaleVariant(SSIMVariant):

    name = "MultiScale Structural Similarity Index Measure"

    def _call_alg(self, image: Tensor, other: Tensor) -> numpy.ndarray:
        return tensorflow.image.ssim_multiscale(
                                image,
                                other,
                                max_val=255, filter_size=11,
                                filter_sigma=1.5, k1=0.01, k2=0.03)    
