import numpy
import tensorflow
from tensorflow import Tensor
from typing import Tuple, List
from src.variants.ssim import SSIMVariant


class SSIMMultiScaleVariant(SSIMVariant):

    name = "MultiScale Structural Similarity Index Measure"

    def _call_alg(self, image: Tensor, other: Tensor) -> float:

        return tensorflow.image.ssim_multiscale(
                                image,
                                other,
                                **self._alg_params).numpy()

    def _dynamic_find_image_match(
        self,
        small_img: numpy.ndarray,
        big_img: numpy.ndarray,
        s1: int, s2: int,
        step: int, window: int   
    ) -> float:
        max_score = 0
        for i in range(0, big_img.shape[1]-window, step):
            score = self._call_alg(
                tensorflow.expand_dims(small_img, axis=0),
                tensorflow.expand_dims(big_img[s1+i:s2+i, s1+i:s2+i], axis=0))[0]
            max_score = numpy.max((score, max_score))
        return max_score
 
    def _match_images(self, image: Tensor, other: Tensor) -> Tuple[float, List[int]]:
        if image.shape[1] > other.shape[1]:
            max_img = image.numpy()[0]
            min_img = other.numpy()[0]
        else:
            min_img = image.numpy()[0]
            max_img = other.numpy()[0]
        filter_size = min_img.shape[1]
        step = 1
        score = self._dynamic_find_image_match(
            min_img, max_img, 0, filter_size, step, filter_size)
        return score, []

