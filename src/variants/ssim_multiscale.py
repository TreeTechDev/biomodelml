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
        l21: int, l22: int,
        c21: int, c22: int,
        step: int, window: int   
    ) -> float:
        s1_same_score = False
        c21 = min(c21, big_img.shape[1]-window)
        c22 = min(c22, big_img.shape[1])
        l21 = min(l21, big_img.shape[0]-small_img.shape[0])
        l22 = min(l22, big_img.shape[0])

        if (c22 == big_img.shape[1]) or (l22 == big_img.shape[0]):
            s1_same_score = True            
        else:
            s1 = self._dynamic_find_image_match(
                    small_img, big_img, l21+step, l22+step, c21+step, c22+step, step, window)

        score = self._call_alg(
            tensorflow.expand_dims(small_img, axis=0),
            tensorflow.expand_dims(big_img[l21:l22, c21:c22], axis=0))[0]

        if s1_same_score:
            s1 = score
        return numpy.max((s1, score))
 
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
            min_img, max_img, 0, filter_size, 0, filter_size, step, filter_size)
        return score, []

