import tensorflow
from tensorflow import Tensor
from typing import Tuple, List
from src.structs import ImgDebug
from src.variants.ssim import SSIMVariant
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
                                **self._alg_params).numpy()

    def _dynamic_find_image_match(
        self,
        small_img: Tensor,
        big_img: Tensor,
        s1: int, s2: int,
        step: int, window: int   
    ) -> Tuple[float, List[ImgDebug]]:
        max_score = 0
        debugs = list()
        for i in range(0, big_img.shape[1]-window+1, step):
            score = self._call_alg(
                tensorflow.expand_dims(small_img, axis=0),
                tensorflow.expand_dims(big_img[s1+i:s2+i, s1+i:s2+i], axis=0))[0]
            max_score = max((score, max_score))
            debugs.append(
                ImgDebug(str(score), str(s1+i), str(s1+i), str(s2+i), str(s2+i), str(big_img.shape[1])))
        return max_score, debugs
 
    def _match_images(self, image: Tensor, other: Tensor) -> Tuple[float, List[ImgDebug]]:
        if image.shape[1] > other.shape[1]:
            max_img = image[0]
            min_img = other[0]
        else:
            min_img = image[0]
            max_img = other[0]
        filter_size = min_img.shape[1]
        step = 1
        return self._dynamic_find_image_match(
            min_img, max_img, 0, filter_size, step, filter_size)
