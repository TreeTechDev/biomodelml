import tensorflow
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import cpu_count
from tensorflow import Tensor
from src.variants.ssim_base import SSIMVariant, MAX_POSSIBLE_SCORE
from src.structs import ImgMap, ImgDebug


class UnrestrictedSSIMVariant(SSIMVariant):

    name = "Unrestricted Sliced Structural Similarity Index Measure"

    def __init__(self, fasta_file: str = None, sequence_type: str = None, image_folder: str = "", **alg_params):
        super().__init__(fasta_file, sequence_type, image_folder, **alg_params)
        self._executor = ThreadPoolExecutor(max_workers=cpu_count()*10)
    

    def _greedy_find_image_match(
        self,
        small_img: Tensor,
        big_img: Tensor,
        l1: int, l2: int,
        c1: int, c2: int,
        max_score: int, last_line: int, step: int
    ):
        c2 = min(c2, big_img.shape[1])
        l2 = min(l2, big_img.shape[0])

        if small_img.shape[0] > big_img[l1:l2, :].shape[0]:
            return max_score, last_line

        score = self._call_alg(
            tensorflow.expand_dims(small_img, axis=0),
            tensorflow.expand_dims(big_img[l1:l2, c1:c2], axis=0))
        if score == MAX_POSSIBLE_SCORE:
            return score, l1
        elif score > max_score:
            max_score = score
            last_line = l1
        if (c2 == big_img.shape[1]) or (l2 == big_img.shape[0]):
            return max_score, last_line

        return self._greedy_find_image_match(
            small_img, big_img, l1+step, l2+step, c1+step, c2+step, max_score, last_line, step)


    def _find_best_col(self, min_img: Tensor, max_img: Tensor, mask_size: int, filter_size: int, step: int) -> ImgMap:
        last_line = 0
        scores = list()
        debugs = list()
        start_score = 0

        for j in range(0, mask_size - filter_size + 1, step):
            min_col_limit = min(j+filter_size, min_img.shape[1])
            last_score, last_line = self._greedy_find_image_match(
                min_img[:, j:min_col_limit], max_img,
                0, mask_size,
                j, j+filter_size,
                start_score, last_line, step
            )
            scores.append(last_score)
            debugs.append(
                ImgDebug(str(last_score), str(last_line+j), str(last_line),
                str(last_line+filter_size+j), str(mask_size+last_line), str(max_img.shape[1])))
        return ImgMap(debugs=debugs, scores=scores)