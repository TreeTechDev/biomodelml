import os
import pandas
import numpy
import tensorflow
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import cpu_count
from typing import List, Tuple
from tensorflow import Tensor
from src.context import RecursionContext
from src.variants.variant import Variant
from src.structs import DistanceStruct, ImgMap


MAX_POSSIBLE_SCORE = 1.0
DEFAULT_PARAMS = dict(
    filter_sigma=1.5,
    k1=0.01,
    k2=0.03,
    filter_size=11,
    max_val=255
)


class SSIMVariant(Variant):

    name = "Structural Similarity Index Measure"

    def __init__(self, fasta_file: str, sequence_type: str, image_folder: str, **alg_params):
        super().__init__(fasta_file, sequence_type)
        self._image_folder = image_folder
        self._alg_params = DEFAULT_PARAMS
        self._alg_params.update(alg_params)
        self.filter_size = self._alg_params["filter_size"]
        self._executor = ThreadPoolExecutor(max_workers=cpu_count())
    
    def _call_alg(self, image: Tensor, other: Tensor) -> float:
        return tensorflow.image.ssim(
                                image,
                                other,
                                **self._alg_params)[0].numpy()

    def _read_image(self, img_name: str) -> Tensor:
        img = tensorflow.expand_dims(
            tensorflow.image.decode_image(
                tensorflow.io.read_file(
                    os.path.join(
                        self._image_folder, img_name)), channels=3), axis=0)
        if tensorflow.math.equal(img[:,:,:,2], img[:,:,:,1], img[:,:,:,0]).numpy().all():
            img = tensorflow.image.rgb_to_grayscale(img)
        return img

    def _greedy_find_image_match(
        self,
        small_img: numpy.ndarray,
        big_img: numpy.ndarray,
        c11: int, c12: int,
        l21: int, l22: int,
        c21: int, c22: int,
        max_score: int, last_line: int, step: int
    ):
        if (c22 > big_img.shape[1]) or (l22 > big_img.shape[0]):
            return max_score, last_line
        
        score = self._call_alg(
            tensorflow.expand_dims(small_img[:, c11:c12], axis=0),
            tensorflow.expand_dims(big_img[l21:l22, c21:c22], axis=0))
        if score == MAX_POSSIBLE_SCORE:
            return score, l21
        elif score > max_score:
            return self._greedy_find_image_match(
                small_img, big_img, c11, c12, l21+step, l22+step, c21+step, c22+step, score, l21, step)
        else:
            return self._greedy_find_image_match(
                small_img, big_img, c11, c12, l21+step, l22+step, c21+step, c22+step, max_score, last_line, step)


    def _find_best_col(self, min_img, max_img, mask_size, filter_size, step) -> ImgMap:
        last_line = 0
        scores = list()
        positions = list()
        start_score = 0
        for j in range(0, mask_size - filter_size, step):
            window = min(mask_size, j+filter_size)
            last_score, last_line = self._greedy_find_image_match(
                min_img, max_img,
                j, window,
                last_line, mask_size+last_line,
                last_line+j, last_line+window,
                start_score, last_line, step
            )
            scores.append(last_score)
            positions.append(
                (str(last_score), str(last_line+j), str(last_line),
                str(last_line+window), str(mask_size+last_line), str(max_img.shape[1])))
        return ImgMap(positions=positions, scores=scores)


    def _match_images(self, image: Tensor, other: Tensor) -> Tuple[float, List[int]]:
        if image.shape[1] > other.shape[1]:
            max_img = image.numpy()[0]
            min_img = other.numpy()[0]
        else:
            min_img = image.numpy()[0]
            max_img = other.numpy()[0]

        mask_size = min_img.shape[1]
        filter_size = self.filter_size
        step = filter_size // 2
        default_args = (mask_size, filter_size, step)
        img_map = self._find_best_col(min_img, max_img, *default_args)
        return numpy.median(img_map.scores), img_map.positions
    
    def _compare(self, img1: str, img2: str) -> Tuple[float, List[int]]:
        img = self._read_image(img1)
        other = self._read_image(img2)
        with RecursionContext():
            results, p = self._match_images(img, other)
            positions = [img1, img2, p] if p else []
        return results, positions     

    def build_matrix(self) -> DistanceStruct:
        files = os.listdir(self._image_folder)
        indexes = {".".join(img.split('.')[:-1]): img.split('.')[-1] for img in files}
        diff = set(self._names).difference(set(indexes.keys()))
        if diff:
            raise IOError(f"Sequences without image created: {diff}")
        files = []
        for i in self._names:
            if indexes.get(i):
                files.append(f"{i}.{indexes.get(i)}")
        indexes = self._names
        df = pandas.DataFrame(index=indexes, columns=indexes)
        last_ids = list()
        positions = list()
        for idx, img1 in enumerate(files):
            results = list()
            idx1 = indexes[idx]
            futures = [self._executor.submit(self._compare, img1, img2) for img2 in files[idx:]]
            for future in futures:
                r, p = future.result()
                results.append(r)
                positions += p
            if last_ids:
                df.loc[idx1, indexes[idx:]] = results
                df.loc[idx1, last_ids] = df.loc[last_ids, idx1]
            else:
                df.loc[idx1, :] = results
            last_ids.append(idx1)

        return DistanceStruct(
            names=indexes,
            matrix=1.0-df.to_numpy(numpy.float64),
            img_positions=positions)
