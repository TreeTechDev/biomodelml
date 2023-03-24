import os
import pandas
import numpy
import tensorflow
from tensorflow import Tensor
from typing import Tuple
from src.variants.variant import Variant
from src.structs import DistanceStruct


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

    def _upscale_images(self, image: Tensor, other: Tensor) -> Tuple[Tensor]:
        max_x = image.shape[1] if image.shape[1] > other.shape[1] else other.shape[1]
        max_y = image.shape[2] if image.shape[2] > other.shape[2] else other.shape[2]
        return (
            tensorflow.image.resize(image, (max_x, max_y), tensorflow.image.ResizeMethod.BICUBIC),
            tensorflow.image.resize(other, (max_x, max_y), tensorflow.image.ResizeMethod.BICUBIC)
        )

    # def _match_images(self, image: Tensor, other: Tensor) -> Tuple[Tensor]:
    #     if image.shape[1] > other.shape[1]:
    #         max_img = image.numpy().squeeze(0)
    #         min_img = other.numpy().squeeze(0)
    #     else:
    #         min_img = image.numpy().squeeze(0)
    #         max_img = other.numpy().squeeze(0)
        
    #     w, h, _ = min_img.shape
        
    #     res = cv2.matchTemplate(max_img, min_img, cv2.TM_CCOEFF_NORMED)
    #     _, _, _, max_loc = cv2.minMaxLoc(res)
    #     max_cropped = max_img[max_loc[1]:max_loc[1]+w, max_loc[0]:max_loc[0]+h, :]
    #     return (
    #         tensorflow.expand_dims(min_img, axis=0),
    #         tensorflow.expand_dims(max_cropped, axis=0)
    #     )

    def calc_alg(self, img_name1: str, img_name2: str) -> float:
        img, other = self._upscale_images(
                            self._read_image(img_name1),
                            self._read_image(img_name2)
                        )
        return self._call_alg(img, other)       

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
        for idx, img1 in enumerate(files):
            idx1 = indexes[idx]
            results = list()
            for img2 in files[idx:]:
                results.append(
                    self.calc_alg(img1, img2)
                )
            if last_ids:
                df.loc[idx1, indexes[idx:]] = results
                df.loc[idx1, last_ids] = df.loc[last_ids, idx1]
            else:
                df.loc[idx1, :] = results
            last_ids.append(idx1)

        return DistanceStruct(names=indexes, matrix=1.0-df.to_numpy(numpy.float64))
