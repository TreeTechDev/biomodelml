import os
import pandas
import numpy
import tensorflow
from tensorflow import Tensor
from typing import Tuple
from src.variants.variant import Variant
from src.structs import DistanceStruct



class SSIMVariant(Variant):

    name = "Structural Similarity Index Measure"

    def __init__(self, fasta_file: str, sequence_type: str, image_folder: str):
        super().__init__(fasta_file, sequence_type)
        self._image_folder = image_folder
    
    def _call_alg(self, image: Tensor, other: Tensor) -> numpy.ndarray:
        return tensorflow.image.ssim(
                                image,
                                other,
                                max_val=255, filter_size=11,
                                filter_sigma=1.5, k1=0.01, k2=0.03)[0].numpy()

    def _read_image(self, img_name: str) -> Tensor:
        return tensorflow.expand_dims(
            tensorflow.image.decode_image(
                tensorflow.io.read_file(
                    os.path.join(
                        self._image_folder, img_name))), axis=0)

    def _upscale_images(self, image: Tensor, other: Tensor) -> Tuple[Tensor]:
        max_x = image.shape[1] if image.shape[1] > other.shape[1] else other.shape[1]
        max_y = image.shape[2] if image.shape[2] > other.shape[2] else other.shape[2]
        return (
            tensorflow.image.resize(image, (max_x, max_y), tensorflow.image.ResizeMethod.BICUBIC),
            tensorflow.image.resize(other, (max_x, max_y), tensorflow.image.ResizeMethod.BICUBIC)
        )

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
                img, other = self._upscale_images(
                    self._read_image(img1),
                    self._read_image(img2)
                )
                results.append(
                    self._call_alg(img, other)
                )
            if last_ids:
                df.loc[idx1, indexes[idx:]] = results
                df.loc[idx1, last_ids] = df.loc[last_ids, idx1]
            else:
                df.loc[idx1, :] = results
            last_ids.append(idx1)

        return DistanceStruct(names=indexes, matrix=1.0-df.to_numpy(numpy.float64))
