import os
import pandas
import numpy
from PIL import Image
from sewar.full_ref import uqi
from multiprocessing import cpu_count
from multiprocessing.dummy import Pool
from src.variants.variant import Variant
from src.structs import DistanceStruct


class UQIVariant(Variant):

    name = "Universal Quality Index"

    def __init__(self, fasta_file, image_folder: str):
        super().__init__(fasta_file)
        self._image_folder = image_folder
    
    def _read_image(self, img_name: str) -> numpy.ndarray:
        with Image.open(
                os.path.join(
                    self._image_folder, img_name)
        ) as img:
            return numpy.asarray(img)        

    def _calc_uqi(self, img_name1: str, img_name2: str) -> float:
        return 1.0 - abs(uqi(
            self._read_image(img_name1), self._read_image(img_name2)))

    def build_matrix(self) -> DistanceStruct:
        files = os.listdir(self._image_folder)
        indexes = [img.split('.')[0] for img in files]
        threads = cpu_count() - 3
        diff = set(self._names).difference(set(indexes))
        if diff:
            raise IOError(f"Sequences without image created: {diff}")

        df = pandas.DataFrame(index=indexes, columns=indexes)
        last_ids = list()
        for idx, img1 in enumerate(files):
            idx1 = img1.split('.')[0]
            with Pool(threads) as pool:
                result = pool.starmap(
                    self._calc_uqi,
                    [(img1, img2) for img2 in files[idx:]]
                )
            if last_ids:
                df.loc[idx1, indexes[idx:]] = result
                df.loc[idx1, last_ids] = df.loc[last_ids, idx1]
            else:
                df.loc[idx1, :] = result
            last_ids.append(idx1)

        return DistanceStruct(names=indexes, matrix=df.to_numpy())
