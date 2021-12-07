import os
import pandas
import tensorflow
from src.variants.variant import Variant
from src.structs import DistanceStruct



class SSIMMultiScaleVariant(Variant):

    name = "MultiScale Structural Similarity Index Measure"

    def __init__(self, fasta_file, image_folder: str):
        super().__init__(fasta_file)
        self._image_folder = image_folder
    
    def _read_image(self, img_name: str) -> tensorflow.Tensor:
        return tensorflow.expand_dims(
            tensorflow.image.decode_image(
                tensorflow.io.read_file(
                    os.path.join(
                        self._image_folder, img_name))), axis=0)

    def build_matrix(self) -> DistanceStruct:
        files = os.listdir(self._image_folder)
        indexes = [img.split('.')[0] for img in files]
        diff = set(self._names).difference(set(indexes))
        if diff:
            raise IOError(f"Sequences without image created: {diff}")

        df = pandas.DataFrame(index=indexes, columns=indexes)
        last_ids = list()
        for idx, img1 in enumerate(files):
            idx1 = img1.split('.')[0]
            results = list()
            for img2 in files[idx, :]:
                results.append(
                    tensorflow.image.ssim_multiscale(
                        self._read_image(img1),
                        self._read_image(img2),
                        max_val=255, filter_size=11,
                        filter_sigma=1.5, k1=0.01, k2=0.03)[0].numpy()
                )
            if last_ids:
                df.loc[idx1, indexes[idx:]] = results
                df.loc[idx1, last_ids] = df.loc[last_ids, idx1]
            else:
                df.loc[idx1, :] = results
            last_ids.append(idx1)

        return DistanceStruct(names=indexes, matrix=1-df.to_numpy())
