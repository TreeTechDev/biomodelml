import pandas
import numpy
from functools import lru_cache
from src.variants.variant import Variant
from src.structs import DistanceStruct
from src.variants.deep_search.feature_extractor import FeatureExtractor
from src.variants.deep_search.indexer import Indexer

class DeepSearchVariant(Variant):

    name = "Deep Search with Annoy"

    def __init__(self, fasta_file: str, sequence_type: str, image_folder: str):
        super().__init__(fasta_file, sequence_type)
        self._image_folder = image_folder
        self._input_shape = (2000, 2000, 3)

    def calc_alg(self, img_name1: str, img_name2: str) -> float:
        distance = self.build_matrix()
        id1 = distance.names.index(img_name1)
        id2 = distance.names.index(img_name2)
        return distance.matrix[id1, id2]

    @lru_cache(maxsize=None)
    def build_matrix(self) -> DistanceStruct:
        features = FeatureExtractor(self._input_shape)
        indexer = Indexer(self._image_folder, self._names, features)
        names = [".".join(name.split("/")[-1].split(".")[:-1]) for name in indexer.image_list]
        diff = set(self._names).difference(set(names))
        if diff:
            raise IOError(f"Sequences without image created: {diff}")
        data = indexer.build()
        df = pandas.DataFrame()
        for i in range(len(self._names)):
            index_list = indexer.search_by_item(i)
            images = data.iloc[index_list[0]]['images_paths'].to_list()
            df = pandas.concat([df,
                pandas.DataFrame(
                    zip(images, index_list[1]),
                    columns=["gene", images[0]]
                ).set_index("gene")
            ], axis=1)
        df = df.sort_index(axis=1).sort_index(axis=0)

        return DistanceStruct(names=names, matrix=df.to_numpy(numpy.float64))
