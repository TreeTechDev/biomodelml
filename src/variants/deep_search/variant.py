import numpy
from functools import cache
from typing import List
from src.variants.variant import Variant
from src.structs import DistanceStruct
from src.variants.deep_search.feature_extractor import FeatureExtractor
from src.variants.deep_search.indexer import Indexer

class DeepSearchVariant(Variant):

    name = "Deep Search with Annoy"

    def __init__(self, fasta_file: str = None, sequence_type: str = None, image_folder: str = ""):
        super().__init__(fasta_file, sequence_type)
        self._image_folder = image_folder
        self._input_shape = (2000, 2000, 3)
        self._sequence_type = sequence_type
    
    @cache
    def _build_once(self):
        self._names = sorted(self._names)
        features = FeatureExtractor(self._input_shape)
        self.indexer = Indexer(self._image_folder, self._names, features)
        self.indexer.load_or_build()

    def calc_alg(self, img_name1: str, img_name2: str) -> float:
        self._build_once()
        img1_idx = self._names.index(img_name1)
        img2_idx = self._names.index(img_name2)
        result = self.indexer.get_distance(img1_idx, img2_idx)
        print(f"{img_name1} and {img_name2} done with score {result}")
        return result
    
    def load_or_build(self, names):
        self._names = names
        self._build_once()
        return self

    def from_name_list(self, name_list = List[str], **kwargs):
        return self

    def build_matrix(self) -> DistanceStruct:
        self.load_or_build(self._names)
        names = [".".join(name.split("/")[-1].split(".")[:-1]) for name in self.indexer.image_list]
        diff = set(self._names).difference(set(names))
        if diff:
            raise IOError(f"Sequences without image created: {diff}")
        matrix = numpy.zeros((len(self._names), len(self._names)))
        for i in range(len(self._names)):
            for j in range(i, len(self._names)):
                matrix[i,j] = matrix[j,i] = self.indexer.get_distance(i, j)

        return DistanceStruct(names=names, matrix=matrix)
