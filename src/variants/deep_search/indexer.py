import os
import pandas
import numpy
from typing import List
from src.variants.deep_search.feature_extractor import FeatureExtractor
from annoy import AnnoyIndex


class Indexer:
    distance_type = "euclidean"

    def __init__(self, image_folder: str, seq_names: str, feature_extractor: FeatureExtractor):
        self.image_list = []
        self._index = None
        for path in os.listdir(image_folder):
            if path.split(".")[-1] in ("png", "jpg", "jpeg") and ".".join(path.split(".")[:-1]) in seq_names:
                self.image_list.append(os.path.join(image_folder, path))
        self._feature_extractor = feature_extractor

    def _feature_extraction(self) -> numpy.ndarray:
        image_data = pandas.DataFrame()
        image_data['images_paths'] = self.image_list
        f_data = self._feature_extractor.get_feature(self.image_list)
        image_data['features']  = f_data
        image_data = image_data.dropna().reset_index(drop=True)
        return image_data
    
    def _construct(self, data: numpy.ndarray):
        len_index = len(data['features'][0])
        self._index = AnnoyIndex(len_index, self.distance_type)
        for i, v in zip(data.index, data['features']):
            self._index.add_item(i, v)
        trees = len(data["features"])
        self._index.build(trees)

    def build(self) -> numpy.ndarray:
        data = self._feature_extraction()
        self._construct(data)
        return data
    
    def search_by_item(self, item: int) -> List[float]:
        if not self._index:
            raise IOError("Run build first")
        return self._index.get_nns_by_item(
            item, -1, search_k=-1, include_distances=True)
