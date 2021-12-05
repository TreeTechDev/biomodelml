import numpy
from typing import Tuple, List
from tensorflow.keras.applications.vgg16 import VGG16, preprocess_input
from tensorflow.keras.preprocessing import image
from PIL import Image

class FeatureExtractor:
    def __init__(self, input_shape: Tuple[int]):
        self._input_shape = input_shape
        self.model = VGG16(
            include_top=False, input_shape=self._input_shape, weights='imagenet')

    def extract(self, img: Image) -> numpy.ndarray:
        img = img.resize((self._input_shape[0], self._input_shape[1]))
        img = img.convert('RGB')
        x = image.img_to_array(img)
        x = numpy.expand_dims(x, axis=0)
        x = preprocess_input(x)
        feature = self.model.predict(x)[0]
        return (feature / numpy.linalg.norm(feature)).reshape((feature.size, 1))

    def get_feature(self, image_data: List[str]):
        self.image_data = image_data 
        features = []
        for img_path in self.image_data:
            try:
                feature = self.extract(img=Image.open(img_path))
                features.append(feature)
            except Exception as e:
                raise e
        return features