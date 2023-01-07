import os
import tensorflow
import numpy
import cv2
import itertools
import pickle
from pathlib import Path
from multiprocessing import Pool, cpu_count
from typing import List

folder = Path("/tmp")
def read_image(img_folder: str, img_name: str):
    im = cv2.imread(os.path.join(
                    img_folder, img_name))
    path = folder / img_folder / img_name
    path.parent.mkdir(exist_ok=True, parents=True)
    with open(path, "wb") as f:
        numpy.save(f, cv2.cvtColor(im, cv2.COLOR_BGR2RGB))   # BGR -> RGB
    return path

def read_all_images(folders: List[str]):
    img_dict = dict()
    for folder in folders:
        for f in os.listdir(folder):
            img_dict[f"{folder.split('/')[-2]}_{f}"] = read_image(folder, f)
    return img_dict

items = read_all_images([
    "data/images/orthologs_androglobin/full/",
    "data/images/orthologs_cytoglobin/full/",
    "data/images/orthologs_hemoglobin_beta/full/",
    "data/images/orthologs_myoglobin/full/",
    "data/images/orthologs_neuroglobin/full/"
])

from src.variants.ssim import SSIMVariant, DEFAULT_PARAMS, RecursionContext


class SSIMSearch(SSIMVariant):
    def __init__(self, **alg_params):
        self._alg_params = DEFAULT_PARAMS
        self._alg_params.update(alg_params)
        self.filter_size = self._alg_params["filter_size"]

ssim = SSIMSearch(filter_size=7, filter_sigma=0.7)

def _metric(img_a, img_b):
    with open(img_a, "rb") as f_a:
        image_a = tensorflow.expand_dims(tensorflow.convert_to_tensor(numpy.array(numpy.load(f_a))), axis=0)
    with open(img_b, "rb") as f_b:
        image_b = tensorflow.expand_dims(tensorflow.convert_to_tensor(numpy.array(numpy.load(f_b))), axis=0)
    with RecursionContext():
        return ssim._match_images(image_a, image_b)[0]

def _fit(item_a, item_b):
    print(str(item_a),str(item_b))
    return _metric(item_a, item_b)


class AdaptedNearestNeighbors():
    def __init__(self, nproc=cpu_count()):
        self._items = {}
        self._i_and_d = {}
        self._nproc = nproc

    def fit(self, items: dict):
        self._items = items
        self._i_and_d = {str(k): dict() for k in items.values()}
        print(f"training for {len(items)}")
        item_by_item = list(itertools.combinations(items.values(), 2))
        with Pool(self._nproc) as pool:
            result = pool.starmap(
                _fit,
                item_by_item)
        for i, (item_a, item_b) in enumerate(item_by_item):
            item_a = str(item_a)
            item_b = str(item_b)
            self._i_and_d[item_a][item_b] = self._i_and_d[item_b][item_a] = result[i]
    
    def predict(self, item):
        if self._i_and_d.get(item):
            close_item = item
        else:  
            closest = 1
            close_item = None
            for i in self._items:
                prediction = _metric(self._items[i], item)
                if prediction <= closest:
                    closest = prediction
                    close_item = str(self._items[i])
                    if closest == 0:
                        break
        return dict(sorted(self._i_and_d[close_item].items(), key=lambda item: item[1]))

ann = AdaptedNearestNeighbors()

if Path("backup.pkl").exists():
    with open("backup.pkl", 'rb') as f:
        ann = pickle.load(f)
else:
    ann.fit(items)
    with open("backup.pkl", 'wb') as f:
        pickle.dump(ann, f)

    with open("p.txt", 'w') as f:
        f.write(str(ann._i_and_d))

r = ann.predict(read_image("data/images/orthologs_hemoglobin_beta/full/", "Carlito_syrichta_ENSTSYP00000007411_Csyr.png"))

print(r)
with open("r.txt", 'w') as f:
    f.write(str(r))
