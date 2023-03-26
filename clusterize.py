import os
import itertools
import pickle
from pathlib import Path
from multiprocessing import Pool, cpu_count
from typing import List
from src.variants.ssim import SSIMVariant, DEFAULT_PARAMS
from src.variants.ssim_multiscale import SSIMMultiScaleVariant
from src.variants.uqi import UQIVariant


def read_all_images(folders: List[str]):
    img_dict = dict()
    for folder in folders:
        for f in os.listdir(folder):
            img_dict[f"{folder.split('/')[-2]}_{f}"] = os.path.join(folder, f)
    return img_dict

items = read_all_images([
    "data/images/indelible/full/",
    "data/images/orthologs_androglobin/full/",
    "data/images/orthologs_cytoglobin/full/",
    "data/images/orthologs_hemoglobin_beta/full/",
    "data/images/orthologs_myoglobin/full/",
    "data/images/orthologs_neuroglobin/full/"
])

class SSIMSearch(SSIMVariant):
    def __init__(self):
        self._image_folder = ""
        self._alg_params = DEFAULT_PARAMS
        self._names = [".".join(name.split("/")[-1].split(".")[:-1]) for name in items.values()]


class SSIMMultiScaleSearch(SSIMMultiScaleVariant):
    def __init__(self):
        self._image_folder = ""
        self._alg_params = DEFAULT_PARAMS
        self._names = [".".join(name.split("/")[-1].split(".")[:-1]) for name in items.values()]

class UQISearch(UQIVariant):
    def __init__(self):
        self._image_folder = ""
        self._names = [".".join(name.split("/")[-1].split(".")[:-1]) for name in items.values()]

ssim = SSIMSearch()
msssim = SSIMMultiScaleSearch()
uqi = UQISearch()

ALGORITMS = (ssim, msssim, uqi)


def _metric(img_a, img_b):
    results = dict()
    for alg in ALGORITMS:
        results[alg.name] = alg.calc_alg(img_a, img_b)
    return results

def _fit(item_a, item_b):
    print(str(item_a),str(item_b))
    return _metric(item_a, item_b)


class AdaptedNearestNeighbors():
    def __init__(self, algs, nproc=cpu_count()):
        self._items = {}
        self._algs = algs
        self._i_and_d = {}
        self._nproc = nproc

    def fit(self, items: dict):
        self._items = items
        self._i_and_d = {a.name:{str(k): dict() for k in items.values()} for a in self._algs}
        print(f"training for {len(items)} items and {len(self._algs)} algoritms")
        item_by_item = list(itertools.combinations(items.values(), 2))
        with Pool(self._nproc) as pool:
            result = pool.starmap(
                _fit,
                item_by_item)
        for i, (item_a, item_b) in enumerate(item_by_item):
            item_a = str(item_a)
            item_b = str(item_b)
            for alg in self._algs:
                self._i_and_d[alg.name][item_a][item_b] = self._i_and_d[alg.name][item_b][item_a] = result[i][alg.name]
    

ann = AdaptedNearestNeighbors(ALGORITMS)

if not Path("data/backup_cluster.pkl").exists():
    ann.fit(items)
    with open("data/backup_cluster.pkl", 'wb') as f:
        pickle.dump(ann, f)

    with open("data/cluster_sim.pkl", 'wb') as f:
        pickle.dump(ann._i_and_d, f)

with open("data/final_cluster.csv", "w") as f:
    f.write("Algoritm, Name, Family, Right, Total\n")
with open("data/cluster_sim.pkl", "rb") as f:
    all_hash = pickle.load(f)
    for alg, results in all_hash.items():
        for k, v in results.items():
            family = k.split("/")[-3]
            name = k.split("/")[-1]
            is_right = 0
            total = 0
            sorted_items = sorted(v.items(), key=lambda item: item[1], reverse=True) 
            for i, v in sorted_items:
                if family == i.split("/")[-3]:
                    total += 1
            for i, v in sorted_items[:total]:
                if family == i.split("/")[-3]:
                    is_right += 1
            result = f"{alg},{name},{family},{is_right},{total}\n"
            print(result)
            with open("data/final_cluster.csv", "a") as f:
                f.write(result) 
