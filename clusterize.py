import os
import itertools
import pickle
from pathlib import Path
from multiprocessing import Pool, cpu_count
from typing import List
from src.variants.resized_ssim import ResizedSSIMVariant
from src.variants.resized_ssim_multiscale import ResizedSSIMMultiScaleVariant
from src.variants.windowed_ssim_multiscale import WindowedSSIMMultiScaleVariant
from src.variants.greedy_ssim import GreedySSIMVariant
from src.variants.unrestricted_ssim import UnrestrictedSSIMVariant
from src.variants.uqi import UQIVariant

TYPES_DICT = dict(P=dict(), N=dict())

def read_all_images(folders: List[str]):
    img_dict = TYPES_DICT.copy()
    for folder in folders:
        for t in img_dict.keys():
            formated_folder = folder.format(t)
            for f in os.listdir(formated_folder):
                img_dict[t][f"{formated_folder.split('/')[-2]}_{f}"] = os.path.join(formated_folder, f)
    return img_dict

items = read_all_images([
    "data/images/{}/indelible/full/",
    "data/images/{}/orthologs_androglobin/full/",
    "data/images/{}/orthologs_cytoglobin/full/",
    "data/images/{}/orthologs_hemoglobin_beta/full/",
    "data/images/{}/orthologs_myoglobin/full/",
    "data/images/{}/orthologs_neuroglobin/full/"
])


item_list = []
for types in TYPES_DICT.keys():
    item_list += [".".join(name.split("/")[-1].split(".")[:-1]) for name in items[types].values()]

rssim = ResizedSSIMVariant
rmsssim = ResizedSSIMMultiScaleVariant
wmsssim = WindowedSSIMMultiScaleVariant
gssim = GreedySSIMVariant
ussim = UnrestrictedSSIMVariant
uqi = UQIVariant

ALGORITMS = (
    ResizedSSIMVariant,
    ResizedSSIMMultiScaleVariant,
    WindowedSSIMMultiScaleVariant,
    GreedySSIMVariant,
    UnrestrictedSSIMVariant,
    UQIVariant
)


def _metric(img_a, img_b, types):
    results = dict()
    for alg in ALGORITMS:
        results[alg.name] = alg.from_name_list(item_list, sequence_type=types).calc_alg(img_a, img_b)
    return results

def _fit(item_a, item_b, types):
    print(str(item_a),str(item_b))
    return _metric(item_a, item_b, types)


class AdaptedNearestNeighbors():
    def __init__(self, algs, nproc=cpu_count()):
        self._items = {}
        self._algs = algs
        self._i_and_d = {}
        self._nproc = nproc

    def fit(self, items: dict):
        self._items = items
        self._i_and_d = dict(
            P={a.name:{str(k): dict() for k in items["P"].values()} for a in self._algs},
            N={a.name:{str(k): dict() for k in items["N"].values()} for a in self._algs})
        print(f"training for {len(items['P'])} protein and {len(items['N'])} nucleotide sequences with {len(self._algs)} algoritms")
        for types in ["P", "N"]:
            print(f"doing for {types}")
            item_by_item = [(ia, ib, types) for ia, ib in itertools.combinations(items[types].values(), 2)]
            with Pool(self._nproc) as pool:
                result = pool.starmap(
                    _fit,
                    item_by_item)
            for i, (item_a, item_b, types) in enumerate(item_by_item):
                item_a = str(item_a)
                item_b = str(item_b)
                for alg in self._algs:
                    self._i_and_d[types][alg.name][item_a][item_b] = self._i_and_d[types][alg.name][item_b][item_a] = result[i][alg.name]
    

ann = AdaptedNearestNeighbors(ALGORITMS)

if not Path("data/backup_cluster.pkl").exists():
    ann.fit(items)
    with open("data/backup_cluster.pkl", 'wb') as f:
        pickle.dump(ann, f)

    with open("data/cluster_sim.pkl", 'wb') as f:
        pickle.dump(ann._i_and_d, f)

with open("data/final_cluster.csv", "w") as f:
    f.write("Type,Algoritm,Name,Family,Right,Total\n")
with open("data/cluster_sim.pkl", "rb") as f:
    all_hash = pickle.load(f)
    for t in all_hash:
        for alg in all_hash[t]:
            for k, v in all_hash[t][alg].items():
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
                result = f"{t},{alg},{name},{family},{is_right},{total}\n"
                print(result)
                with open("data/final_cluster.csv", "a") as f:
                    f.write(result) 
