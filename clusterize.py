import os
import sys
import itertools
import pickle
from pathlib import Path
from multiprocessing import Pool, cpu_count
from typing import List
from signal import signal, SIGABRT, SIGSTOP, SIGILL, SIGINT, SIGSEGV, SIGTERM
from src.variants.deep_search.variant import DeepSearchVariant
from src.variants.resized_ssim import ResizedSSIMVariant
from src.variants.resized_ssim_multiscale import ResizedSSIMMultiScaleVariant
from src.variants.windowed_ssim_multiscale import WindowedSSIMMultiScaleVariant
from src.variants.greedy_ssim import GreedySSIMVariant
from src.variants.unrestricted_ssim import UnrestrictedSSIMVariant
from src.variants.uqi import UQIVariant

TYPES_DICT = dict(P=dict(), N=dict())
PARTIAL = True
IMG_FOLDER = "data/images"

def read_all_images(folders: List[str]):
    img_dict = TYPES_DICT.copy()
    for folder in folders:
        for t in img_dict.keys():
            formated_folder = folder.format(t)
            for f in os.listdir(formated_folder):
                img_dict[t][f"{formated_folder.split('/')[-2]}_{f}"] = os.path.join(formated_folder, f)
    return img_dict

items = read_all_images([
    IMG_FOLDER + "/{}/indelible/full/",
    IMG_FOLDER + "/{}/orthologs_androglobin/full/",
    IMG_FOLDER + "/{}/orthologs_cytoglobin/full/",
    IMG_FOLDER + "/{}/orthologs_hemoglobin_beta/full/",
    IMG_FOLDER + "/{}/orthologs_myoglobin/full/",
    IMG_FOLDER + "/{}/orthologs_neuroglobin/full/"
])


item_list = []
for types in TYPES_DICT.keys():
    item_list += [".".join(name.split("/")[-1].split(".")[:-1]) for name in items[types].values()]

ALGORITMS = (
    DeepSearchVariant(image_folder=IMG_FOLDER).load_or_build(item_list),
    ResizedSSIMVariant,
    ResizedSSIMMultiScaleVariant,
    WindowedSSIMMultiScaleVariant,
    GreedySSIMVariant,
    UnrestrictedSSIMVariant,
    UQIVariant
)


def _metric(img_a, img_b, types, alg_name):
    results = dict()
    for alg in ALGORITMS:
        if alg.name == alg_name:
            results[alg.name] = alg.from_name_list(item_list, sequence_type=types).calc_alg(img_a, img_b)
    return results

def _fit(item_a, item_b, types, alg_name):
    print(str(item_a),str(item_b), alg_name)
    return _metric(item_a, item_b, types, alg_name)


def save_checkpoint(ann):
    all_hash = ann._i_and_d
    with open("data/cluster_sim.pkl", 'wb') as f:
        pickle.dump(all_hash, f)

    with open("data/backup_cluster.pkl", 'wb') as f:
        pickle.dump(ann, f)

def checkpoint(ann):
    def internal_check(*args):
        save_checkpoint(ann)
        sys.exit(0)
    return internal_check


class AdaptedNearestNeighbors():
    def __init__(self, algs, i_and_d={}, nproc=cpu_count()):
        self._items = {}
        self._algs = algs
        self._i_and_d = i_and_d
        self._nproc = nproc

    def fit(self, items: dict):
        self._items = items
        self._i_and_d = self._i_and_d or dict(
            P={a.name:{str(k): dict() for k in items["P"].values()} for a in self._algs},
            N={a.name:{str(k): dict() for k in items["N"].values()} for a in self._algs})
        print(f"training for {len(items['P'])} protein and {len(items['N'])} nucleotide sequences with {len(self._algs)} algoritms")
        for types in ["P", "N"]:
            print(f"doing for {types}")
            item_by_item = []
            for ia, ib in itertools.combinations(items[types].values(), 2):
                for alg in self._algs:
                    if not self._i_and_d[types].get(alg.name, {}).get(ia, {}).get(ib, {}):
                        item_by_item.append((ia, ib, types, alg.name))
            print(f"filtered and combined pairwise now just with {len(item_by_item)} combinations for {types}")
            with Pool(self._nproc) as pool:
                result = pool.starmap(
                    _fit,
                    item_by_item)
            for i, (item_a, item_b, types, alg_name) in enumerate(item_by_item):
                item_a = str(item_a)
                item_b = str(item_b)
                self._i_and_d[types][alg_name][item_a][item_b] = self._i_and_d[types][alg_name][item_b][item_a] = result[i][alg_name]
    


if PARTIAL and Path("data/cluster_sim.pkl").exists():
    with open("data/cluster_sim.pkl", "rb") as f:
        all_hash = pickle.load(f)
        ann = AdaptedNearestNeighbors(ALGORITMS, all_hash)
else:
    ann = AdaptedNearestNeighbors(ALGORITMS)

for sig in (SIGABRT, SIGILL, SIGINT, SIGSEGV, SIGTERM):
    signal(sig, checkpoint(ann))

ann.fit(items)
all_hash = ann._i_and_d

checkpoint(ann)()

with open("data/final_cluster.csv", "w") as f:
    f.write("Type,Algorithm,Name,Family,Right,Total\n")
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
