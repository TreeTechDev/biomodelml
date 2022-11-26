#!/usr/bin/env python
# coding: utf-8

import sys
import optuna
from ete3 import Tree
from pathlib import Path
from src.variants.control import ControlVariant
from src.variants.ssim_multiscale import SSIMVariant
from src.experiment import Experiment

SEED = 42

def control():
    data_path = Path(sys.argv[1])
    seq = sys.argv[2]
    fasta_file = data_path / f"{seq}.fasta.sanitized"
    experiment = Experiment(
        data_path,
        ControlVariant(fasta_file, "N")
    ).run()
    tree_struct = experiment._trees[0]
    newick_ctl = tree_struct.tree.to_newick(
        labels=tree_struct.distances.names,
        include_distance=False)    
    return Tree(newick_ctl, format=1)

CONTROL_TREE = control()


def objective(trial):
    data_path = Path(sys.argv[1])
    seq = sys.argv[2]

    params = dict(
        filter_sigma=trial.suggest_float("filter_sigma", 0.1, 1.5, step=0.1),
        filter_size=trial.suggest_int("filter_size", 3, 15)
    )
    fasta_file = data_path / f"{seq}.fasta.sanitized"
    experiment = Experiment(
        data_path,
        SSIMVariant(
            fasta_file,
            "N",
            data_path / "images" / seq / "full",
            **params
        )
    ).run()
    
    tree_struct = experiment._trees[0]
    newick_alg = tree_struct.tree.to_newick(
        labels=tree_struct.distances.names,
        include_distance=False)
    tree = Tree(newick_alg, format=1)
    return CONTROL_TREE.compare(tree, unrooted=True)["norm_rf"]

if __name__ == "__main__":

    seq = sys.argv[2]
    study = optuna.create_study(
        storage=f"sqlite:////app/data/{seq}.db",
        study_name="bioinfo",
        load_if_exists= True,
        direction="minimize",
        sampler=optuna.samplers.TPESampler(seed=SEED),
        pruner=optuna.pruners.MedianPruner(n_warmup_steps=10)
    )
    try:
        study.optimize(objective, n_trials=200)
    except:
        print("Best Params")
        print(study.best_params)
        print("Best Values")
        print(study.best_value)