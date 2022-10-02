#!/usr/bin/env python
# coding: utf-8

import sys
import optuna
import pandas
from ete3 import Tree
from pathlib import Path
from src.variants.control import ControlVariant
from src.variants.ssim import SSIMVariant
from src.experiment import Experiment

DEFAULT_ALG = "Structural Similarity Index Measure"
CONTROL_ALG = "Control with Clustal Omega"
SEED = 42


def objective(trial):
    data_path = Path(sys.argv[1])
    seq = sys.argv[2]

    params = dict(
        filter_sigma=trial.suggest_float("filter_sigma", 0.0, 2.0),
        k1=trial.suggest_float("k1", 0.0, 0.1),
        k2=trial.suggest_float("k2", 0.0, 0.1),
        filter_size=trial.suggest_int("filter_size", 3, 12)
    )
    fasta_file = data_path / f"{seq}.fasta.sanitized"
    experiment = Experiment(
        data_path,
        ControlVariant(fasta_file, "N"),
        SSIMVariant(
            fasta_file,
            "N",
            data_path / "images" / seq / "full",
            **params
        )
    ).run()
    for tree_struct in experiment._trees:
        if tree_struct.name == DEFAULT_ALG:
            newick_alg = tree_struct.tree.to_newick(
                labels=tree_struct.distances.names,
                include_distance=False)
        elif tree_struct.name == CONTROL_ALG:
            newick_ctl = tree_struct.tree.to_newick(
                labels=tree_struct.distances.names,
                include_distance=False)    
    control_tree = Tree(newick_ctl, format=1)
    tree = Tree(newick_alg, format=1)
    return control_tree.compare(tree, unrooted=True)["norm_rf"]

if __name__ == "__main__":

    seq = sys.argv[2]
    study = optuna.create_study(
        storage=f"sqlite:////app/{seq}.db",
        study_name="bioinfo",
        load_if_exists= True,
        direction="minimize",
        sampler=optuna.samplers.TPESampler(seed=SEED),
        pruner=optuna.pruners.MedianPruner(n_warmup_steps=10)
    )
    study.optimize(objective, n_trials=500)
    print("Best Params")
    print(study.best_params)
    print("Best Values")
    print(study.best_value)