#!/usr/bin/env python
# coding: utf-8

import os
import sys
import ppscore
import pandas
import numpy
from scipy.stats import spearmanr
from ete3 import Tree

#  all_channels = ("full", "red", "green", "blue", "red_green", "red_blue", "green_blue", "gray_r", "gray_g", "gray_b", "gray_max", "gray_mean")


def read_and_compare(tree_path: str, dataset: str, channel: str = "full"):
    result_dict = {}
    path = f"{tree_path}/{channel}/{dataset}"
    control_df = pandas.read_csv(f"{path}/Control with Clustal Omega.csv", index_col=0)
    for file in os.listdir(path):
        if file.endswith(".csv"):
            basename = ".".join(file.join(".")[0:-1])
            df = pandas.read_csv(f"{path}/{file}", index_col=0)
            result_dict[basename] = round(numpy.sqrt(numpy.sum((control_df.values - df.values)**2)), 4)
    result_df = pandas.DataFrame(result_dict, index=[dataset])
    result_df.rename(columns={"MultiScale Structural Similarity Index Measure": channel}, inplace=True)
    result_df.index.name = "dataset"
    return result_df


def euclidean(tree_path: str, dataset: str):
    channels = ("full", "gray_mean")
    datasets = [dataset] if dataset else ["orthologs_hemoglobin_beta", "orthologs_myoglobin", "orthologs_neuroglobin", "orthologs_cytoglobin", "orthologs_androglobin"]
    print(f"Euclidean Distance from distance matrices for datasets: {','.join(datasets)}")
    sum_dfs = []
    for channel in channels:
        dfs = []
        for dataset in datasets:
            dfs.append(read_and_compare(tree_path, dataset, channel))
        sum_dfs.append(pandas.concat(dfs))
    return pandas.concat(sum_dfs, axis=1).T.reset_index().drop_duplicates().set_index("index").T


def read_and_linear_correlate(tree_path: str, dataset: str, channel: str = "full"):
    result_dict = {}
    path = f"{tree_path}/{channel}/{dataset}"
    control_df = pandas.read_csv(f"{path}/Control with Clustal Omega.csv", index_col=0)
    for file in os.listdir(path):
        if file.endswith(".csv"):
            basename = ".".join(file.split(".")[0:-1])
            df = pandas.read_csv(f"{path}/{file}", index_col=0)
            result_dict[basename] = round(spearmanr(control_df.values.flatten(), df.values.flatten())[0], 4)
    result_df = pandas.DataFrame(result_dict, index=[dataset])
    result_df.rename(columns={"MultiScale Structural Similarity Index Measure": channel}, inplace=True)
    result_df.index.name = "dataset"
    return result_df


def linear_correlation(tree_path: str, dataset: str):
    channels = ("full", "gray_mean")
    datasets = [dataset] if dataset else ["orthologs_hemoglobin_beta", "orthologs_myoglobin", "orthologs_neuroglobin", "orthologs_cytoglobin", "orthologs_androglobin"]
    print(f"Linear Correlation from distance matrices for datasets: {','.join(datasets)}")
    sum_dfs = []
    for channel in channels:
        dfs = []
        for dataset in datasets:
            dfs.append(read_and_linear_correlate(tree_path, dataset, channel))
        sum_dfs.append(pandas.concat(dfs))
    return pandas.concat(sum_dfs, axis=1).T.reset_index().drop_duplicates().set_index("index").T


def read_and_correlate(tree_path: str, dataset: str, channel: str = "full"):
    result_dict = {}
    path = f"{tree_path}/{channel}/{dataset}"
    control_df = pandas.read_csv(f"{path}/Control with Clustal Omega.csv", index_col=0)
    for file in os.listdir(path):
        if file.endswith(".csv"):
            basename = ".".join(file.split(".")[0:-1])
            df = pandas.read_csv(f"{path}/{file}", index_col=0)
            result_dict[basename] = round(ppscore.score(
                pandas.DataFrame({"x": control_df.values.flatten(), "y": df.values.flatten()}), "x", "y")["ppscore"], 4)
    result_df = pandas.DataFrame(result_dict, index=[dataset])
    result_df.rename(columns={"MultiScale Structural Similarity Index Measure": channel}, inplace=True)
    result_df.index.name = "dataset"
    return result_df


def pps(tree_path: str, dataset: str):
    channels = ("full", "gray_mean")
    datasets = [dataset] if dataset else ["orthologs_hemoglobin_beta", "orthologs_myoglobin", "orthologs_neuroglobin", "orthologs_cytoglobin", "orthologs_androglobin"]
    print(f"Predictive Power Score from distance matrices for datasets: {','.join(datasets)}")
    sum_dfs = []
    for channel in channels:
        dfs = []
        for dataset in datasets:
            dfs.append(read_and_correlate(tree_path, dataset, channel))
        sum_dfs.append(pandas.concat(dfs))
    return pandas.concat(sum_dfs, axis=1).T.reset_index().drop_duplicates().set_index("index").T


def read_and_tree_compare(tree_path: str, dataset: str, channel: str = "full"):
    result_dict = {}
    path = f"{tree_path}/{channel}/{dataset}"
    control_tree = Tree(f"{path}/Control with Clustal Omega.nw", format=1)
    for file in os.listdir(path):
        if file.endswith(".nw"):
            basename = ".".join(file.split(".")[0:-1])
            tree = Tree(f"{path}/{file}", format=1)
            result = control_tree.compare(tree, unrooted=True)
            result_dict[basename] = round(result["norm_rf"], 4)
    result_df = pandas.DataFrame(result_dict, index=[dataset])
    result_df.rename(columns={"MultiScale Structural Similarity Index Measure": channel}, inplace=True)
    result_df.index.name = "dataset"
    return result_df


# ## Robinson-foulds distance from generated Trees
def robinson_foulds(tree_path: str, dataset: str):
    channels = ("full", "gray_mean")
    datasets = [dataset] if dataset else ["orthologs_hemoglobin_beta", "orthologs_myoglobin", "orthologs_neuroglobin", "orthologs_cytoglobin", "orthologs_androglobin"]
    print(f"Robinson-foulds distance from generated Trees for datasets: {','.join(datasets)}")
    sum_dfs = []
    for channel in channels:
        dfs = []
        for dataset in datasets:
            dfs.append(read_and_tree_compare(tree_path, dataset, channel))
        sum_dfs.append(pandas.concat(dfs))
    return pandas.concat(sum_dfs, axis=1).T.reset_index().drop_duplicates().set_index("index").T


def read_and_tree_compare_branches(tree_path: str, dataset: str, channel: str = "full"):
    result_dict = {}
    path = f"{tree_path}/{channel}/{dataset}"
    control_tree = Tree(f"{path}/Control with Clustal Omega.nw", format=1)
    for file in os.listdir(path):
        if file.endswith(".nw"):
            basename = ".".join(file.split(".")[0:-1])
            tree = Tree(f"{path}/{file}", format=1)
            result = control_tree.compare(tree, unrooted=True)
            result_dict[basename] = 1.0 - round(result["source_edges_in_ref"], 4)
    result_df = pandas.DataFrame(result_dict, index=[dataset])
    result_df.rename(columns={"MultiScale Structural Similarity Index Measure": channel}, inplace=True)
    result_df.index.name = "dataset"
    return result_df


def branch_score(tree_path: str, dataset: str = None):
    channels = ("full", "gray_mean")
    datasets = [dataset] if dataset else ["orthologs_hemoglobin_beta", "orthologs_myoglobin", "orthologs_neuroglobin", "orthologs_cytoglobin", "orthologs_androglobin"]
    print(f"Compatibility branch score from generated Trees for datasets: {','.join(datasets)}")
    sum_dfs = []
    for channel in channels:
        dfs = []
        for dataset in datasets:
            dfs.append(read_and_tree_compare_branches(tree_path, dataset, channel))
        sum_dfs.append(pandas.concat(dfs))
    return pandas.concat(sum_dfs, axis=1).T.reset_index().drop_duplicates().set_index("index").T


def main(tree_path, dataset):
    print(euclidean(tree_path, dataset))
    print(linear_correlation(tree_path, dataset))
    print(pps(tree_path, dataset))
    print(robinson_foulds(tree_path, dataset))
    print(branch_score(tree_path, dataset))


if __name__ == "__main__":
    tree_path = sys.argv[1]
    dataset = sys.argv[2]
    main(tree_path, dataset)

