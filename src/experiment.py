from __future__ import annotations
from os import path
from typing import Iterable
from src.variants.variant import Variant
from src.structs import TreeStruct
from biotite.sequence import phylo
from matplotlib import pyplot
from Bio import Phylo
from io import StringIO


class Experiment:
    def __init__(self, output_path: str, *variants: Iterable[Variant]):
        self._output_path = output_path
        self._variants = variants
        self._trees = []
    
    def run(self) -> Experiment:
        self._trees = []
        for variant in self._variants:
            distances = variant.build_matrix()
            tree = phylo.neighbor_joining(distances.matrix)
            self._trees.append(
                TreeStruct(name=variant.name, names=distances.names, tree=tree))
        return self
    
    def save(self):
        fig = pyplot.figure(figsize=(12.0, 12.0))
        for tree_struct in self._trees:
            fig.suptitle(tree_struct.name, fontsize=16)
            ax = fig.gca()
            ax.axis("off")
            t = Phylo.read(
                StringIO(tree_struct.tree.to_newick()), "newick", values_are_confidence=True)
            t.ladderize()
            Phylo.draw(
                t,
                show_confidence=True,
                axes=ax,
                do_show=False,
                label_func=lambda clade: "" if not clade.name else tree_struct.names[int(clade.name)],
                branch_labels=lambda clade: "" if not clade.name else "{:.2f}".format(clade.confidence)
            )
            pyplot.savefig(
                path.join(self._output_path, f"{tree_struct.name}.png"))
