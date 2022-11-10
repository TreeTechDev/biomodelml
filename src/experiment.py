from __future__ import annotations
import traceback
import pandas
from pathlib import Path
from typing import Iterable
from src.variants.variant import Variant
from src.structs import TreeStruct
from biotite.sequence import phylo
from matplotlib import pyplot
from Bio import Phylo
from io import StringIO
from matplotlib.figure import Figure

class Experiment:
    def __init__(self, output_path: Path, *variants: Iterable[Variant]):
        self._output_path = output_path
        self._variants = variants
        self._trees = []

    def run(self) -> Experiment:
        self._trees = []
        for variant in self._variants:
            try:
                distances = variant.build_matrix()
                tree = phylo.neighbor_joining(distances.matrix)
                self._trees.append(
                    TreeStruct(name=variant.name, distances=distances, tree=tree))
            except Exception:
                print(traceback.format_exc())
        return self

    def _save_distance_matrix(self, tree_struct: TreeStruct):
        pandas.DataFrame(
            data=tree_struct.distances.matrix,
            index=tree_struct.distances.names,
            columns=tree_struct.distances.names).to_csv(
                self._output_path / f"{tree_struct.name}.csv"
            )

    def _save_img_positions(self, tree_struct: TreeStruct):
        if tree_struct.distances.img_positions:
            pos = "img1,img2,score,start_col,start_line,stop_col,stop_line,max_size\n"
            for items in tree_struct.distances.img_positions:
                img1, img2, positions = items
                for position in positions:
                    pos += f"{img1},{img2},{','.join(position)}\n"
            with open(self._output_path / f"{tree_struct.name}.map", "w") as f:
                f.write(pos)

    def _save_align(self, tree_struct: TreeStruct):
        if tree_struct.distances.align:
            align = ""
            for i, seq in enumerate(tree_struct.distances.align.get_gapped_sequences()):
                align += f">{tree_struct.distances.names[i]}\n"
                align += f"{seq}\n"
            with open(self._output_path / f"{tree_struct.name}.fasta", "w") as f:
                f.write(align)

    def _save_newick_tree(self, tree_struct: TreeStruct):
        newick = tree_struct.tree.to_newick(
            labels=tree_struct.distances.names, include_distance=False)
        with open(self._output_path / f"{tree_struct.name}.nw", "w") as f:
            f.write(newick)

    def _save_plot_tree(self, fig: Figure, tree_struct: TreeStruct):
        fig.suptitle(tree_struct.name, fontsize=16)
        ax = fig.gca()
        ax.axis("off")
        newick = tree_struct.tree.to_newick(include_distance=False)
        t = Phylo.read(StringIO(newick), "newick")
        t.ladderize()
        Phylo.draw(
            t,
            show_confidence=False,
            axes=ax,
            do_show=False,
            label_func=lambda clade: "" if not clade.name else tree_struct.distances.names[int(
                clade.name)],
            branch_labels=lambda clade: "" if not clade.name else "{:.2f}".format(
                clade.confidence) if clade.confidence else ""
        )
        pyplot.savefig(self._output_path / f"{tree_struct.name}.png")
        ax.clear()

    def save(self):
        fig = pyplot.figure(figsize=(12.0, 12.0))
        for tree_struct in self._trees:
            self._save_img_positions(tree_struct)
            self._save_align(tree_struct)
            self._save_newick_tree(tree_struct)
            self._save_distance_matrix(tree_struct)
            self._save_plot_tree(fig, tree_struct)
