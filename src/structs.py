import numpy
from dataclasses import dataclass
from typing import List
from biotite.sequence.phylo import Tree


@dataclass
class DistanceStruct:
    names: List[str]
    matrix: numpy.ndarray


@dataclass
class TreeStruct:
    name: str
    names: List[str]
    tree: Tree