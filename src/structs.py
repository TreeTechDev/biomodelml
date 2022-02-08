import numpy
from dataclasses import dataclass
from typing import List, Optional
from biotite.sequence.phylo import Tree
from biotite.sequence.align import Alignment


@dataclass
class DistanceStruct:
    names: List[str]
    matrix: numpy.ndarray
    align: Optional[Alignment] = None


@dataclass
class TreeStruct:
    name: str
    distances: DistanceStruct
    tree: Tree