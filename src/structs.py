import numpy
from dataclasses import dataclass
from typing import List, Optional, Any
from biotite.sequence.phylo import Tree
from biotite.sequence.align import Alignment

@dataclass
class ImgMap:
    positions: List[int]
    scores: List[float]

@dataclass
class DistanceStruct:
    names: List[str]
    matrix: numpy.ndarray
    align: Optional[Alignment] = None
    img_positions: Optional[List[List[Any]]] = None


@dataclass
class TreeStruct:
    name: str
    distances: DistanceStruct
    tree: Tree


@dataclass
class SeqTypeStruct:
    N: List[str]
    P: List[str]