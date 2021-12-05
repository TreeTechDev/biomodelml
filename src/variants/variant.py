from Bio import SeqIO
from src.structs import DistanceStruct


class Variant:
    name = "variant"
    
    def __init__(self, fasta_file: str):
        seqs = SeqIO.parse(fasta_file, "fasta")
        self._names = []
        self._sequences = []
        for s in seqs:
            self._names.append(s.description)
            self._sequences.append(s.seq)


    def build_matrix(self) -> DistanceStruct:
        raise NotImplementedError()