from Bio import SeqIO
from src.structs import DistanceStruct


class Variant:
    name = "variant"
    protein_type = "P"
    nucleotide_type = "N"
    def __init__(self, fasta_file: str, sequence_type):
        if sequence_type not in [self.protein_type, self.nucleotide_type]:
            raise IOError(
                f"Sequence must be a protein or nucleotide and is: {sequence_type}")
        seqs = SeqIO.parse(fasta_file, "fasta")
        self._names = []
        self._sequences = []
        self._fasta_file = fasta_file
        self._sequence_type = sequence_type
        for s in seqs:
            self._names.append(s.description)
            self._sequences.append(s.seq)


    def build_matrix(self) -> DistanceStruct:
        raise NotImplementedError()