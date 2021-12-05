from src.variants.variant import Variant
from src.structs import DistanceStruct
from biotite import sequence
from biotite.application import clustalo
from biotite.sequence import align


class ControlVariant(Variant):

    name = "Control with Clustal Omega"

    def __init__(self, fasta_file):
        super().__init__(fasta_file)
        self._sequences = list(map(sequence.NucleotideSequence, self._sequences))

    def build_matrix(self) -> DistanceStruct:
        alignment = clustalo.ClustalOmegaApp.align(self._sequences)
        distances = 1 - align.get_pairwise_sequence_identity(
            alignment, mode="shortest"
        )
        return DistanceStruct(names=self._names, matrix=distances)
