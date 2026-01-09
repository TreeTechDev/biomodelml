import numpy
from src.variants.control import ControlVariant
from src.structs import DistanceStruct
from biotite.sequence import align, Sequence


class SmithWatermanVariant(ControlVariant):

    name = "Local with Smith–Waterman"


    def _get_substitution_matrix(self) -> align.SubstitutionMatrix:
        if self._sequence_type == self.nucleotide_type:
            return align.SubstitutionMatrix.std_nucleotide_matrix()
        return align.SubstitutionMatrix.std_protein_matrix()

    def call_alg(self, seq: Sequence, seq_compare: Sequence) -> float:
        matrix = self._get_substitution_matrix()
        max_score = matrix.score_matrix().max()*max(len(seq), len(seq_compare))
        return 1 - align.align_optimal(
            seq, seq_compare, matrix, local=True)[0].score / max_score

    def build_matrix(self) -> DistanceStruct:
        distances = numpy.empty(
            (len(self._sequences), len(self._sequences)), numpy.float64)
        for i, seq in enumerate(self._sequences):
            for j, seq_compare in enumerate(self._sequences[i:]):
                score = max(self.call_alg(seq, seq_compare), 0)
                distances[i, j+i] = score
                distances[j+i, i] = score
        # Ensure strictly positive distances
        epsilon = 1e-8
        distances = numpy.clip(distances, epsilon, None)
        numpy.fill_diagonal(distances, 0.0)
        distances = (distances + distances.T) / 2
        return DistanceStruct(
            names=self._names, matrix=distances)
