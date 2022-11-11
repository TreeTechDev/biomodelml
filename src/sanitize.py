import os
from Bio import SeqIO
from biotite.sequence import NucleotideSequence, ProteinSequence
from src.structs import SeqTypeStruct


NUCLEOTIDE_SYMBOLS = NucleotideSequence.alphabet_unamb.get_symbols() + ["U"]
PROTEIN_SYMBOLS = ProteinSequence.alphabet.get_symbols()[:20]
SEQ_TYPES = SeqTypeStruct(N=NUCLEOTIDE_SYMBOLS, P=PROTEIN_SYMBOLS)


def remove_unrelated_sequences(seq_path: str, seq_type):
    if seq_type not in SEQ_TYPES.__dataclass_fields__:
        raise IOError("Only sequence type N or P accepted")
    if os.path.exists(seq_path + ".sanitized"):
        print(f"already sanitized file {seq_path + '.sanitized'}")
        return
    with open(seq_path) as handle:
        sequences = SeqIO.parse(handle, "fasta")
        sanitized_seqs = []
        for s in sequences:
            alphabet = set(s.seq)
            if alphabet.issubset(getattr(SEQ_TYPES, seq_type)):
                sanitized_seqs.append(s)

    print(f"writing {len(sanitized_seqs)} sequences")
    SeqIO.write(sanitized_seqs, seq_path + ".sanitized", "fasta")