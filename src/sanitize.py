import os
from Bio import SeqIO
from biotite.sequence import NucleotideSequence, ProteinSequence, AlphabetError
from src.structs import SeqTypeStruct


NUCLEOTIDE_SYMBOLS = NucleotideSequence.alphabet_unamb.get_symbols() + ["U"]
PROTEIN_SYMBOLS = ProteinSequence.alphabet.get_symbols()[:20]
SEQ_TYPES = SeqTypeStruct(N=NUCLEOTIDE_SYMBOLS, P=PROTEIN_SYMBOLS)


def convert_and_remove_unrelated_sequences(seq_path: str, seq_type):
    if seq_type not in SEQ_TYPES.__dataclass_fields__:
        raise IOError("Only sequence type N or P accepted")
    with open(seq_path) as handle:
        sequences = SeqIO.parse(handle, "fasta")
        sanitized_seqs = []
        for s in sequences:
            alphabet = set(s.seq)
            if alphabet.issubset(getattr(SEQ_TYPES, seq_type)):
                if seq_type == "P":
                    try:
                        NucleotideSequence(s.seq, False)
                        t = s.translate(stop_symbol="")
                        t.description = s.description
                        t.id = s.id
                        print(f"Sequence {t.description} translated")
                    except AlphabetError:
                        pass
                else:
                    t = s
                sanitized_seqs.append(t)
            else:
                print(f"Sequence {s.description} removed")

    print(f"writing {len(sanitized_seqs)} sequences")
    SeqIO.write(sanitized_seqs, f"{seq_path}.{seq_type}.sanitized", "fasta")