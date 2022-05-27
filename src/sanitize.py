from Bio import SeqIO
from biotite.sequence import NucleotideSequence, ProteinSequence


NUCLEOTIDE_SYMBOLS = NucleotideSequence.alphabet_unamb.get_symbols() + ["U"]
PROTEIN_SYMBOLS = ProteinSequence.alphabet.get_symbols()[:20]

SEQ_TYPES = {
    "N": NUCLEOTIDE_SYMBOLS,
    "P": PROTEIN_SYMBOLS
}

def remove_unrelated_sequences(seq_path: str, seq_type):

    with open(seq_path) as handle:
        sequences = SeqIO.parse(handle, "fasta")
        sanitized_seqs = []
        for s in sequences:
            alphabet = set(s.seq)
            if alphabet.issubset(SEQ_TYPES[seq_type]):
                sanitized_seqs.append(s)

    print(f"writing {len(sanitized_seqs)} sequences")
    SeqIO.write(sanitized_seqs, seq_path + ".sanitized", "fasta")