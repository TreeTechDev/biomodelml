import sys
from Bio import SeqIO
from biotite.sequence import NucleotideSequence, ProteinSequence

nucleotide_symbols = NucleotideSequence.alphabet_unamb.get_symbols() + ["U"]
protein_symbols = ProteinSequence.alphabet.get_symbols()[:20]

seq_types = {
    "N": nucleotide_symbols,
    "P": protein_symbols
}

def main(seq_path: str, seq_type):


    with open(seq_path) as handle:
        sequences = SeqIO.parse(handle, "fasta")
        sanitized_seqs = []
        for s in sequences:
            alphabet = set(s.seq)
            if alphabet.issubset(seq_types[seq_type]):
                sanitized_seqs.append(s)

    print(f"writing {len(sanitized_seqs)} sequences")
    SeqIO.write(sanitized_seqs, seq_path + ".sanitized", "fasta")


if __name__ == "__main__":
    seq_path = sys.argv[1]
    seq_type = sys.argv[2]
    main(seq_path, seq_type)