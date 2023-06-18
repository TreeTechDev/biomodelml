import sys
from Bio import SeqIO
from biotite.sequence import NucleotideSequence, ProteinSequence, AlphabetError

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

    print(f"Writing {len(sanitized_seqs)} sequences")
    SeqIO.write(sanitized_seqs, f"{seq_path}.{seq_type}.sanitized", "fasta")


if __name__ == "__main__":
    seq_path = sys.argv[1]
    seq_type = sys.argv[2]
    main(seq_path, seq_type)