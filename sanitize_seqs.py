import sys
from src.sanitize import remove_unrelated_sequences


if __name__ == "__main__":
    seq_path = sys.argv[1]
    seq_type = sys.argv[2]
    remove_unrelated_sequences(seq_path, seq_type)