import sys
from src.sanitize import convert_and_remove_unrelated_sequences


if __name__ == "__main__":
    seq_path = sys.argv[1]
    seq_type = sys.argv[2]
    convert_and_remove_unrelated_sequences(seq_path, seq_type)