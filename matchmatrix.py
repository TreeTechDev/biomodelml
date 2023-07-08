import sys
import shutil
import os
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
from src.matrices import save_image_by_matrices


def main(fasta_file: str, output_path: str, seq_type: str):

    max_window = 255
    procs = os.cpu_count()

    with open(fasta_file, "r") as handle:
        sequences = SeqIO.parse(handle, "fasta")
        print(f"File with {len(list(sequences))} sequences")

    with open(fasta_file) as handle:
        sequences = SeqIO.parse(handle, "fasta")
        to_run = []
        for s in sequences:
            to_run.append(
                (s.description, s.description, s.seq, s.seq, max_window, output_path, seq_type)
            )
        print(f"starting to build image matrix for {len(to_run)} sequences")
        try:
            with ThreadPoolExecutor(max_workers=procs) as pool:
                futures = [pool.submit(save_image_by_matrices, *data) for data in to_run]
            [f.result() for f in futures]
        except Exception as e:
            shutil.rmtree(outpath)
            raise e

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    outpath = sys.argv[2]
    seq_type = sys.argv[3]
    outpath = os.path.join(outpath, fasta_file.split(".")[0].split("/")[-1])
    os.makedirs(outpath, exist_ok=True)
    main(fasta_file, outpath, seq_type)
