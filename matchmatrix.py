import os
import sys
from multiprocessing import Pool, cpu_count
from Bio import SeqIO
from src.matrices import save_image_by_matrices


def main(fasta_file: str, output_path: str):

    max_window = 255
    procs = cpu_count()

    with open(fasta_file, "r") as handle:
        sequences = SeqIO.parse(handle, "fasta")
        print(f"File with {len(list(sequences))} sequences")

    with open(fasta_file) as handle:
        sequences = SeqIO.parse(handle, "fasta")
        to_run = []
        for s in sequences:
            if not os.path.exists(os.path.join(output_path, "full", f"{s.description}.png")):
                to_run.append(
                    (s.description, s.description, s.seq, s.seq, max_window, output_path)
                )
        print(f"starting to build image matrix for {len(to_run)} sequences")

        with Pool(procs) as pool:
            pool.starmap(
                save_image_by_matrices,
                to_run
            )


if __name__ == "__main__":
    fasta_file = sys.argv[1]
    outpath = sys.argv[2]
    outpath = os.path.join(outpath, fasta_file.split(".")[0].split("/")[-1])
    os.makedirs(outpath, exist_ok=True)
    main(fasta_file, outpath)
