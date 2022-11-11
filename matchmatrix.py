import sys
import shutil
from os import path, cpu_count, makedirs
from concurrent.futures import ThreadPoolExecutor
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
            if not path.exists(path.join(output_path, "full", f"{s.description}.png")):
                to_run.append(
                    (s.description, s.description, s.seq, s.seq, max_window, output_path)
                )
        print(f"starting to build image matrix for {len(to_run)} sequences")

        try:
            with ThreadPoolExecutor(max_workers=procs) as pool:
                [pool.submit(save_image_by_matrices, *data) for data in to_run]
        except:
            shutil.rmtree(outpath)

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    outpath = sys.argv[2]
    outpath = path.join(outpath, fasta_file.split(".")[0].split("/")[-1])
    makedirs(outpath, exist_ok=True)
    main(fasta_file, outpath)
