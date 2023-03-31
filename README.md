# bioinfo2.0

In this repository all the results of the Bioinfo2.0 research are gathered, with a special focus on the use of self-comparison matrices.

## Steps to Run

Each of the versions (v1, v2, v2-v3, v3, v4, v5) is on Github in the Bioinfo2.0 repository which can be downloaded on the computer by the command:

    git clone https://github.com/BioBD/bioinfo2.0.git

Once downloaded to your computer, you need Docker, Linux or Mac with GNU programs like Make to reproduce the experiments and Python 3.7+ if you want to run the Notebooks. There is a Makefile with all the commands in each version to reproduce the experiments inside a Docker. To change between versions just run the command:

    git checkout branch-X

Where *X* is the version.

Before doing anything when changing versions you need to rebuild the Docker image with the command:

    make build

There are only significant changes to the Dockerfile between v1 and v2, so no need to re-image switching between newer versions. Another command that already brings the results and data obtained is:

    dvc pull

The script to get the datasets used by Globinas is executed by:

    ./get_globins.sh

## Changes between Versions of Algorithms

### v1

Initial version that generated a first version of self-comparison matrices still quite limited with numbers for exact matches and with the reverse complement of the sequence. A clustering method was used based on the result of the similarity between the matrices, but without much success and with a high computational cost.

### v2

First version to use the concept of images to represent self-comparison data in channels. For each channel R, G, B were used:

- **R**: the sequence against itself
- **G**: the sequence against the reverse complement
- **B**: the sequence against its inverse

Comparison between images of different sizes is done using resize and the best algorithm was MS-SSIM for the task, called Resized MS-SSIM (RMS-SSIM).

### v3

Second version to use the concept of images to represent the self-comparison data in channels. For each channel R, G, B were used:

- **R**: the sequence against itself, including inverses
- **G**: the sequence against the reverse complement
- **B**: where they had no match in any of the previous ones

In this version, it was evaluated that showing what is not a match gives good results and since the removal of the last channel had no impact on anything, it was decided to make it clear to the image comparison algorithms where the differences were and there were no matches for this to be taken into account when evaluating the similarities between the self-comparison images.

This matrix was used until the end of the work to generate features for the image algorithms, having good results with MS-SSIM and SSIM.

Another novelty in this version is the algorithm used through MS-SSIM to place a small image on top of a smaller image walking diagonally until finding the best match, called Windowed MS-SSIM (WMS-SSIM).

### v2-v3
Same algorithm RMS-SSIM as v2, but with newer RGB matrix. 

### v4

The novelty here is the diagonally walking algorithm, but now using SSIM given some technical limitations of MS-SSIM. The idea is to go through the diagonal until finding the best match of a first slice of the small image in the larger image and from that best match to look for the best matches for each new slice, called Greedy Sliced SSIM (GS-SSIM).

### v5

A new version of SSIM, but now without the restriction of looking for the best matches from a previous one. Here the algorithm finds the best matches by looking for slices of the smaller image across the larger image diagonally. The final score is calculated by averaging all the best matches, called Unrestricted Sliced SSIM (US-SSIM).

## Running Experiments

Before run the experiments you should go to the right version branch of this version with:

    git checkout branch-X  # where X is the version
    make build  # just if change from v1 to others
    dvc pull

The experiments described here are just valid for v2, v2-v3, v3, v4 and v5.

### Phylogenetic Dendrograms

This command generate results for dendrograms analysis and similarity metrics by homologues:

    make clean experiments

### Homologues Search and Clusterize 

This command generate results using algoritms to search into multiple sequences for each sequence looking for most similar sequences as homologues:

    make clean cluster

## Algoritms Optimization

There are only algoritms of versions v2-v3, v3, v4 and v5 that were implemented bayesian optimization in its hyperparameters. Before run the optimizations you should run:

    git checkout branch-X  # where X is the version
    make build  # just if change from v1 to others
    dvc pull
    make clean

Each optimization runs for each fasta file and can be one of these, where SEQ is fasta file name:

    SEQ="orthologs_hemoglobin_beta" TYPE="N" make optimize 
    SEQ="orthologs_myoglobin" TYPE="N" make optimize
    SEQ="orthologs_neuroglobin" TYPE="N" make optimize
    SEQ="orthologs_cytoglobin" TYPE="N" make optimize 