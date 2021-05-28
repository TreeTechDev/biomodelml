FROM python:3.8

ENV CONDA_ALWAYS_YES="true" \
    PATH="/root/miniconda3/bin:$PATH"

RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh &&\
    sh Miniconda3-latest-Linux-x86_64.sh -b -p "/root/miniconda3" &&\
    conda config --add channels defaults &&\
    conda config --add channels bioconda &&\
    conda config --add channels conda-forge &&\
    conda install -c bioconda fastp bowtie2 sra-tools fastqc trimmomatic picard samtools bedtools cd-hit emboss
