FROM python:3.8.11-buster

ENV CONDA_ALWAYS_YES="true" \
    PATH="/root/miniconda3/bin:$PATH"

RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh &&\
    sh Miniconda3-latest-Linux-x86_64.sh -b -p "/root/miniconda3" &&\
    conda config --add channels defaults &&\
    conda config --add channels bioconda &&\
    conda config --add channels conda-forge &&\
    conda install -c bioconda fastp bowtie2 sra-tools fastqc trimmomatic picard samtools bedtools cd-hit emboss muscle &&\
    conda install -c conda-forge libiconv &&\
    apt update && apt -y install multiarch-support &&\
    wget http://security.debian.org/debian-security/pool/updates/main/o/openssl/libssl1.0.0_1.0.1t-1+deb8u12_amd64.deb &&\
    dpkg -i libssl1.0.0_1.0.1t-1+deb8u12_amd64.deb