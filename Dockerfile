FROM python:3.8.11-buster

ENV CONDA_ALWAYS_YES="true" \
    PATH="/root/miniconda3/bin:$PATH"

ADD requirements.txt .

RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh &&\
    sh Miniconda3-latest-Linux-x86_64.sh -b -p "/root/miniconda3" &&\
    conda config --add channels defaults &&\
    conda config --add channels bioconda &&\
    conda install -c bioconda clustalo &&\
    pip install -r requirements.txt