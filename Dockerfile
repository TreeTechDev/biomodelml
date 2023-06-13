FROM python:3.8.11-buster

ENV CONDA_ALWAYS_YES="true" \
    PATH="/opt/miniconda3/bin:$PATH" \
    LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/miniconda3/lib/"

ADD requirements.txt .

RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh &&\
    sh Miniconda3-latest-Linux-x86_64.sh -b -p "/opt/miniconda3" &&\
    conda config --add channels bioconda &&\
    apt-get update && apt-get install libgl1-mesa-glx -y &&\
    pip install -r requirements.txt &&\
    pip install pyqt5 lxml six &&\
    conda install -c bioconda clustalo cudnn cudatoolkit blast &&\
    rm -rf /var/lib/apt/lists/* &&\
    conda clean -afy
