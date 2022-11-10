FROM python:3.8.11-buster

ENV CONDA_ALWAYS_YES="true" \
    PATH="/root/miniconda3/bin:$PATH"

ADD requirements.txt .
ADD requirements_test.txt .

RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh &&\
    sh Miniconda3-latest-Linux-x86_64.sh -b -p "/root/miniconda3" &&\
    conda config --add channels bioconda &&\
    conda config --add channels etetoolkit &&\
    apt-get update && apt-get install libgl1-mesa-glx -y &&\
    pip install -r requirements.txt &&\
    pip install -r requirements_test.txt &&\
    pip install pyqt5 lxml six &&\
    pip install --upgrade ete3 &&\
    conda install -c bioconda -c etetoolkit slr clustalo paml phyml muscle iqtree &&\
    ln -s /root/miniconda3/bin/ete3_apps/bin/Slr /root/miniconda3/bin/Slr &&\
    ln -s /root/miniconda3/bin/ete3_apps/bin/phyml /root/miniconda3/bin/phyml &&\
    ete3 build check &&\
    rm -rf /var/lib/apt/lists/* &&\
    conda clean -afy

