FROM python:3.8.11-buster

ENV CONDA_ALWAYS_YES="true" \
    PATH="/opt/miniconda3/bin:$PATH" \
    LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/miniconda3/lib/"

ADD requirements.txt .
ADD requirements_test.txt .

RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh &&\
    sh Miniconda3-latest-Linux-x86_64.sh -b -p "/opt/miniconda3" &&\
    conda config --add channels bioconda &&\
    conda config --add channels etetoolkit &&\
    apt-get update && apt-get install libgl1-mesa-glx -y &&\
    pip install -r requirements.txt &&\
    pip install -r requirements_test.txt &&\
    pip install pyqt5 lxml six &&\
    conda install -c bioconda -c etetoolkit slr clustalo paml phyml muscle iqtree cudnn cudatoolkit &&\
    ln -s /opt/miniconda3/bin/ete3_apps/bin/Slr /opt/miniconda3/bin/Slr &&\
    ln -s /opt/miniconda3/bin/ete3_apps/bin/phyml /opt/miniconda3/bin/phyml &&\
    ete3 build check &&\
    rm -rf /var/lib/apt/lists/* &&\
    conda clean -afy
