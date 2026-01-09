FROM python:3.10-bookworm

ENV CONDA_ALWAYS_YES="true" \
    PATH="/opt/miniforge3/bin:/usr/local/bin:$PATH" \
    LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/miniforge3/lib/"

ADD requirements.txt .
ADD requirements_test.txt .

# Instala dependências de compilação
RUN apt-get update && apt-get install -y \
    libgl1-mesa-glx \
    build-essential \
    autoconf \
    automake \
    pkg-config \
    git \
    libopenblas-dev \
    liblapack-dev \
    liblapacke-dev \
    && rm -rf /var/lib/apt/lists/*

# Compila SLR do código fonte
RUN git clone https://github.com/timmassingham/SLR.git /tmp/slr && \
    cd /tmp/slr && \
    make -f Makefile && \
    cp /tmp/slr/bin/Slr /usr/local/bin/slr && \
    rm -rf /tmp/slr

# Compila PAML do código fonte
RUN git clone https://github.com/abacus-gene/paml.git /tmp/paml && \
    cd /tmp/paml/src && \
    make -j4 && \
    cp /tmp/paml/src/baseml /tmp/paml/src/codeml /tmp/paml/src/evolver \
       /tmp/paml/src/mcmctree /tmp/paml/src/pamp /tmp/paml/src/yn00 /usr/local/bin/ && \
    rm -rf /tmp/paml

# Compila PhyML do código fonte
RUN git clone https://github.com/stephaneguindon/phyml.git /tmp/phyml && \
    cd /tmp/phyml && \
    sh ./autogen.sh && \
    ./configure --enable-phyml && \
    make -j4 && \
    cp /tmp/phyml/src/phyml /usr/local/bin/ && \
    rm -rf /tmp/phyml

# Instala Miniforge (conda-forge)
RUN ARCH=$(uname -m) && \
    if [ "$ARCH" = "aarch64" ] || [ "$ARCH" = "arm64" ]; then \
        MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh"; \
    else \
        MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"; \
    fi && \
    curl -L -O $MINIFORGE_URL && \
    sh Miniforge3-Linux-*.sh -b -p "/opt/miniforge3" && \
    rm Miniforge3-Linux-*.sh && \
    conda config --add channels bioconda && \
    conda install python=3.10 -y && \
    pip install -r requirements.txt && \
    pip install -r requirements_test.txt && \
    pip install lxml six && \
    conda install pyqt clustalo muscle iqtree blast -y && \
    conda clean -afy
