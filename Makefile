.PHONY: clean build test pull push sanitize matches tree run experiments afproject optimize try cluster

.SECONDARY:

NPROCS = $(shell grep -c 'processor' /proc/cpuinfo)
MAKEFLAGS += -j1
APP_DIR=/app
DATA_DIR=$(APP_DIR)/data
FULL_ROOT_DIR=`pwd`
FULL_DATA_DIR=$(FULL_ROOT_DIR)/data

IMG_NAME="dmvieira/bioinfo2.0"

clean:
	mkdir -p $(FULL_DATA_DIR)/images
	mkdir -p $(FULL_DATA_DIR)/trees
	rm -rf $(FULL_DATA_DIR)/images/*
	rm -rf $(FULL_DATA_DIR)/trees/*
	touch $(FULL_DATA_DIR)/.keep
	touch $(FULL_DATA_DIR)/images/.keep
	touch $(FULL_DATA_DIR)/trees/.keep
	rm $(FULL_DATA_DIR)/*.pkl
	rm $(FULL_DATA_DIR)/final_cluster.csv
	rm $(FULL_DATA_DIR)/*.db

build:
	docker build . -t $(IMG_NAME)
	pip install 'dvc<3.0.0' 'dvc-gdrive<3.0.0'

test:
	DOCKER_FLAGS="-it" CMD="pytest $(APP_DIR)" $(MAKE) run-docker

pull:
	docker pull $(IMG_NAME)

push:
	docker push $(IMG_NAME)

run-docker:
	docker run -it $(DOCKER_FLAGS) -v $(FULL_ROOT_DIR):$(APP_DIR) $(IMG_NAME) $(CMD)

sanitize:
	CMD="python $(APP_DIR)/sanitize_seqs.py $(DATA_DIR)/$(SEQ).fasta $(TYPE)" $(MAKE) run-docker

matches:
	CMD="python $(APP_DIR)/matchmatrix.py $(DATA_DIR)/$(SEQ).fasta.$(TYPE).sanitized $(DATA_DIR)/images/$(TYPE)/ $(TYPE)" $(MAKE) run-docker

tree-by-channel:
	-CMD="python $(APP_DIR)/tree_builder.py $(DATA_DIR)/$(SEQ).fasta.$(TYPE).sanitized $(DATA_DIR)/trees/$(TYPE)/$(CHANNEL)/ $(TYPE) $(DATA_DIR)/images/$(TYPE)/$(SEQ)/$(CHANNEL)/" $(MAKE) run-docker

t_%:
	CHANNEL="$*" $(MAKE) tree-by-channel

tree: t_full t_gray_r t_gray_g t_gray_b

run: | sanitize matches tree

optimize: | sanitize matches
	CMD="python $(APP_DIR)/optimize.py $(DATA_DIR) $(SEQ)" $(MAKE) run-docker

try:
	rm -rf $(FULL_DATA_DIR)/images/N/orthologs_neuroglobin/*
	rm -rf $(FULL_DATA_DIR)/trees/N/full/orthologs_neuroglobin/*
	SEQ="orthologs_neuroglobin" TYPE="N"  CHANNEL="full" $(MAKE) sanitize matches tree-by-channel

exp_by_type:
	SEQ="orthologs_hemoglobin_beta" $(MAKE) run
	SEQ="orthologs_myoglobin" $(MAKE) run
	SEQ="orthologs_neuroglobin" $(MAKE) run
	SEQ="orthologs_cytoglobin" $(MAKE) run
	SEQ="orthologs_androglobin" $(MAKE) run
	SEQ="indelible" $(MAKE) run

experiments:
	TYPE="P" $(MAKE) exp_by_type
	TYPE="N" $(MAKE) exp_by_type

afproject:
	TYPE="P" SEQ="ST001" $(MAKE) run
	TYPE="P" SEQ="ST002" $(MAKE) run
	TYPE="P" SEQ="ST003" $(MAKE) run
	TYPE="P" SEQ="ST004" $(MAKE) run
	TYPE="P" SEQ="ST005" $(MAKE) run
	TYPE="P" SEQ="ST007" $(MAKE) run
	TYPE="P" SEQ="ST008" $(MAKE) run
	TYPE="P" SEQ="ST009" $(MAKE) run
	TYPE="P" SEQ="ST010" $(MAKE) run
	TYPE="P" SEQ="ST011" $(MAKE) run
	TYPE="P" SEQ="ST012" $(MAKE) run
	TYPE="N" SEQ="fish_mito" $(MAKE) run


cluster:
	CMD="bash run_blastn.sh 11" DOCKER_FLAGS="-w $(APP_DIR)" $(MAKE) run-docker
	CMD="bash run_blastp.sh 4" DOCKER_FLAGS="-w $(APP_DIR)" $(MAKE) run-docker
	SEQ="orthologs_androglobin" TYPE="P" $(MAKE) sanitize matches
	SEQ="orthologs_cytoglobin" TYPE="P" $(MAKE) sanitize matches
	SEQ="orthologs_myoglobin" TYPE="P" $(MAKE) sanitize matches
	SEQ="orthologs_neuroglobin" TYPE="P" $(MAKE) sanitize matches
	SEQ="orthologs_hemoglobin_beta" TYPE="P" $(MAKE) sanitize matches
	SEQ="indelible" TYPE="P" $(MAKE) sanitize matches
	SEQ="orthologs_androglobin" TYPE="N" $(MAKE) sanitize matches
	SEQ="orthologs_cytoglobin" TYPE="N" $(MAKE) sanitize matches
	SEQ="orthologs_myoglobin" TYPE="N" $(MAKE) sanitize matches
	SEQ="orthologs_neuroglobin" TYPE="N" $(MAKE) sanitize matches
	SEQ="orthologs_hemoglobin_beta" TYPE="N" $(MAKE) sanitize matches
	SEQ="indelible" TYPE="N" $(MAKE) sanitize matches
	CMD="python clusterize.py" DOCKER_FLAGS="-it -w $(APP_DIR)" $(MAKE) run-docker
