.PHONY: clean build push sanitize matches tree run optimize experiments cluster

.SECONDARY:

APP_DIR="/app"
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

push:
	docker push $(IMG_NAME)

run-docker:
	docker run $(DOCKER_FLAGS) -v $(FULL_ROOT_DIR):$(APP_DIR) $(IMG_NAME) $(CMD)

sanitize:
	CMD="python $(APP_DIR)/sanitize_seqs.py $(DATA_DIR)/$(SEQ).fasta $(TYPE)" $(MAKE) run-docker

matches:
	CMD="python $(APP_DIR)/matchmatrix.py $(DATA_DIR)/$(SEQ).fasta.$(TYPE).sanitized $(DATA_DIR)/images/$(TYPE)/ $(TYPE)" $(MAKE) run-docker

tree-by-channel:
	CMD="python $(APP_DIR)/tree_builder.py $(DATA_DIR)/$(SEQ).fasta.$(TYPE).sanitized $(DATA_DIR)/trees/$(TYPE)/$(CHANNEL)/ $(TYPE) $(DATA_DIR)/images/$(TYPE)/$(SEQ)/$(CHANNEL)/" $(MAKE) run-docker

tree:
	CHANNEL="full" $(MAKE) tree-by-channel
	-CHANNEL="red" $(MAKE) tree-by-channel
	-CHANNEL="green" $(MAKE) tree-by-channel
	-CHANNEL="blue" $(MAKE) tree-by-channel
	-CHANNEL="red_green" $(MAKE) tree-by-channel
	-CHANNEL="red_blue" $(MAKE) tree-by-channel
	-CHANNEL="green_blue" $(MAKE) tree-by-channel
	-CHANNEL="gray_r" $(MAKE) tree-by-channel
	-CHANNEL="gray_g" $(MAKE) tree-by-channel
	-CHANNEL="gray_b" $(MAKE) tree-by-channel
	-CHANNEL="gray_max" $(MAKE) tree-by-channel
	-CHANNEL="gray_mean" $(MAKE) tree-by-channel

run: | sanitize matches tree

optimize:
	CMD="python $(APP_DIR)/optimize.py $(DATA_DIR) $(SEQ)" $(MAKE) run-docker

exp_by_type:
	SEQ="orthologs_myoglobin" $(MAKE) run
	SEQ="orthologs_neuroglobin" $(MAKE) run
	SEQ="orthologs_hemoglobin_beta" $(MAKE) run
	SEQ="orthologs_cytoglobin" $(MAKE) run
	SEQ="orthologs_androglobin" $(MAKE) run
	SEQ="indelible" $(MAKE) run

experiments:
	# TYPE="P" $(MAKE) exp_by_type
	TYPE="N" $(MAKE) exp_by_type

cluster:
	CMD="bash run_blast.sh 11" DOCKER_FLAGS="-w $(APP_DIR)" $(MAKE) run-docker
	SEQ="orthologs_androglobin" TYPE="N" $(MAKE) matches
	SEQ="orthologs_cytoglobin" TYPE="N" $(MAKE) matches
	SEQ="orthologs_myoglobin" TYPE="N" $(MAKE) matches
	SEQ="orthologs_neuroglobin" TYPE="N" $(MAKE) matches
	SEQ="orthologs_hemoglobin_beta" TYPE="N" $(MAKE) matches
	CMD="python clusterize.py" DOCKER_FLAGS="-w $(APP_DIR)" $(MAKE) run-docker
