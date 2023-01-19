.PHONY: clean build test pull push sanitize matches tree run experiments optimize try

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

build:
	docker build . -t $(IMG_NAME)
	pip install 'dvc==2.34.0' 'dvc-gdrive==2.19.0'

test:
	DOCKER_FLAGS="-it" CMD="pytest $(APP_DIR)" $(MAKE) run-docker

pull:
	docker pull $(IMG_NAME)

push:
	docker push $(IMG_NAME)

run-docker:
	docker run $(DOCKER_FLAGS) -v $(FULL_ROOT_DIR):$(APP_DIR) $(IMG_NAME) $(CMD)
	sudo chown -R $(USER) $(FULL_ROOT_DIR)

sanitize:
	CMD="python $(APP_DIR)/sanitize_seqs.py $(DATA_DIR)/$(SEQ).fasta $(TYPE)" $(MAKE) run-docker

matches:
	CMD="python $(APP_DIR)/matchmatrix.py $(DATA_DIR)/$(SEQ).fasta.sanitized $(DATA_DIR)/images/" $(MAKE) run-docker

tree-by-channel:
	CMD="python $(APP_DIR)/tree_builder.py $(DATA_DIR)/$(SEQ).fasta.sanitized $(DATA_DIR)/trees/$(CHANNEL) $(TYPE) $(DATA_DIR)/images/$(SEQ)/$(CHANNEL)/" $(MAKE) run-docker

t_%:
	CHANNEL="$*" $(MAKE) tree-by-channel

tree: t_full t_gray_r t_gray_b t_gray_g t_gray_mean

validate:
	CMD="python $(APP_DIR)/validate.py $(DATA_DIR)/trees/ $(SEQ)" $(MAKE) run-docker

run: | sanitize matches tree validate

e_%: 
	SEQ="orthologs_$*" TYPE="N" $(MAKE) run

experiments: e_hemoglobin_beta e_myoglobin e_neuroglobin e_cytoglobin #e_androglobin

optimize: | sanitize matches
	CMD="python $(APP_DIR)/optimize.py $(DATA_DIR) $(SEQ)" $(MAKE) run-docker

try:
	rm -rf $(FULL_DATA_DIR)/images/orthologs_neuroglobin/*
	SEQ="orthologs_neuroglobin" TYPE="N"  CHANNEL="full" $(MAKE) sanitize matches tree-by-channel
