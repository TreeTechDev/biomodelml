.PHONY: clean build test pull push sanitize matches tree run experiments try

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

build:
	docker build . -t $(IMG_NAME)

test:
	CMD="pytest $(APP_DIR)" $(MAKE) run-docker

pull:
	docker pull $(IMG_NAME)

push:
	docker push $(IMG_NAME)

run-docker:
	docker run -it -v $(FULL_ROOT_DIR):$(APP_DIR) $(IMG_NAME) $(CMD)

iqtree:
	CMD="iqtree -s '$(DATA_DIR)/trees/orthologs_cytoglobin/Control with Clustal Omega.fasta' -t '$(DATA_DIR)/trees/orthologs_cytoglobin/Control with Clustal Omega.nw'" $(MAKE) run-docker

evol:
	CMD="ete3 evol -t '$(DATA_DIR)/trees/orthologs_cytoglobin/MultiScale Structural Similarity Index Measure.nw' --alg '$(DATA_DIR)/trees/orthologs_cytoglobin/Control with Clustal Omega.fasta' --models M2 M1 b_free b_neut --leaves --tests b_free,b_neut --cpu 4" $(MAKE) run-docker

compare:
	CMD="ete3 compare -t '$(DATA_DIR)/trees/orthologs_cytoglobin/MultiScale Structural Similarity Index Measure.nw' -r '$(DATA_DIR)/trees/orthologs_cytoglobin/MultiScale Structural Similarity Index Measure.nw' --unrooted" $(MAKE) run-docker

view:
	CMD="ete3 view -t '$(DATA_DIR)/trees/orthologs_cytoglobin/Control with Clustal Omega.nw'" $(MAKE) run-docker

# example: SEQ=MIDORI_LONGEST_NUC_GB246_A6_RAW TYPE=N make sanitize
sanitize:
	CMD="python $(APP_DIR)/sanitize_seqs.py $(DATA_DIR)/$(SEQ).fasta $(TYPE)" $(MAKE) run-docker

matches:
	@CMD="python $(APP_DIR)/matchmatrix.py $(DATA_DIR)/$(SEQ).fasta.sanitized $(DATA_DIR)/images/" $(MAKE) run-docker

tree-by-channel:
	CMD="python $(APP_DIR)/tree_builder.py $(DATA_DIR)/$(SEQ).fasta.sanitized $(DATA_DIR)/trees/$(CHANNEL) $(TYPE) $(DATA_DIR)/images/$(SEQ)/$(CHANNEL)/" $(MAKE) run-docker

tree:
	CHANNEL="full" $(MAKE) tree-by-channel
	CHANNEL="red" $(MAKE) tree-by-channel
	CHANNEL="green" $(MAKE) tree-by-channel
	CHANNEL="blue" $(MAKE) tree-by-channel
	CHANNEL="red_green" $(MAKE) tree-by-channel
	CHANNEL="red_blue" $(MAKE) tree-by-channel
	CHANNEL="green_blue" $(MAKE) tree-by-channel
	CHANNEL="gray_r" $(MAKE) tree-by-channel
	CHANNEL="gray_g" $(MAKE) tree-by-channel
	CHANNEL="gray_b" $(MAKE) tree-by-channel
	CHANNEL="gray_max" $(MAKE) tree-by-channel
	CHANNEL="gray_mean" $(MAKE) tree-by-channel

validate:
	CMD="python $(APP_DIR)/validate.py $(DATA_DIR)/trees/ $(SEQ)" $(MAKE) run-docker

run: | pull sanitize matches tree validate

experiments:
	SEQ="orthologs_hemoglobin_beta" TYPE="N" $(MAKE) run
	SEQ="orthologs_myoglobin" TYPE="N" $(MAKE) run
	SEQ="orthologs_neuroglobin" TYPE="N" $(MAKE) run
	SEQ="orthologs_cytoglobin" TYPE="N" $(MAKE) run
	SEQ="orthologs_androglobin" TYPE="N" $(MAKE) run
	SEQ="indelible" TYPE="N" $(MAKE) run

try:
	rm -rf $(FULL_DATA_DIR)/images/orthologs_neuroglobin/*
	SEQ="orthologs_neuroglobin" TYPE="N"  CHANNEL="full" $(MAKE) sanitize matches tree-by-channel validate