.PHONY: clean build sanitize matches tree run

.SECONDARY:

APP_DIR="/app"
DATA_DIR=$(APP_DIR)/data
FULL_ROOT_DIR=`pwd`
FULL_DATA_DIR=$(FULL_ROOT_DIR)/data

IMG_NAME="bioinfo2.0"

clean:
	rm -rf $(FULL_DATA_DIR)/*
	mkdir -p $(FULL_DATA_DIR)/images
	mkdir -p $(FULL_DATA_DIR)/trees
	touch $(FULL_DATA_DIR)/.keep
	touch $(FULL_DATA_DIR)/images/.keep
	touch $(FULL_DATA_DIR)/trees/.keep

build:
	docker build . -t $(IMG_NAME)

run-docker:
	docker run -it -v $(FULL_ROOT_DIR):$(APP_DIR) $(IMG_NAME) $(CMD)

evol:
	CMD="ete3 evol -t '$(DATA_DIR)/trees/orthologs_cytoglobin/Control with Clustal Omega.nw' --alg '$(DATA_DIR)/trees/orthologs_cytoglobin/Control with Clustal Omega.fasta' --models M2 M1 b_free b_neut --leaves --tests b_free,b_neut --cpu 4" $(MAKE) run-docker

compare:
	CMD="ete3 compare -t '$(DATA_DIR)/trees/orthologs_cytoglobin/Control with Clustal Omega.nw' -r '$(DATA_DIR)/trees/orthologs_cytoglobin/MultiScale Structural Similarity Index Measure.nw' --unrooted" $(MAKE) run-docker

# example: SEQ=MIDORI_LONGEST_NUC_GB246_A6_RAW TYPE=N make sanitize
sanitize:
	CMD="python $(APP_DIR)/sanitize_seqs.py $(DATA_DIR)/$(SEQ).fasta $(TYPE)" $(MAKE) run-docker

matches:
	CMD="python $(APP_DIR)/matchmatrix.py $(DATA_DIR)/$(SEQ).fasta.sanitized $(DATA_DIR)/images/" $(MAKE) run-docker

tree:
	CMD="python $(APP_DIR)/tree_builder.py $(DATA_DIR)/$(SEQ).fasta.sanitized $(DATA_DIR)/trees/ $(TYPE) $(DATA_DIR)/images/$(SEQ)/" $(MAKE) run-docker

run: | sanitize matches tree