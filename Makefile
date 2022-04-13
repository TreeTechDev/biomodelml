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
	touch $(FULL_DATA_DIR)/trees/full/.keep
	touch $(FULL_DATA_DIR)/trees/red/.keep
	touch $(FULL_DATA_DIR)/trees/green/.keep
	touch $(FULL_DATA_DIR)/trees/blue/.keep

build:
	docker build . -t $(IMG_NAME)

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
	CMD="python $(APP_DIR)/matchmatrix.py $(DATA_DIR)/$(SEQ).fasta.sanitized $(DATA_DIR)/images/" $(MAKE) run-docker

tree:
	CMD="python $(APP_DIR)/tree_builder.py $(DATA_DIR)/$(SEQ).fasta.sanitized $(DATA_DIR)/trees/full $(TYPE) $(DATA_DIR)/images/$(SEQ)/full/" $(MAKE) run-docker
	CMD="python $(APP_DIR)/tree_builder.py $(DATA_DIR)/$(SEQ).fasta.sanitized $(DATA_DIR)/trees/red $(TYPE) $(DATA_DIR)/images/$(SEQ)/red/" $(MAKE) run-docker
	CMD="python $(APP_DIR)/tree_builder.py $(DATA_DIR)/$(SEQ).fasta.sanitized $(DATA_DIR)/trees/green $(TYPE) $(DATA_DIR)/images/$(SEQ)/green/" $(MAKE) run-docker
	CMD="python $(APP_DIR)/tree_builder.py $(DATA_DIR)/$(SEQ).fasta.sanitized $(DATA_DIR)/trees/blue $(TYPE) $(DATA_DIR)/images/$(SEQ)/blue/" $(MAKE) run-docker

run: | sanitize matches tree
