.PHONY: clean build workflow reports

.SECONDARY:

DATA_DIR="/data"
TRIMM_DIR=$(DATA_DIR)/trimmomatic/adapters
FULL_DATA_DIR=`pwd`$(DATA_DIR)

IMG_NAME="bioinfo2.0"

clean:
	rm -rf $(FULL_DATA_DIR)/*
	mkdir -p $(FULL_DATA_DIR)/matrix
	mkdir -p $(FULL_DATA_DIR)/images
	mkdir -p $(FULL_DATA_DIR)/trees
	touch $(FULL_DATA_DIR)/.keep
	touch $(FULL_DATA_DIR)/matrix/.keep
	touch $(FULL_DATA_DIR)/images/.keep
	touch $(FULL_DATA_DIR)/trees/.keep

build:
	docker build . -t $(IMG_NAME)

run-docker:
	docker run -it -v $(FULL_DATA_DIR):$(DATA_DIR) $(IMG_NAME) $(CMD)

run:
	CMD="muscle -in $(DATA_DIR)/Ecoli_K12_MG1655.5UTR.mRNA.seq.cdhit -html -out $(DATA_DIR)/align.html" $(MAKE) run-docker