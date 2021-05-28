.PHONY: build workflow

DATA_DIR="/data"

FULL_DATA_DIR=`pwd`$(DATA_DIR)

IMG_NAME="bioinfo2.0"

build:
	@docker build . -t $(IMG_NAME)

run-docker:
	@docker run -it -v $(FULL_DATA_DIR):$(DATA_DIR) $(IMG_NAME) $(CMD)

fastq-dump:
	@CMD="fastq-dump --split-3 --gzip --outdir $(DATA_DIR) SRR390728 -v" $(MAKE) run-docker

workflow: | fastq-dump
