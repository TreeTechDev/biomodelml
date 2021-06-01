.PHONY: build workflow

DATA_DIR="/data"
TRIMM_DIR=$(DATA_DIR)/trimmomatic/adapters
FULL_DATA_DIR=`pwd`$(DATA_DIR)

IMG_NAME="bioinfo2.0"

build:
	@docker build . -t $(IMG_NAME)

run-docker:
	@docker run -it -v $(FULL_DATA_DIR):$(DATA_DIR) $(IMG_NAME) $(CMD)

%.fastq.gz:
	@-rm $(FULL_DATA_DIR)/SRR8173221_1.fastq.gz $(FULL_DATA_DIR)/SRR8173221_2.fastq.gz
	@CMD="fastq-dump --split-3 --gzip --outdir $(DATA_DIR) SRR8173221 -v" $(MAKE) run-docker

%/trimmomatic:
	@-rm $(FULL_DATA_DIR)/trimmomatic
	@CMD="git clone https://github.com/timflutre/trimmomatic $(DATA_DIR)/trimmomatic" $(MAKE) run-docker

%/ALL_PE.fa: data/trimmomatic
	@-rm $(FULL_DATA_DIR)/ALL_PE.fa
	@CMD="cat $(TRIMM_DIR)/NexteraPE-PE.fa $(TRIMM_DIR)/TruSeq2-PE.fa $(TRIMM_DIR)/TruSeq3-PE-2.fa $(TRIMM_DIR)/TruSeq3-PE.fa > ALL_PE.fa" $(MAKE) run-docker
	@mv ALL_PE.fa $(FULL_DATA_DIR)/ALL_PE.fa

%/trim.stats: data/SRR8173221_1.fastq.gz data/SRR8173221_2.fastq.gz data/ALL_PE.fa
	@-rm $(FULL_DATA_DIR)/trim.stats 
	@CMD="trimmomatic PE -threads 4 -phred33 -trimlog $(DATA_DIR)/trim.log -summary $(DATA_DIR)/trim.stats -quiet $(DATA_DIR)/SRR8173221_1.fastq.gz $(DATA_DIR)/SRR8173221_2.fastq.gz $(DATA_DIR)/SRR8173221_1.paired.fastq.gz $(DATA_DIR)/SRR8173221_1.unpaired.fastq.gz $(DATA_DIR)/SRR8173221_2.paired.fastq.gz $(DATA_DIR)/SRR8173221_2.unpaired.fastq.gz ILLUMINACLIP:$(DATA_DIR)/ALL_PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20" $(MAKE) run-docker

workflow: data/trim.stats
