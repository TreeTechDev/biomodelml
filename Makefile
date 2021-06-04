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

%.html: data/trim.stats
	@-rm $(FULL_DATA_DIR)/*.html
	@CMD="fastqc -q -t 4 $(DATA_DIR)/SRR8173221_1.paired.fastq.gz $(DATA_DIR)/SRR8173221_1.unpaired.fastq.gz $(DATA_DIR)/SRR8173221_2.paired.fastq.gz $(DATA_DIR)/SRR8173221_2.unpaired.fastq.gz" $(MAKE) run-docker

%.fna:
	@bash gdrive_download.sh https://drive.google.com/file/d/1_o7WjLJ7DLRPLLJUBjSofj46SxbDAbNP/view

%.bt2: data/Escherichia_coli_str_K-12_substr_MG1655_complete_genome_NCBI_Reference_Sequence_NC_000913.3.fna
	@-rm $(FULL_DATA_DIR)/*.bt2
	@CMD="bowtie2-build --quiet $(DATA_DIR)/Escherichia_coli_str_K-12_substr_MG1655_complete_genome_NCBI_Reference_Sequence_NC_000913.3.fna $(DATA_DIR)/Ecoli_K12_MG1655" $(MAKE) run-docker

%.sam: data/Ecoli_K12_MG1655.1.bt2 data/Ecoli_K12_MG1655.2.bt2 data/Ecoli_K12_MG1655.3.bt2 data/Ecoli_K12_MG1655.4.bt2 data/Ecoli_K12_MG1655.rev.1.bt2 data/Ecoli_K12_MG1655.rev.2.bt2
	@CMD="bowtie2 --quiet --rg-id --very-sensitive -x $(DATA_DIR)/Ecoli_K12_MG1655 -1 $(DATA_DIR)/SRR8173221_1.paired.fastq.gz -2 $(DATA_DIR)/SRR8173221_2.paired.fastq.gz -S $(DATA_DIR)/Ecoli_K12_MG1655.sam" $(MAKE) run-docker

%/picard.log: data/Ecoli_K12_MG1655.sam
	@CMD="picard ValidateSamFile --INPUT $(DATA_DIR)/Ecoli_K12_MG1655.sam --OUTPUT $(DATA_DIR)/picard.log --REFERENCE_SEQUENCE $(DATA_DIR)/Escherichia_coli_str_K-12_substr_MG1655_complete_genome_NCBI_Reference_Sequence_NC_000913.3.fna --MODE SUMMARY" $(MAKE) run-docker

%/Ecoli_K12_MG1655.bam: data/picard.log
	@CMD="samtools view --threads 4 --output-fmt bam $(DATA_DIR)/Ecoli_K12_MG1655.sam -o $(DATA_DIR)/Ecoli_K12_MG1655.bam" $(MAKE) run-docker

%/Ecoli_K12_MG1655.sort.bam: data/Ecoli_K12_MG1655.bam
	@CMD="samtools sort --threads 4 --output-fmt bam $(DATA_DIR)/Ecoli_K12_MG1655.bam -o $(DATA_DIR)/Ecoli_K12_MG1655.sort.bam" $(MAKE) run-docker

%/Ecoli_K12_MG1655.sort.bam.index: data/Ecoli_K12_MG1655.sort.bam
	@CMD="samtools index -b -@ 4 $(DATA_DIR)/Ecoli_K12_MG1655.sort.bam $(DATA_DIR)/Ecoli_K12_MG1655.sort.bam.index" $(MAKE) run-docker

%/Escherichia_coli_str_K-12_substr_MG1655_complete_genome_NCBI_Reference_Sequence_NC_000913.3.gff3:
	@bash gdrive_download.sh https://drive.google.com/file/d/1Wqf9zHLd6bOc_eN3nbs8y1L3P1JDvZFe/view

%/Ecoli_K12_MG1655.CDS.gff3: data/Escherichia_coli_str_K-12_substr_MG1655_complete_genome_NCBI_Reference_Sequence_NC_000913.3.gff3
	@-CMD='grep -P "NC_000913\.3\tRefSeq\tCDS" $(DATA_DIR)/Escherichia_coli_str_K-12_substr_MG1655_complete_genome_NCBI_Reference_Sequence_NC_000913.3.gff3 > Ecoli_K12_MG1655.CDS.gff3' $(MAKE) run-docker
	@mv Ecoli_K12_MG1655.CDS.gff3 $(FULL_DATA_DIR)/Ecoli_K12_MG1655.CDS.gff3

%/Escherichia_coli_str_K-12_substr_MG1655_complete_genome_NCBI_Reference_Sequence_NC_000913.3.fna.fai: data/Escherichia_coli_str_K-12_substr_MG1655_complete_genome_NCBI_Reference_Sequence_NC_000913.3.fna
	@CMD="samtools faidx $(DATA_DIR)/Escherichia_coli_str_K-12_substr_MG1655_complete_genome_NCBI_Reference_Sequence_NC_000913.3.fna" $(MAKE) run-docker

%/ecoli.genome: data/Escherichia_coli_str_K-12_substr_MG1655_complete_genome_NCBI_Reference_Sequence_NC_000913.3.fna.fai
	@bash make_genome.sh $(FULL_DATA_DIR)/Escherichia_coli_str_K-12_substr_MG1655_complete_genome_NCBI_Reference_Sequence_NC_000913.3.fna.fai $(FULL_DATA_DIR)/ecoli.genome

%/Ecoli_K12_MG1655.CDS.UTR.gff: data/ecoli.genome data/Ecoli_K12_MG1655.CDS.gff3
	@CMD="bedtools flank -s -l 10 -r 50 -i $(DATA_DIR)/Ecoli_K12_MG1655.CDS.gff3 -g $(DATA_DIR)/ecoli.genome > Ecoli_K12_MG1655.CDS.UTR.gff" $(MAKE) run-docker
	@mv Ecoli_K12_MG1655.CDS.UTR.gff $(FULL_DATA_DIR)/Ecoli_K12_MG1655.CDS.UTR.gff

%/Ecoli_K12_MG1655.5UTR3.gff: data/Ecoli_K12_MG1655.CDS.UTR.gff data/Ecoli_K12_MG1655.CDS.gff3
	@CMD="bedtools subtract -s -F 1.0 -a $(DATA_DIR)/Ecoli_K12_MG1655.CDS.UTR.gff -b $(DATA_DIR)/Ecoli_K12_MG1655.CDS.gff3 > Ecoli_K12_MG1655.5UTR3.gff" $(MAKE) run-docker
	@mv Ecoli_K12_MG1655.5UTR3.gff $(FULL_DATA_DIR)/Ecoli_K12_MG1655.5UTR3.gff

%/Ecoli_K12_MG1655.5UTR3.cov: data/Ecoli_K12_MG1655.sort.bam data/Ecoli_K12_MG1655.5UTR3.gff
	@sort -k 4,4n $(FULL_DATA_DIR)/Ecoli_K12_MG1655.5UTR3.gff > $(FULL_DATA_DIR)/Ecoli_K12_MG1655.sort.5UTR3.bed
	@CMD="bedtools coverage -sorted -s -f 0.30 -a $(DATA_DIR)/Ecoli_K12_MG1655.sort.5UTR3.bed -b $(DATA_DIR)/Ecoli_K12_MG1655.sort.bam > Ecoli_K12_MG1655.5UTR3.cov" $(MAKE) run-docker
	@mv Ecoli_K12_MG1655.5UTR3.cov $(FULL_DATA_DIR)/Ecoli_K12_MG1655.5UTR3.cov

workflow: data/Ecoli_K12_MG1655.sort.bam.index data/SRR8173221_1.paired_fastqc.html data/SRR8173221_1.unpaired_fastqc.html data/SRR8173221_2.paired_fastqc.html data/SRR8173221_2.unpaired_fastqc.html data/Ecoli_K12_MG1655.5UTR3.cov
