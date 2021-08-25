# Variables
DATE := $(shell date "+%y%m%d")
RANDOM := $(shell bash -c 'echo $$RANDOM')
IMAGE_NAME = lungcancer:latest
SHELL := /bin/bash
CPUS = 1
MEMS = 100
PARALLEL = 2
PWD := $(shell pwd)
TOOLS =
DOCKER = docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) $(IMAGE_NAME)
ARRIBA_DOCKER = docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) --workdir $(addprefix /,$(@D)) $(IMAGE_NAME)
PYTHON = python3 -B
R = Rscript --vanilla
CHROMOSOMES = chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
PAIRED_NORMAL = $(shell echo $* | grep --perl-regexp --only-match "^(cn)?\d+")N
PAIRED_PRIMARY = $(shell echo $* | grep --perl-regexp --only-match "^(cn)?\d+")P$(shell echo $* | grep --perl-regexp --only-match "\d$$" > /dev/null && echo "1")

# Options
VOLUME_OPTS = --volume $(abspath Output):/Output --volume $(abspath Data):/Data:ro --volume $(abspath Python):/Python:ro --volume $(abspath R):/R:ro --volume $(abspath Bash):/Bash:ro --tmpfs $(TMPFS)
RUN_OPTS = --tty --cpus="$(CPUS)" --memory="$(MEMS)G"
JAVA_OPTS = -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Xmx$(MEMS)g
GATK = gatk --java-options "$(JAVA_OPTS)"
PICARD = java $(JAVA_OPTS) -jar /Tools/picard.jar
VCF2MAF = perl /Tools/vcf2maf-1.6.19/vcf2maf.pl
MUTENRICHER = $(PYTHON) /Tools/MutEnricher-1.3.1/mutEnricher.py
TMPFS = /tmpfs
CIBERSORTx_TOKEN = --username "jwlee230@unist.ac.kr" --token "23aabc1d937bb0ae7fe0d4bf50c04b0c"
PYCLONE = /root/anaconda2/envs/pyclone/bin/PyClone

# Reference Files
REFERENCE_DIR = Data/hg38
REFERENCE_FASTA = $(REFERENCE_DIR)/Homo_sapiens_assembly38.fasta
REFERENCE_KNOWN_VCF = $(sort $(wildcard $(REFERENCE_DIR)/*.vcf))
REFERENCE_DBSNP_VCF = $(REFERENCE_DIR)/dbsnp_146.hg38.vcf
REFERENCE_1000G_VCF = $(REFERENCE_DIR)/1000G_phase1.snps.high_confidence.hg38.vcf
REFERENCE_HAPMAP_VCF = $(REFERENCE_DIR)/hapmap_3.3.hg38.vcf.gz
REFERENCE_OMNI_VCF = $(REFERENCE_DIR)/1000G_omni2.5.hg38.vcf.gz
REFERENCE_SIZE_TSV = $(REFERENCE_DIR)/hg38.chrom.sizes
REFERENCE_CENTROMERE = Data/UCSC_hg38_centromere.txt
LIBRARY_BED = Data/S04380110_Regions.bed
LIBRARY_BED_GZ = $(LIBRARY_BED).gz
GENCODE_GTF = Data/GENCODE/gencode.v38.annotation.gtf
TREMBL_TSV = Data/Homo_sapiens.GRCh38.104.uniprot.tsv.gz
GENE_CENSUS_CSV = Data/cancer_gene_census.csv
LUNG_SINGLECELL = Output/CIBERSORTx/GSE131907_Lung_Cancer_raw_UMI_matrix.txt
KNOWN_GENE_GTF = Data/hg38.gtf.gz
ARRIBA_STAR_DIR = Tools/arriba_v2.1.0/STAR_index_hg38_RefSeq_hg38
ARRIBA_DATABASE_DIR = Tools/arriba_v2.1.0/database
ARRIBA_BLACKLIST_TSV = $(ARRIBA_DATABASE_DIR)/blacklist_hg38_GRCh38_v2.1.0.tsv.gz
ARRIBA_KNOWN_FUSIONS_TSV = $(ARRIBA_DATABASE_DIR)/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz
ARRIBA_PROTEIN_DOMAIN_GFF3 = $(ARRIBA_DATABASE_DIR)/protein_domains_hg38_GRCh38_v2.1.0.gff3
FUSIONCATCHER_DATA_DIR = opt/fusioncatcher/data/current

# Sample IDs
PATIENT_ID = $(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+"))
WES_ID = $(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+[PCADMN]\d?"))
WES_NORMAL_ID = $(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+N"))
WES_TUMOR_ID = $(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+[PCADM]\d?"))
WTS_ID = $(sort $(shell ls -1 Data/WTS/*.fastq.gz | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+[PCADMN]"))
WTS_NORMAL_ID = $(sort $(shell ls -1 Data/WTS/*.fastq.gz | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+N"))
WTS_TUMOR_ID = $(sort $(shell ls -1 Data/WTS/*.fastq.gz | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+[PCADM]"))
BAM_ID = $(sort $(shell ls -1 Data/BAM/*.bam | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+[PCADMN](_?\d)?"))
BAM_NORMAL_ID = $(sort $(shell ls -1 Data/BAM/*.bam | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+N(_?\d)?"))
BAM_TUMOR_ID = $(sort $(shell ls -1 Data/BAM/*.bam | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+[PCADM](_?\d)?"))

# SQC / ADC
SQC_PATIENTS = 1 2 4 8 12 14 15 33 34 38 48 49 53 57 58 63 64 cn1 cn2 cn3 cn4 cn5 cn6 cn7 cn8 cn10 cn11 cn13 cn15 cn16 cn17 cn18 cn21 cn23 cn26 cn27 cn30 cn31 cn32 cn33 cn34 cn35 cn36 cn37 cn38 cn39 cn40 cn41 cn42 cn43 cn44 cn45 cn47 cn48 cn49 cn50 cn51 cn52 cn53 cn54 cn55 cn56 cn57 cn60 cn61 cn62 cn63 cn64 cn65 cn66 cn67 cn68 cn69 cn70 cn71 cn72 cn73 cn74 cn75 cn76 cn77 cn78 cn79
SQC_WES_SAMPLES = 8P1 8C1 12P1 12C1 14P1 14C1 15P1 15C1 33D1 34P1 34C1 38P1 38C1 48P1 48A1 48A2 48A3 49P1 49C1 53P1 53A1 57P1 57C1 57D1 57D2 58P1 58A1 58A2 58A3 58A4 63P1 63C1 63C2 64P1 64D1 cn1P cn1C cn1D cn2P cn2C cn3P cn3C cn4P cn4C cn5P cn5C cn7P cn7C cn10P cn10C cn11P cn11C cn13P cn13C cn15P cn15C cn16P cn16C cn17P cn17C cn18P cn18C cn21P cn21C cn23P cn23C cn26P cn26C cn27P cn27C cn30P cn30C cn31P cn31C cn32P cn32C cn33P cn33C cn34P cn34C cn35P cn35C cn36P cn36C cn37P cn37C cn38P cn38C cn39P cn39C cn40P cn40C cn41P cn41C cn42P cn42C cn43P cn43C cn44P cn44C cn45P cn45C cn47P cn47C cn48P cn48C cn49P cn49C cn50P cn50C cn51P cn51C cn52P cn52C cn53P cn53C cn54P cn54C cn55P cn55C cn56P cn56C cn57P cn57C cn60P cn60C cn61P cn61C cn62P cn62C cn63P cn63C cn64P cn64C cn65P cn65C cn66P cn66C cn67P cn67C cn68P cn68C cn69P cn69C cn70P cn70C cn71P cn71C cn72P cn72C cn73P cn73C cn74P cn74C cn75P cn75C cn76P cn76C cn78P cn78C cn79P cn79C
SQC_WTS_SAMPLES = cn1P cn1C cn1D cn1N cn2P cn2C cn3P cn3C cn3N cn4P cn4N cn5P cn5C cn5N cn6P cn6C cn6N cn7P cn7C cn8P cn8C cn10P cn10C cn16P cn16C cn17P cn17C cn23P cn23C cn30P cn30C cn32P cn32C cn38P cn38C cn39P cn39C cn40P cn40C cn44P cn44C cn49P cn49C cn52P cn52C cn52N cn53P cn53C cn54P cn54C cn55N cn57P cn57C cn57N cn60N cn62N cn63P cn63C cn63N cn64P cn64C cn64N cn65P cn65C cn66P cn66C cn66N cn67P cn67C cn67N cn68N cn69P cn69C cn69N cn70P cn70C cn74P cn74C cn75P cn75C cn75N cn76P cn76C cn77P cn77D cn78P cn78C cn79P cn79C
ADC_PATIENTS = 39 51 52 54 55 61 62 65 cn80 cn90 cn93 cn95 cn96 cn97 cn99 cn101 cn102 cn103 cn105 cn106 cn107 cn108 cn109 cn110 cn111 cn112 cn113 cn115 cn116
ADC_WES_SAMPLES = 39P1 39A1 39A2 51P1 51A1 52P1 52M1 54P1 54A1 55P1 55A1 55A2 61P1 61A1 61A2 62P1 62A1 62A2 62A3 62A4 65P1 65A1 cn93P cn93C cn95P cn95A cn97P cn97C cn101P cn101C cn105P cn105C cn106P cn106C cn107P cn107C cn107A1 cn108P cn108C cn110P cn110C cn112P cn112C
ADC_WTS_SAMPLES = cn90P cn90C cn90N cn95N cn96N cn99N cn101N cn102P cn102C cn102N cn103N cn107P cn107C cn109N cn110P cn110C cn111N cn112N cn113P cn113C cn113N cn115P cn115A cn115N cn116N

# General
all:
.PHONY: all

log Output Docker/Tools:
	mkdir -p $@

# Tools
Docker/Tools/Python.tar.xz: | Docker/Tools
	wget "https://www.python.org/ftp/python/3.9.1/Python-3.9.1.tar.xz" -O $@
TOOLS += Docker/Tools/Python.tar.xz

Docker/Tools/samtools.tar.bz2: | Docker/Tools
	wget "https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2" -O $@
TOOLS += Docker/Tools/samtools.tar.bz2

Docker/Tools/gatk.zip: | Docker/Tools
	wget "https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip" -O $@
TOOLS += Docker/Tools/gatk.zip

Docker/Tools/bwa.tar.bz2: | Docker/Tools
	wget "https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2" -O $@
TOOLS += Docker/Tools/bwa.tar.bz2

Docker/Tools/picard.jar: | Docker/Tools
	wget "https://github.com/broadinstitute/picard/releases/download/2.24.1/picard.jar" -O $@
TOOLS += Docker/Tools/picard.jar

Docker/Tools/FastQC.zip: | Docker/Tools
	wget "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip" -O $@
TOOLS += Docker/Tools/FastQC.zip

Docker/Tools/STAR.tar.gz: | Docker/Tools
	wget "https://github.com/alexdobin/STAR/archive/2.7.8a.tar.gz" -O $@
TOOLS += Docker/Tools/STAR.tar.gz

Docker/Tools/Bowtie2.zip: | Docker/Tools
	wget "https://github.com/BenLangmead/bowtie2/releases/download/v2.4.2/bowtie2-2.4.2-linux-x86_64.zip" -O $@
TOOLS += Docker/Tools/Bowtie2.zip

Docker/Tools/vcf2maf.tar.gz: | Docker/Tools
	wget "https://github.com/mskcc/vcf2maf/archive/v1.6.19.tar.gz" -O $@
TOOLS += Docker/Tools/vcf2maf.tar.gz

Docker/Tools/VEP.zip: | Docker/Tools
	wget "https://github.com/Ensembl/ensembl-vep/archive/refs/tags/release/103.zip" -O $@
TOOLS += Docker/Tools/VEP.zip

Docker/Tools/BigFile.tar.gz: | Docker/Tools
	wget "https://github.com/ucscGenomeBrowser/kent/archive/v335_base.tar.gz" -O $@
TOOLS += Docker/Tools/BigFile.tar.gz

Docker/Tools/MutEnricher.tar.gz: | Docker/Tools
	wget "https://github.com/asoltis/MutEnricher/archive/refs/tags/v1.3.1.tar.gz" -O $@
TOOLS += Docker/Tools/MutEnricher.tar.gz

Docker/Tools/bcftools.tar.bz2: | Docker/Tools
	wget "https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2" -O $@
TOOLS += Docker/Tools/bcftools.tar.bz2

Docker/Tools/RSEM.tar.gz: | Docker/Tools
	wget "https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz" -O $@
TOOLS += Docker/Tools/RSEM.tar.gz

Docker/Tools/Lemon.tar.gz: | Docker/Tools
	wget "http://lemon.cs.elte.hu/pub/sources/lemon-1.3.tar.gz" -O $@
TOOLS += Docker/Tools/Lemon.tar.gz

Docker/Tools/Miniconda.sh: | Docker/Tools
	wget "https://repo.anaconda.com/archive/Anaconda2-2019.10-Linux-x86_64.sh" -O $@
TOOLS += Docker/Tools/Miniconda.sh

Docker/Tools/arriba.tar.gz: | Docker/Tools
	wget "https://github.com/suhrig/arriba/releases/download/v2.1.0/arriba_v2.1.0.tar.gz" -O $@
TOOLS += Docker/Tools/arriba.tar.gz

Docker/Tools/FusionCatcher.py: | Docker/Tools
	wget "http://sf.net/projects/fusioncatcher/files/bootstrap.py" -O $@
TOOLS += Docker/Tools/FusionCatcher.py

Docker/Tools/Star-Fusion.tar.gz: | Docker/Tools
	wget "https://github.com/STAR-Fusion/STAR-Fusion/releases/download/STAR-Fusion-v1.10.1/STAR-Fusion.v1.10.1.tar.gz" -O $@
TOOLS += Docker/Tools/Star-Fusion.tar.gz

# Docker
build.log: Docker/Dockerfile Docker/requirements.txt $(TOOLS) | log Output
	rm -f $@
	docker images | grep $(IMAGE_NAME) && docker rmi $(IMAGE_NAME) || true
	docker build --rm --tag $(IMAGE_NAME) $(dir $(word 1,$^)) | tee $@

build: build.log
.PHONY: build

interactive: build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) --interactive $(IMAGE_NAME) /bin/bash || true
.PHONY: interactive

delete: build.log
	docker rmi $(IMAGE_NAME)
	rm -rfv build.log
.PHONY: delete

stop:
	docker rm $(CONTAINER_NAME)
.PHONY: stop

# SGE
run: | log Output build.log
	echo "make -j $(PARALLEL) -C $(PWD) latest" > tmp.sh
	qsub -cwd -l h_vmem=$(shell echo "$(MEMS) * $(PARALLEL)" | bc)G -m abe -M "230@fumire.moe" -N Lung_$(DATE) -pe smp $(shell echo "$(CPUS) * $(PARALLEL)" | bc) -o $(abspath log) -e $(abspath log) tmp.sh
.PHONY: run

# Actual
latest: step03 step04 step05 step06 step07 step08 step09 step10 step11 step12 step13 step14 step15 step16 step17 step18 step19 step21 step22 step23 step24 step25
.PHONY: latest

# Step 01 (FastQC for WES)
Output/FastQC/WES:
	mkdir -p $@

Output/FastQC/WES/%_fastqc.zip: Data/RawData/%.fastq.gz | Output/FastQC/WES build.log
	$(DOCKER) fastqc --outdir $(addprefix /,$(@D)) --format fastq --noextract $(sort $(addprefix /,$^)) --dir $(TMPFS) --threads $(CPUS) 1> $@.stdout 2> $@.stderr

Output/FastQC/FastQC_WES.pdf: Python/aggregate_fastqc.py $(addprefix Output/FastQC/WES/,$(addsuffix _fastqc.zip,$(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | cut --delimiter . --fields 1)))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --threshold 3 1> $@.stdout 2> $@.stderr

step01: $(addprefix Output/FastQC/WES/,$(addsuffix _fastqc.zip,$(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | cut --delimiter . --fields 1)))) Output/FastQC/FastQC_WES.pdf
.PHONY: step01

# Step 02 (FastQC for WTS)
Output/FastQC/WTS:
	mkdir -p $@

Output/FastQC/WTS/%_fastqc.zip: Data/WTS/%.fastq.gz | Output/FastQC/WTS build.log
	$(DOCKER) fastqc --outdir $(addprefix /,$(@D)) --format fastq --noextract $(sort $(addprefix /,$^)) --dir $(TMPFS) --threads $(CPUS) 1> $@.stdout 2> $@.stderr

Output/FastQC/FastQC_WTS.pdf: Python/aggregate_fastqc.py $(addprefix Output/FastQC/WTS/,$(addsuffix _fastqc.zip,$(sort $(shell ls -1 Data/WTS/*.fastq.gz | cut --delimiter / --fields 3 | cut --delimiter . --fields 1)))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --threshold 5 1> $@.stdout 2> $@.stderr

step02: $(addprefix Output/FastQC/WTS/,$(addsuffix _fastqc.zip,$(sort $(shell ls -1 Data/WTS/*.fastq.gz | cut --delimiter / --fields 3 | cut --delimiter . --fields 1)))) Output/FastQC/FastQC_WTS.pdf
.PHONY: step02

# Step 03 (WES: fq.gz >> sam)
Output/WES:
	mkdir -p $@

Output/WES/%.BWA.sam: $(REFERENCE_FASTA) Data/RawData/%_*.fastq.gz | Output/WES build.log
	$(DOCKER) bwa mem -M -t $(CPUS) -R "@RG\tID:$(@F)\tPL:ILLUMINA\tLB:$(@F)\tSM:$(@F)\tCN:UNIST" -v 3 $(addprefix /,$(filter %.fasta,$^)) $(sort $(addprefix /,$(filter %.fastq.gz,$^))) 1> $@ 2> $@.stderr
	$(DOCKER) sed --in-place "/^\[.*/d" $(addprefix /,$@)

Output/WES/%.Bowtie2.sam: $(REFERENCE_FASTA) Data/RawData/%_*.fastq.gz | Output/WES build.log
	$(DOCKER) bowtie2 --threads $(CPUS) --rg-id "$(@F)" --rg "ID:$(@F)" --rg "PL:ILLUMINA" --rg "LB:$(@F)" --rg "SM:$(@F)" --rg "CN:UNIST" --time --qc-filter -x $(basename $(notdir $(addprefix /,$(filter %.fasta,$^)))) -1 $(addprefix /,$(word 1,$(sort $(filter %.fastq.gz,$^)))) -2 $(addprefix /,$(word 2,$(sort $(filter %.fastq.gz,$^)))) -S $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step03: $(addprefix Output/WES/,$(addsuffix .BWA.sam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .Bowtie2.sam,$(WES_ID)))
.PHONY: step03

# Step 04 (WES: sam >> bam)
Output/WES/%.bam: $(REFERENCE_FASTA) Output/WES/%.sam | build.log
	$(DOCKER) samtools view -b -h --threads $(CPUS) --reference $(addprefix /,$(filter %.fasta,$^)) $(addprefix /,$(filter %.sam,$^)) -o $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step04: $(addprefix Output/WES/,$(addsuffix .BWA.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .Bowtie2.bam,$(WES_ID)))
.PHONY: step04

# Step 05 (WES: bam >> Sort.bam)
Output/WES/%.Sort.bam: $(REFERENCE_FASTA) Output/WES/%.bam | build.log
	$(DOCKER) samtools sort -l 9 -m $(shell echo "$(MEMS) / $(CPUS)" | bc)G --threads $(CPUS) --reference $(addprefix /,$(filter %.fasta,$^)) -o $(addprefix /,$@) $(addprefix /,$(filter %.bam,$^)) 1> $@.stdout 2> $@.stderr

Output/WES/%.Sort.bai: Output/WES/%.Sort.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step05: $(addprefix Output/WES/,$(addsuffix .BWA.Sort.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .BWA.Sort.bai,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.bai,$(WES_ID)))
.PHONY: step05

# Step 06 (WES: Sort.bam >> Sort.MarkDuplicate.bam)
Output/WES/%.Sort.dm-tmp.bam: $(REFERENCE_FASTA) Output/WES/%.Sort.bam Output/WES/%.Sort.bai | build.log
	$(DOCKER) $(PICARD) MarkDuplicates --INPUT $(addprefix /,$(filter %.bam,$^)) --REFERENCE_SEQUENCE $(addprefix /,$(filter %.fasta,$^)) --OUTPUT $(addprefix /,$@) --METRICS_FILE $(addprefix /,$@.metrics) --ASSUME_SORTED true --VALIDATION_STRINGENCY SILENT --COMPRESSION_LEVEL 9 1> $@.stdout 2> $@.stderr

Output/WES/%.Sort.MarkDuplicate.bam: Bash/remove_CR.bash Output/WES/%.Sort.dm-tmp.bam | build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/WES/%.Sort.MarkDuplicate.bai: Output/WES/%.Sort.MarkDuplicate.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step06: $(addprefix Output/WES/,$(addsuffix .BWA.Sort.dm-tmp.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.bai,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.dm-tmp.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.bai,$(WES_ID)))
.PHONY: step06

# Step 07 (WES: Sort.MarkDuplicate.bam >> Sort.MarkDuplicate.BQSR.bam)
Output/WES/%.Sort.MarkDuplicate.BQSR.table: $(REFERENCE_FASTA) $(REFERENCE_KNOWN_VCF) Output/WES/%.Sort.MarkDuplicate.bam Output/WES/%.Sort.MarkDuplicate.bai | build.log
	$(DOCKER) $(GATK) BaseRecalibrator --input $(addprefix /,$(filter %.bam,$^)) $(foreach site,$(addprefix /,$(filter %.vcf,$^)),--known-sites $(site)) --output $(addprefix /,$@) --reference $(addprefix /,$(filter %.fasta,$^)) 1> $@.stdout 2> $@.stderr

Output/WES/%.Sort.MarkDuplicate.BQSR.bam: Output/WES/%.Sort.MarkDuplicate.bam Output/WES/%.Sort.MarkDuplicate.BQSR.table $(REFERENCE_FASTA) | build.log
	$(DOCKER) $(GATK) ApplyBQSR --bqsr-recal-file $(addprefix /,$(filter %.table,$^)) --input $(addprefix /,$(filter %.bam,$^)) --output $(addprefix /,$@) --reference $(addprefix /,$(filter %.fasta,$^)) 1> $@.stdout 2> $@.stderr

Output/WES/%.Sort.MarkDuplicate.BQSR.bai: Output/WES/%.Sort.MarkDuplicate.BQSR.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step07: $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.table,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.bai,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.table,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.bai,$(WES_ID)))
.PHONY: step07

clean-step07:
	rm -fv $(filter-out %.Sort.MarkDuplicate.BQSR.bam,$(shell find Output/WES -size +1G))
.PHONY: clean-step07

# Step 08 (Panel of Normal: Run Mutect2 in tumor-only mode for each normal sample)
Output/PanelOfNormal/BWA Output/PanelOfNormal/Bowtie2:
	mkdir -p $@

Output/PanelOfNormal/BWA/%.vcf: $(REFERENCE_FASTA) Output/WES/%.BWA.Sort.MarkDuplicate.BQSR.bam Output/WES/%.BWA.Sort.MarkDuplicate.BQSR.bai | Output/PanelOfNormal/BWA build.log
	$(DOCKER) $(GATK) Mutect2 --reference $(addprefix /,$(filter %.fasta,$^)) --input $(addprefix /,$(filter %.bam,$^)) --output $(addprefix /,$@) --max-mnp-distance 0 1> $@.stdout 2> $@.stderr

Output/PanelOfNormal/Bowtie2/%.vcf: $(REFERENCE_FASTA) Output/WES/%.Bowtie2.Sort.MarkDuplicate.BQSR.bam Output/WES/%.Bowtie2.Sort.MarkDuplicate.BQSR.bai | Output/PanelOfNormal/Bowtie2 build.log
	$(DOCKER) $(GATK) Mutect2 --reference $(addprefix /,$(filter %.fasta,$^)) --input $(addprefix /,$(filter %.bam,$^)) --output $(addprefix /,$@) --max-mnp-distance 0 1> $@.stdout 2> $@.stderr

step08: $(addprefix Output/PanelOfNormal/BWA/,$(addsuffix .vcf,$(WES_NORMAL_ID))) $(addprefix Output/PanelOfNormal/Bowtie2/,$(addsuffix .vcf,$(WES_NORMAL_ID)))
.PHONY: step08

# Step 09 (Panel of Normal: Aggregate VCF files)
Output/PanelOfNormal/BWA.genomicsdb/callset.json: $(REFERENCE_FASTA) $(addprefix Output/PanelOfNormal/BWA/,$(addsuffix .vcf,$(WES_NORMAL_ID))) | build.log
	$(DOCKER) $(GATK) GenomicsDBImport --reference $(addprefix /,$(filter %.fasta,$^)) --genomicsdb-workspace-path $(addprefix /,$(@D)) $(foreach sample,$(addprefix /,$(filter %.vcf,$^)),--variant $(sample)) --reader-threads $(CPUS) --max-num-intervals-to-import-in-parallel $(CPUS) --overwrite-existing-genomicsdb-workspace true --merge-input-intervals true $(foreach chr,$(CHROMOSOMES),--intervals $(chr)) 1> $(@D).stdout 2> $(@D).stderr
	$(DOCKER) chmod --recursive 777 $(addprefix /,$(@D))

Output/PanelOfNormal/Bowtie2.genomicsdb/callset.json: $(REFERENCE_FASTA) $(addprefix Output/PanelOfNormal/Bowtie2/,$(addsuffix .vcf,$(WES_NORMAL_ID))) | build.log
	$(DOCKER) $(GATK) GenomicsDBImport --reference $(addprefix /,$(filter %.fasta,$^)) --genomicsdb-workspace-path $(addprefix /,$(@D)) $(foreach sample,$(addprefix /,$(filter %.vcf,$^)),--variant $(sample)) --reader-threads $(CPUS) --max-num-intervals-to-import-in-parallel $(CPUS) --overwrite-existing-genomicsdb-workspace true --merge-input-intervals true $(foreach chr,$(CHROMOSOMES),--intervals $(chr)) 1> $(@D).stdout 2> $(@D).stderr
	$(DOCKER) chmod --recursive 777 $(addprefix /,$(@D))

Output/PanelOfNormal/PON.BWA.vcf.gz: $(REFERENCE_FASTA) Output/PanelOfNormal/BWA.genomicsdb/callset.json | build.log
	$(DOCKER) $(GATK) CreateSomaticPanelOfNormals --reference $(addprefix /,$(filter %.fasta,$^)) --variant "gendb://$(addprefix /,$(dir $(filter %.json,$^)))" --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/PanelOfNormal/PON.Bowtie2.vcf.gz: $(REFERENCE_FASTA) Output/PanelOfNormal/Bowtie2.genomicsdb/callset.json | build.log
	$(DOCKER) $(GATK) CreateSomaticPanelOfNormals --reference $(addprefix /,$(filter %.fasta,$^)) --variant "gendb://$(addprefix /,$(dir $(filter %.json,$^)))" --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step09: Output/PanelOfNormal/BWA.genomicsdb/callset.json Output/PanelOfNormal/PON.BWA.vcf.gz Output/PanelOfNormal/Bowtie2.genomicsdb/callset.json Output/PanelOfNormal/PON.Bowtie2.vcf.gz
.PHONY: step09

# Step 10 (WES - samtools depth)
Output/WESDepth:
	mkdir -p $@

Output/WESDepth/%.tsv: $(REFERENCE_FASTA) $(LIBRARY_BED) Output/WES/%.Sort.MarkDuplicate.BQSR.bam Output/WES/%.Sort.MarkDuplicate.BQSR.bai | Output/WESDepth build.log
	$(DOCKER) samtools depth -H -b $(addprefix /,$(filter %.bed,$^)) --reference $(addprefix /,$(filter %.fasta,$^)) -o $(addprefix /,$@) $(addprefix /,$(filter %.bam,$^)) 1> $@.stdout 2> $@.stderr

Output/WESDepth/BWA.pdf: Python/analysis_depth.py $(addprefix Output/WESDepth/,$(addsuffix .BWA.tsv,$(WES_ID))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/WESDepth/Bowtie2.pdf: Python/analysis_depth.py $(addprefix Output/WESDepth/,$(addsuffix .Bowtie2.tsv,$(WES_ID))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step10: $(addprefix Output/WESDepth/,$(addsuffix .BWA.tsv,$(WES_ID))) Output/WESDepth/BWA.pdf $(addprefix Output/WESDepth/,$(addsuffix .Bowtie2.tsv,$(WES_ID))) Output/WESDepth/Bowtie2.pdf
.PHONY: step10

# Step 11 (Sequenza: Make reference)
Output/Sequenza:
	mkdir -p $@

Output/Sequenza/hg38.wig.gz: $(REFERENCE_FASTA) | Output/Sequenza build.log
	$(DOCKER) sequenza-utils gc_wiggle --fasta $(addprefix /,$^) -o $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step11: Output/Sequenza/hg38.wig.gz
.PHONY: step11

# Step 12 (Sequenza: Normal-Tumor pairs)
Output/Sequenza/BWA Output/Sequenza/Bowtie2:
	mkdir -p $@

# I know that this recipe contains unnecessary dependencies; but I cannot see the solution.
Output/Sequenza/BWA/%.seqz.gz: $(REFERENCE_FASTA) Output/Sequenza/hg38.wig.gz Output/WES/%.BWA.Sort.MarkDuplicate.BQSR.bam $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.bam,$(WES_NORMAL_ID))) | Output/Sequenza/BWA build.log
	@echo "This step requires multiple CPUs!!"
	$(DOCKER) sequenza-utils bam2seqz --normal "/Output/WES/$(PAIRED_NORMAL).BWA.Sort.MarkDuplicate.BQSR.bam" --tumor "/Output/WES/$*.BWA.Sort.MarkDuplicate.BQSR.bam" --fasta $(addprefix /,$(filter %.fasta,$^)) -gc $(addprefix /,$(filter %.wig.gz,$^)) --output $(addprefix /,$@) --chromosome $(CHROMOSOMES) --parallel $(CPUS) 1> $@.stdout 2> $@.stderr
	cd $(@D) && gzip --decompress --stdout --keep $(foreach chr,$(CHROMOSOMES),$*_$(chr).seqz.gz) | gawk '{if (NR!=1 && $$1 != "chromosome") {print $$0}}' | bgzip --threads $(CPUS) > $(@F)
	$(DOCKER) rm $(addprefix /$(@D)/,$(foreach chr,$(CHROMOSOMES),$*_$(chr).seqz.gz))

Output/Sequenza/Bowtie2/%.seqz.gz: $(REFERENCE_FASTA) Output/Sequenza/hg38.wig.gz Output/WES/%.Bowtie2.Sort.MarkDuplicate.BQSR.bam $(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.bam,$(WES_NORMAL_ID))) | Output/Sequenza/Bowtie2 build.log
	@echo "This step requires multiple CPUs!!"
	$(DOCKER) sequenza-utils bam2seqz --normal "/Output/WES/$(PAIRED_NORMAL).BWA.Sort.MarkDuplicate.BQSR.bam" --tumor "/Output/WES/$*.BWA.Sort.MarkDuplicate.BQSR.bam" --fasta $(addprefix /,$(filter %.fasta,$^)) -gc $(addprefix /,$(filter %.wig.gz,$^)) --output $(addprefix /,$@) --chromosome $(CHROMOSOMES) --parallel $(CPUS) 1> $@.stdout 2> $@.stderr
	cd $(@D) && gzip --decompress --stdout --keep $(foreach chr,$(CHROMOSOMES),$*_$(chr).seqz.gz) | gawk '{if (NR!=1 && $$1 != "chromosome") {print $$0}}' | bgzip --threads $(CPUS) > $(@F)
	$(DOCKER) rm $(addprefix /$(@D)/,$(foreach chr,$(CHROMOSOMES),$*_$(chr).seqz.gz))

Output/Sequenza/BWA/%.seqz.gz.tbi: Output/Sequenza/BWA/%.seqz.gz | build.log
	$(DOCKER) tabix --force --sequence 1 --begin 2 --end 2 --skip-lines 1 $(addprefix /,$<) 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2/%.seqz.gz.tbi: Output/Sequenza/Bowtie2/%.seqz.gz | build.log
	$(DOCKER) tabix --force --sequence 1 --begin 2 --end 2 --skip-lines 1 $(addprefix /,$<) 1> $@.stdout 2> $@.stderr

step12: $(addprefix Output/Sequenza/BWA/,$(addsuffix .seqz.gz,$(WES_TUMOR_ID))) $(addprefix Output/Sequenza/BWA/,$(addsuffix .seqz.gz.tbi,$(WES_TUMOR_ID))) $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix .seqz.gz,$(WES_TUMOR_ID))) $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix .seqz.gz.tbi,$(WES_TUMOR_ID)))
.PHONY: step12

# Step 13 (Sequenza: Draw the plots)
Output/Sequenza/BWA/%.tar.gz: R/run_sequenza.R Output/Sequenza/BWA/%.seqz.gz Output/Sequenza/BWA/%.seqz.gz.tbi | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --sampleid $* --input $(addprefix /,$(filter %.seqz.gz,$^)) --outdir $(addprefix /,$(@D)) --parallel $(CPUS) 1> $@.stdout 2> $@.stderr
	cd $(abspath $(@D)) && tar -zcvf $(@F) $*_* && rm -fv $*_*

Output/Sequenza/Bowtie2/%.tar.gz: R/run_sequenza.R Output/Sequenza/Bowtie2/%.seqz.gz Output/Sequenza/Bowtie2/%.seqz.gz.tbi | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --sampleid $* --input $(addprefix /,$(filter %.seqz.gz,$^)) --outdir $(addprefix /,$(@D)) --parallel $(CPUS) 1> $@.stdout 2> $@.stderr
	cd $(abspath $(@D)) && tar -zcvf $(@F) $*_* && rm -fv $*_*

Output/Sequenza/BWA-sequenza-SQC.pdf: Python/aggregate_sequenza.py $(addprefix Output/Sequenza/BWA/,$(addsuffix .tar.gz,$(SQC_WES_SAMPLES))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-sequenza-ADC.pdf: Python/aggregate_sequenza.py $(addprefix Output/Sequenza/BWA/,$(addsuffix .tar.gz,$(ADC_WES_SAMPLES))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-SQC.pdf: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/BWA/,$(addsuffix .tar.gz,$(SQC_WES_SAMPLES))) Data/hg38/hg38.chrom.sizes | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-ADC.pdf: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/BWA/,$(addsuffix .tar.gz,$(ADC_WES_SAMPLES))) Data/hg38/hg38.chrom.sizes | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-simple-SQC.pdf: Python/aggregate_sequenza_simple.py $(addprefix Output/Sequenza/BWA/,$(addsuffix .tar.gz,$(SQC_WES_SAMPLES))) Data/hg38/hg38.chrom.sizes | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-simple-ADC.pdf: Python/aggregate_sequenza_simple.py $(addprefix Output/Sequenza/BWA/,$(addsuffix .tar.gz,$(ADC_WES_SAMPLES))) Data/hg38/hg38.chrom.sizes | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-sequenza-SQC.pdf: Python/aggregate_sequenza.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix .tar.gz,$(SQC_WES_SAMPLES))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-sequenza-ADC.pdf: Python/aggregate_sequenza.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix .tar.gz,$(ADC_WES_SAMPLES))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-genome-SQC.pdf: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix .tar.gz,$(SQC_WES_SAMPLES))) Data/hg38/hg38.chrom.sizes | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-genome-ADC.pdf: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix .tar.gz,$(ADC_WES_SAMPLES))) Data/hg38/hg38.chrom.sizes | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-simple-SQC.pdf: Python/aggregate_sequenza_simple.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix .tar.gz,$(SQC_WES_SAMPLES))) Data/hg38/hg38.chrom.sizes | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-simple-ADC.pdf: Python/aggregate_sequenza_simple.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix .tar.gz,$(ADC_WES_SAMPLES))) Data/hg38/hg38.chrom.sizes | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step13: $(addprefix Output/Sequenza/BWA/,$(addsuffix .tar.gz,$(WES_TUMOR_ID))) Output/Sequenza/BWA-sequenza-SQC.pdf Output/Sequenza/BWA-sequenza-ADC.pdf Output/Sequenza/BWA-genome-SQC.pdf Output/Sequenza/BWA-genome-ADC.pdf Output/Sequenza/BWA-simple-SQC.pdf Output/Sequenza/BWA-simple-ADC.pdf $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix .tar.gz,$(WES_TUMOR_ID))) Output/Sequenza/Bowtie2-sequenza-SQC.pdf Output/Sequenza/Bowtie2-sequenza-ADC.pdf Output/Sequenza/Bowtie2-genome-SQC.pdf Output/Sequenza/Bowtie2-genome-ADC.pdf Output/Sequenza/Bowtie2-simple-SQC.pdf Output/Sequenza/Bowtie2-simple-ADC.pdf
.PHONY: step13

# Step 14 (Mutect2)
Output/Mutect2/BWA Output/Mutect2/Bowtie2:
	mkdir -p $@

# I know that this recipe contains unnecessary dependencies; but I cannot see the solution.
Output/Mutect2/BWA/%.BWA.vcf: $(REFERENCE_FASTA) Output/PanelOfNormal/PON.BWA.vcf.gz Output/WES/%.BWA.Sort.MarkDuplicate.BQSR.bam $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.bam,$(WES_NORMAL_ID))) | Output/Mutect2/BWA build.log
	$(DOCKER) $(GATK) Mutect2 --reference $(addprefix /,$(filter %.fasta,$^)) --input "/Output/WES/$*.BWA.Sort.MarkDuplicate.BQSR.bam" --input "/Output/WES/$(PAIRED_NORMAL).BWA.Sort.MarkDuplicate.BQSR.bam" --output $(addprefix /,$@) --normal-sample "$(PAIRED_NORMAL).BWA.sam" --panel-of-normals $(addprefix /,$(filter %.vcf.gz,$^)) 1> $@.stdout 2> $@.stderr

# I know that this recipe contains unnecessary dependencies; but I cannot see the solution.
Output/Mutect2/Bowtie2/%.Bowtie2.vcf: $(REFERENCE_FASTA) Output/PanelOfNormal/PON.Bowtie2.vcf.gz Output/WES/%.Bowtie2.Sort.MarkDuplicate.BQSR.bam $(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.bam,$(WES_NORMAL_ID))) | Output/Mutect2/Bowtie2 build.log
	$(DOCKER) $(GATK) Mutect2 --reference $(addprefix /,$(filter %.fasta,$^)) --input "/Output/WES/$*.Bowtie2.Sort.MarkDuplicate.BQSR.bam" --input "/Output/WES/$(PAIRED_NORMAL).Bowtie2.Sort.MarkDuplicate.BQSR.bam" --output $(addprefix /,$@) --normal-sample "$(PAIRED_NORMAL).Bowtie2.sam" --panel-of-normals $(addprefix /,$(filter %.vcf.gz,$^)) 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA/%.BWA.filter.vcf: $(REFERENCE_FASTA) Output/Mutect2/BWA/%.BWA.vcf | build.log
	$(DOCKER) $(GATK) FilterMutectCalls --reference $(addprefix /,$(filter %.fasta,$^)) --variant $(addprefix /,$(filter %.vcf,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2/%.Bowtie2.filter.vcf: $(REFERENCE_FASTA) Output/Mutect2/Bowtie2/%.Bowtie2.vcf | build.log
	$(DOCKER) $(GATK) FilterMutectCalls --reference $(addprefix /,$(filter %.fasta,$^)) --variant $(addprefix /,$(filter %.vcf,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA/%.BWA.PASS.vcf: Bash/select_pass.bash Output/Mutect2/BWA/%.BWA.filter.vcf | build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2/%.Bowtie2.PASS.vcf: Bash/select_pass.bash Output/Mutect2/Bowtie2/%.Bowtie2.filter.vcf | build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA/%.BWA.PASS.maf: $(REFERENCE_FASTA) Output/Mutect2/BWA/%.BWA.PASS.vcf | build.log
	$(DOCKER) $(VCF2MAF) --vep-path "/Tools/ensembl-vep-release-103" --vep-data "/Data" --vep-forks $(CPUS) --ncbi-build "GRCh38" --input-vcf $(addprefix /,$(filter %.vcf,$^)) --output $(addprefix /,$@) --tumor-id "$*.BWA.sam" --normal-id "$(PAIRED_NORMAL).BWA.sam" --ref-fasta $(addprefix /,$(filter %.fasta,$^)) --tmp-dir $(TMPFS) 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2/%.Bowtie2.PASS.maf: $(REFERENCE_FASTA) Output/Mutect2/Bowtie2/%.Bowtie2.PASS.vcf | build.log
	$(DOCKER) $(VCF2MAF) --vep-path "/Tools/ensembl-vep-release-103" --vep-data "/Data" --vep-forks $(CPUS) --ncbi-build "GRCh38" --input-vcf $(addprefix /,$(filter %.vcf,$^)) --output $(addprefix /,$@) --tumor-id "$*.Bowtie2.sam" --normal-id "$(PAIRED_NORMAL).Bowtie2.sam" --ref-fasta $(addprefix /,$(filter %.fasta,$^)) --tmp-dir $(TMPFS) 1> $@.stdout 2> $@.stderr

step14: $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.vcf,$(WES_TUMOR_ID))) $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.filter.vcf,$(WES_TUMOR_ID))) $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.vcf,$(WES_TUMOR_ID))) $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.vcf,$(WES_TUMOR_ID))) $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.filter.vcf,$(WES_TUMOR_ID))) $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.vcf,$(WES_TUMOR_ID))) $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID)))
.PHONY: step14

# Step 15 (HaplotypeCaller)
Output/HaplotypeCaller/BWA Output/HaplotypeCaller/Bowtie2:
	mkdir -p $@

Output/HaplotypeCaller/BWA/%.g.vcf.gz: $(REFERENCE_FASTA) Output/WES/%.BWA.Sort.MarkDuplicate.BQSR.bam $(REFERENCE_DBSNP_VCF) | Output/HaplotypeCaller/BWA build.log
	$(DOCKER) $(GATK) HaplotypeCaller --reference $(addprefix /,$(filter %.fasta,$^)) --input $(addprefix /,$(filter %.bam,$^)) --output $(addprefix /,$@) --emit-ref-confidence GVCF --dbsnp $(addprefix /,$(filter %.vcf,$^)) --native-pair-hmm-threads $(CPUS) --native-pair-hmm-use-double-precision true 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/Bowtie2/%.g.vcf.gz: $(REFERENCE_FASTA) Output/WES/%.Bowtie2.Sort.MarkDuplicate.BQSR.bam $(REFERENCE_DBSNP_VCF) | Output/HaplotypeCaller/Bowtie2 build.log
	$(DOCKER) $(GATK) HaplotypeCaller --reference $(addprefix /,$(filter %.fasta,$^)) --input $(addprefix /,$(filter %.bam,$^)) --output $(addprefix /,$@) --emit-ref-confidence GVCF --dbsnp $(addprefix /,$(filter %.vcf,$^)) --native-pair-hmm-threads $(CPUS) --native-pair-hmm-use-double-precision true 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/BWA/HaplotypeCaller.genomicsdb/callset.json: $(REFERENCE_FASTA) $(addprefix Output/HaplotypeCaller/BWA/,$(addsuffix .g.vcf.gz,$(WES_TUMOR_ID))) | build.log
	$(DOCKER) $(GATK) GenomicsDBImport --reference $(addprefix /,$(filter %.fasta,$^)) --genomicsdb-workspace-path $(addprefix /,$(@D)) $(foreach sample,$(addprefix /,$(filter %.g.vcf.gz,$^)),--variant $(sample)) --reader-threads $(CPUS) --max-num-intervals-to-import-in-parallel $(CPUS) --overwrite-existing-genomicsdb-workspace true --merge-input-intervals true $(foreach chr,$(CHROMOSOMES),--intervals $(chr)) 1> $(@D).stdout 2> $(@D).stderr
	$(DOCKER) chmod --recursive 777 $(addprefix /,$(@D))

Output/HaplotypeCaller/Bowtie2/HaplotypeCaller.genomicsdb/callset.json: $(REFERENCE_FASTA) $(addprefix Output/HaplotypeCaller/Bowtie2/,$(addsuffix .g.vcf.gz,$(WES_TUMOR_ID))) | build.log
	$(DOCKER) $(GATK) GenomicsDBImport --reference $(addprefix /,$(filter %.fasta,$^)) --genomicsdb-workspace-path $(addprefix /,$(@D)) $(foreach sample,$(addprefix /,$(filter %.g.vcf.gz,$^)),--variant $(sample)) --reader-threads $(CPUS) --max-num-intervals-to-import-in-parallel $(CPUS) --overwrite-existing-genomicsdb-workspace true --merge-input-intervals true $(foreach chr,$(CHROMOSOMES),--intervals $(chr)) 1> $(@D).stdout 2> $(@D).stderr
	$(DOCKER) chmod --recursive 777 $(addprefix /,$(@D))

Output/HaplotypeCaller/%/HaplotypeCaller.vcf.gz: $(REFERENCE_FASTA) Output/HaplotypeCaller/%/HaplotypeCaller.genomicsdb/callset.json $(REFERENCE_DBSNP_VCF) | build.log
	$(DOCKER) $(GATK) GenotypeGVCFs --reference $(addprefix /,$(filter %.fasta,$^)) --variant "gendb://$(addprefix /,$(dir $(filter %.json,$^)))" --output $(addprefix /,$@) --dbsnp $(addprefix /,$(filter %.vcf,$^)) 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/%/HaplotypeCaller.VQSR.recal Output/HaplotypeCaller/%/HaplotypeCaller.VQSR.tranches: $(REFERENCE_FASTA) $(REFERENCE_HAPMAP_VCF) $(REFERENCE_OMNI_VCF) $(REFERENCE_1000G_VCF) $(REFERENCE_DBSNP_VCF) Output/HaplotypeCaller/%/HaplotypeCaller.vcf.gz | build.log
	$(DOCKER) $(GATK) VariantRecalibrator --reference $(addprefix /,$(filter %.fasta,$^)) --output $(addprefix /,$(basename $@).recal) --tranches-file $(addprefix /,$(basename $@).tranches) --variant $(addprefix /,$(lastword $^)) --use-annotation MQ --use-annotation MQRankSum --use-annotation ReadPosRankSum --use-annotation FS --use-annotation SOR --mode BOTH --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $(addprefix /,$(REFERENCE_HAPMAP_VCF)) --resource:omni,known=false,training=true,truth=false,prior=12.0 $(addprefix /,$(REFERENCE_OMNI_VCF)) --resource:1000G,known=false,training=true,truth=false,prior=10.0 $(addprefix /,$(REFERENCE_1000G_VCF)) --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $(addprefix /,$(REFERENCE_DBSNP_VCF)) 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/%/HaplotypeCaller.VQSR.vcf: $(REFERENCE_FASTA) Output/HaplotypeCaller/%/HaplotypeCaller.vcf.gz Output/HaplotypeCaller/%/HaplotypeCaller.VQSR.recal Output/HaplotypeCaller/%/HaplotypeCaller.VQSR.tranches | build.log
	$(DOCKER) $(GATK) ApplyVQSR --reference $(addprefix /,$(filter %.fasta,$^)) --variant $(addprefix /,$(filter %.vcf.gz,$^)) --recal-file $(addprefix /,$(filter %.recal,$^)) --tranches-file $(addprefix /,$(filter %.tranches,$^)) --output $(addprefix /,$@) --truth-sensitivity-filter-level 99.0 --mode BOTH 1> $@.stdout 2> $@.stderr

step15-BWA: $(addprefix Output/HaplotypeCaller/BWA/,$(addsuffix .g.vcf.gz,$(WES_TUMOR_ID))) Output/HaplotypeCaller/BWA/HaplotypeCaller.genomicsdb/callset.json Output/HaplotypeCaller/BWA/HaplotypeCaller.vcf.gz Output/HaplotypeCaller/BWA/HaplotypeCaller.VQSR.recal Output/HaplotypeCaller/BWA/HaplotypeCaller.VQSR.vcf
step15-Bowtie2: $(addprefix Output/HaplotypeCaller/Bowtie2/,$(addsuffix .g.vcf.gz,$(WES_TUMOR_ID))) Output/HaplotypeCaller/Bowtie2/HaplotypeCaller.genomicsdb/callset.json Output/HaplotypeCaller/Bowtie2/HaplotypeCaller.vcf.gz Output/HaplotypeCaller/Bowtie2/HaplotypeCaller.VQSR.recal Output/HaplotypeCaller/Bowtie2/HaplotypeCaller.VQSR.vcf
step15: step15-BWA step15-Bowtie2
.PHONY: step15-BWA step15-Bowtie2 step15
# TODO: wallace >> /BiO/Research/UNIST-Pancreas_Liver_Metastasis-SMC-2018-05/Workspace/seunghoon/hg38/haplotypecaller

# Step 16 (MutEnricher)
Output/MutEnricher Output/MutEnricher/BWA Output/MutEnricher/BWA/ADC Output/MutEnricher/BWA/SQC Output/MutEnricher/Bowtie2 Output/MutEnricher/Bowtie2/ADC Output/MutEnricher/Bowtie2/SQC:
	mkdir -p $@

Output/MutEnricher/BWA/%.hg38_multianno.vcf: Output/Mutect2/BWA/%.BWA.PASS.vcf | build.log Output/MutEnricher/BWA
	$(DOCKER) table_annovar.pl $(addprefix /,$^) /Tools/annovar/hg38 --verbose --buildver "hg38" --remove --vcfinput --protocol "refGene" --operation g --nastring "." --outfile $(addprefix /,$(@D)/$*) --thread $(CPUS) --maxgenethread $(CPUS) --tempdir $(TMPFS) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/BWA/%.vcf.gz: Bash/prepare_VCF_annotation.bash Output/MutEnricher/BWA/%.hg38_multianno.vcf | build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@)

Output/MutEnricher/BWA/%.vcf.gz.tbi: Output/MutEnricher/BWA/%.vcf.gz | build.log
	$(DOCKER) tabix --preset "vcf" --force $(addprefix /,$^) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/BWA/SQC.tsv: Python/make_list_MutEnricher.py $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz,$(SQC_WES_SAMPLES))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz.tbi,$(SQC_WES_SAMPLES))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$< $(filter %.vcf.gz,$^) $@) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/BWA/ADC.tsv: Python/make_list_MutEnricher.py $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz,$(ADC_WES_SAMPLES))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz.tbi,$(ADC_WES_SAMPLES))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$< $(filter %.vcf.gz,$^) $@) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/Bowtie2/%.hg38_multianno.vcf: Output/Mutect2/Bowtie2/%.Bowtie2.PASS.vcf | build.log Output/MutEnricher/BWA
	$(DOCKER) table_annovar.pl $(addprefix /,$^) /Tools/annovar/hg38 --verbose --buildver "hg38" --remove --vcfinput --protocol "refGene" --operation g --nastring "." --outfile $(addprefix /,$(@D)/$*) --thread $(CPUS) --maxgenethread $(CPUS) --tempdir $(TMPFS) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/Bowtie2/%.vcf.gz: Bash/prepare_VCF_annotation.bash Output/MutEnricher/Bowtie2/%.hg38_multianno.vcf | build.log Output/MutEnricher/Bowtie2
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/Bowtie2/%.vcf.gz.tbi: Output/MutEnricher/Bowtie2/%.vcf.gz | build.log
	$(DOCKER) tabix --preset "vcf" --force $(addprefix /,$^) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/Bowtie2/SQC.tsv: Python/make_list_MutEnricher.py $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .vcf.gz,$(SQC_WES_SAMPLES))) $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .vcf.gz.tbi,$(SQC_WES_SAMPLES))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$< $(filter %.vcf.gz,$^) $@) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/Bowtie2/ADC.tsv: Python/make_list_MutEnricher.py $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .vcf.gz,$(ADC_WES_SAMPLES))) $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .vcf.gz.tbi,$(ADC_WES_SAMPLES))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$< $(filter %.vcf.gz,$^) $@) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/gene_covariates.txt: $(KNOWN_GENE_GTF) $(REFERENCE_FASTA) | build.log
	$(DOCKER) $(PYTHON) /Tools/MutEnricher-1.3.1/utilities/get_gene_covariates.py --outname $(addprefix /,$@) --processors $(CPUS) $(addprefix /,$^) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/BWA/%/MutEnricher_gene_enrichments.txt Output/MutEnricher/BWA/%/MutEnricher_hotspot.txt Output/MutEnricher/BWA/%/MutEnricher_gene_hotspot_Fisher_enrichments.txt Output/MutEnricher/BWA/%/MutEnricher_gene_data.pkl Output/MutEnricher/BWA/%/MutEnricher.log: Output/MutEnricher/gene_covariates.txt $(KNOWN_GENE_GTF) Output/MutEnricher/BWA/%.tsv | build.log Output/MutEnricher/BWA/%
	$(DOCKER) $(MUTENRICHER) coding --outdir $(addprefix /,$(@D)) --prefix "MutEnricher" --anno-type "annovar-refGene" --exome-only --by-contig --processors $(CPUS) --covariates-file $(addprefix /,$(filter %.txt,$^)) $(addprefix /,$(filter %.gtf.gz,$^)) $(addprefix /,$(filter %.tsv,$^)) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/Bowtie2/%/MutEnricher_gene_enrichments.txt Output/MutEnricher/Bowtie2/%/MutEnricher_hotspot.txt Output/MutEnricher/Bowtie2/%/MutEnricher_gene_hotspot_Fisher_enrichments.txt Output/MutEnricher/Bowtie2/%/MutEnricher_gene_data.pkl Output/MutEnricher/Bowtie2/%/MutEnricher.log: Output/MutEnricher/gene_covariates.txt $(KNOWN_GENE_GTF) Output/MutEnricher/Bowtie2/%.tsv | build.log Output/MutEnricher/Bowtie2/%
	$(DOCKER) $(MUTENRICHER) coding --outdir $(addprefix /,$(@D)) --prefix "MutEnricher" --anno-type "annovar-refGene" --exome-only --by-contig --processors $(CPUS) --covariates-file $(addprefix /,$(filter %.txt,$^)) $(addprefix /,$(filter %.gtf.gz,$^)) $(addprefix /,$(filter %.tsv,$^)) 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-SQC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(SQC_WES_SAMPLES))) Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-SQC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(SQC_WES_SAMPLES))) Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-ADC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(ADC_WES_SAMPLES))) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-ADC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(ADC_WES_SAMPLES))) Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step16-preprocess: Output/MutEnricher/gene_covariates.txt $(addprefix Output/MutEnricher/BWA/,$(addsuffix .hg38_multianno.vcf,$(SQC_WES_SAMPLES))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz,$(SQC_WES_SAMPLES))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz.tbi,$(SQC_WES_SAMPLES))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .hg38_multianno.vcf,$(ADC_WES_SAMPLES))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz,$(ADC_WES_SAMPLES))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz.tbi,$(ADC_WES_SAMPLES))) Output/MutEnricher/BWA/SQC.tsv Output/MutEnricher/BWA/ADC.tsv $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .hg38_multianno.vcf,$(SQC_WES_SAMPLES))) $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .vcf.gz,$(SQC_WES_SAMPLES))) $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .hg38_multianno.vcf,$(ADC_WES_SAMPLES))) $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .vcf.gz.tbi,$(SQC_WES_SAMPLES))) $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .vcf.gz,$(ADC_WES_SAMPLES))) $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .vcf.gz.tbi,$(ADC_WES_SAMPLES))) Output/MutEnricher/Bowtie2/SQC.tsv Output/MutEnricher/Bowtie2/ADC.tsv
step16-MutEnricher: Output/MutEnricher/BWA/SQC/MutEnricher.log Output/MutEnricher/BWA/ADC/MutEnricher.log Output/MutEnricher/Bowtie2/ADC/MutEnricher.log Output/MutEnricher/Bowtie2/SQC/MutEnricher.log
step16-plot: Output/Mutect2/BWA-SQC.pdf Output/Mutect2/BWA-ADC.pdf Output/Mutect2/Bowtie2-SQC.pdf Output/Mutect2/Bowtie2-ADC.pdf
step16: step16-preprocess step16-MutEnricher step16-plot
.PHONY: step16-preprocess step16-MutEnricher step16-plot step16

# Step 17 (RSEM: Prepare reference genome)
Output/RSEM:
	mkdir -p $@

Output/RSEM/hg38.gtf.gz: | Output/RSEM
	wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz" -O $@
	touch $@

Output/RSEM/hg38.gtf: Output/RSEM/hg38.gtf.gz | build.log
	$(DOCKER) gzip --decompress --keep --force $(addprefix /,$<) 1> $@.stdout 2> $@.stderr

Output/RSEM/hg38.seq Output/RSEM/hg38.grp: $(REFERENCE_FASTA) Output/RSEM/hg38.gtf | build.log
	$(DOCKER) rsem-prepare-reference --gtf $(addprefix /,$(filter %.gtf,$^)) --bowtie2 --star --num-threads $(CPUS) $(addprefix /,$(filter %.fasta,$^)) /$(@D)/"hg38" 1> $@.stdout 2> $@.stderr

step17: Output/RSEM/hg38.gtf Output/RSEM/hg38.seq
.PHONY: step17

# Step 18 (RSEM: Calculate expression)
Output/RSEM/Bowtie2 Output/RSEM/STAR Output/RSEM/DEG:
	mkdir -p $@

Output/RSEM/Bowtie2/%.genome.bam Output/RSEM/Bowtie2/%.genes.results: Data/WTS/%_*.fastq.gz Output/RSEM/hg38.grp | Output/RSEM/Bowtie2 build.log
	$(DOCKER) rsem-calculate-expression --num-threads $(CPUS) --estimate-rspd --ci-memory $(shell echo "$(MEMS) * 1024" | bc) --bowtie2 --output-genome-bam --paired-end $(addprefix /,$(sort $(filter %.fastq.gz,$^)) $(basename $(filter %.grp,$^)) $(@D)/$*) 1> $@.stdout 2> $@.stderr

Output/RSEM/STAR/%.genome.bam Output/RSEM/STAR/%.genes.results: Data/WTS/%_*.fastq.gz Output/RSEM/hg38.grp | Output/RSEM/STAR build.log
	$(DOCKER) rsem-calculate-expression --num-threads $(CPUS) --estimate-rspd --ci-memory $(shell echo "$(MEMS) * 1024" | bc) --star --star-gzipped-read-file --output-genome-bam --paired-end $(addprefix /,$(sort $(filter %.fastq.gz,$^)) $(basename $(filter %.grp,$^)) $(@D)/$*) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/SQC.Bowtie2.Normal-Dysplasia.tsv Output/RSEM/DEG/SQC.Bowtie2.Normal-Dysplasia.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(SQC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Normal" "Dysplasia" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/SQC.Bowtie2.Normal-CIS.tsv Output/RSEM/DEG/SQC.Bowtie2.Normal-CIS.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(SQC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Normal" "CIS+AIS" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/SQC.Bowtie2.Normal-Primary.tsv Output/RSEM/DEG/SQC.Bowtie2.Normal-Primary.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(SQC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Normal" "Primary" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/SQC.Bowtie2.Dysplasia-CIS.tsv Output/RSEM/DEG/SQC.Bowtie2.Dysplasia-CIS.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(SQC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Dysplasia" "CIS+AIS" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/SQC.Bowtie2.Dysplasia-Primary.tsv Output/RSEM/DEG/SQC.Bowtie2.Dysplasia-Primary.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(SQC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Dysplasia" "Primary" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/SQC.Bowtie2.CIS-Primary.tsv Output/RSEM/DEG/SQC.Bowtie2.CIS-Primary.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(SQC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "CIS+AIS" "Primary" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.Bowtie2.Normal-AAH.tsv Output/RSEM/DEG/ADC.Bowtie2.Normal-AAH.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Normal" "AAH" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.Bowtie2.Normal-AIS.tsv Output/RSEM/DEG/ADC.Bowtie2.Normal-AIS.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Normal" "CIS+AIS" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.Bowtie2.Normal-MIA.tsv Output/RSEM/DEG/ADC.Bowtie2.Normal-MIA.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Normal" "MIA" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.Bowtie2.Normal-Primary.tsv Output/RSEM/DEG/ADC.Bowtie2.Normal-Primary.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Normal" "Primary" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.Bowtie2.AAH-AIS.tsv Output/RSEM/DEG/ADC.Bowtie2.AAH-CIS.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "AAH" "CIS+AIS" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.Bowtie2.AAH-MIA.tsv Output/RSEM/DEG/ADC.Bowtie2.AAH-MIA.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "AAH" "MIA" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.Bowtie2.AAH-Primary.tsv Output/RSEM/DEG/ADC.Bowtie2.AAH-Primary.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "AAH" "Primary" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.Bowtie2.AIS-MIA.tsv Output/RSEM/DEG/ADC.Bowtie2.AIS-MIA.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "CIS+AIS" "MIA" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.Bowtie2.AIS-Primary.tsv Output/RSEM/DEG/ADC.Bowtie2.AIS-Primary.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "CIS+AIS" "Primary" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.Bowtie2.MIA-Primary.tsv Output/RSEM/DEG/ADC.Bowtie2.MIA-Primary.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "MIA" "Primary" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/SQC.STAR.Normal-Dysplasia.tsv Output/RSEM/DEG/SQC.STAR.Normal-Dysplasia.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(SQC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Normal" "Dysplasia" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/SQC.STAR.Normal-CIS.tsv Output/RSEM/DEG/SQC.STAR.Normal-CIS.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(SQC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Normal" "CIS+AIS" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/SQC.STAR.Normal-Primary.tsv Output/RSEM/DEG/SQC.STAR.Normal-Primary.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(SQC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Normal" "Primary" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/SQC.STAR.Dysplasia-CIS.tsv Output/RSEM/DEG/SQC.STAR.Dysplasia-CIS.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(SQC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Dysplasia" "CIS+AIS" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/SQC.STAR.Dysplasia-Primary.tsv Output/RSEM/DEG/SQC.STAR.Dysplasia-Primary.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(SQC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Dysplasia" "Primary" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/SQC.STAR.CIS-Primary.tsv Output/RSEM/DEG/SQC.STAR.CIS-Primary.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(SQC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "CIS+AIS" "Primary" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.STAR.Normal-AAH.tsv Output/RSEM/DEG/ADC.STAR.Normal-AAH.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Normal" "AAH" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.STAR.Normal-AIS.tsv Output/RSEM/DEG/ADC.STAR.Normal-AIS.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Normal" "CIS+AIS" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.STAR.Normal-MIA.tsv Output/RSEM/DEG/ADC.STAR.Normal-MIA.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Normal" "MIA" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.STAR.Normal-Primary.tsv Output/RSEM/DEG/ADC.STAR.Normal-Primary.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "Normal" "Primary" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.STAR.AAH-AIS.tsv Output/RSEM/DEG/ADC.STAR.AAH-AIS.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "AAH" "CIS+AIS" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.STAR.AAH-MIA.tsv Output/RSEM/DEG/ADC.STAR.AAH-MIA.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "AAH" "MIA" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.STAR.AAH-Primary.tsv Output/RSEM/DEG/ADC.STAR.AAH-Primary.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "AAH" "Primary" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.STAR.AIS-MIA.tsv Output/RSEM/DEG/ADC.STAR.AIS-MIA.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "CIS+AIS" "MIA" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.STAR.AIS-Primary.tsv Output/RSEM/DEG/ADC.STAR.AIS-Primary.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "CIS+AIS" "Primary" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.STAR.MIA-Primary.tsv Output/RSEM/DEG/ADC.STAR.MIA-Primary.coldata: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --target "MIA" "Primary" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.DEG.tsv: R/run_DESeq.R Output/RSEM/DEG/%.tsv Output/RSEM/DEG/%.coldata | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --readcount $(addprefix /,$(filter %.tsv,$^)) --coldata $(addprefix /,$(filter %.coldata,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.Bowtie2.Aggregate.tsv: Python/aggregate_RSEM2.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/SQC.Bowtie2.Aggregate.tsv: Python/aggregate_RSEM2.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(SQC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/ADC.STAR.Aggregate.tsv: Python/aggregate_RSEM2.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/SQC.STAR.Aggregate.tsv: Python/aggregate_RSEM2.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(SQC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step18-preparation: $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(WTS_ID))) $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(WTS_ID)))
step18-Bowtie2: Output/RSEM/DEG/SQC.Bowtie2.Normal-Dysplasia.tsv Output/RSEM/DEG/SQC.Bowtie2.Normal-CIS.tsv Output/RSEM/DEG/SQC.Bowtie2.Normal-Primary.tsv Output/RSEM/DEG/SQC.Bowtie2.Dysplasia-CIS.tsv Output/RSEM/DEG/SQC.Bowtie2.Dysplasia-Primary.tsv Output/RSEM/DEG/SQC.Bowtie2.CIS-Primary.tsv Output/RSEM/DEG/ADC.Bowtie2.Normal-AAH.tsv Output/RSEM/DEG/ADC.Bowtie2.Normal-AIS.tsv Output/RSEM/DEG/ADC.Bowtie2.Normal-MIA.tsv Output/RSEM/DEG/ADC.Bowtie2.Normal-Primary.tsv Output/RSEM/DEG/ADC.Bowtie2.AAH-AIS.tsv Output/RSEM/DEG/ADC.Bowtie2.AAH-MIA.tsv Output/RSEM/DEG/ADC.Bowtie2.AAH-Primary.tsv Output/RSEM/DEG/ADC.Bowtie2.AIS-MIA.tsv Output/RSEM/DEG/ADC.Bowtie2.AIS-Primary.tsv Output/RSEM/DEG/ADC.Bowtie2.MIA-Primary.tsv
step18-Bowtie2-DEG: Output/RSEM/DEG/SQC.Bowtie2.Normal-Dysplasia.DEG.tsv Output/RSEM/DEG/SQC.Bowtie2.Normal-CIS.DEG.tsv Output/RSEM/DEG/SQC.Bowtie2.Normal-Primary.DEG.tsv Output/RSEM/DEG/SQC.Bowtie2.Dysplasia-CIS.DEG.tsv Output/RSEM/DEG/SQC.Bowtie2.Dysplasia-Primary.DEG.tsv Output/RSEM/DEG/SQC.Bowtie2.CIS-Primary.DEG.tsv Output/RSEM/DEG/ADC.Bowtie2.Normal-AAH.DEG.tsv Output/RSEM/DEG/ADC.Bowtie2.Normal-AIS.DEG.tsv Output/RSEM/DEG/ADC.Bowtie2.Normal-Primary.DEG.tsv Output/RSEM/DEG/ADC.Bowtie2.AAH-AIS.DEG.tsv Output/RSEM/DEG/ADC.Bowtie2.AAH-Primary.DEG.tsv Output/RSEM/DEG/ADC.Bowtie2.AIS-Primary.DEG.tsv
step18-STAR: Output/RSEM/DEG/SQC.STAR.Normal-Dysplasia.tsv Output/RSEM/DEG/SQC.STAR.Normal-CIS.tsv Output/RSEM/DEG/SQC.STAR.Normal-Primary.tsv Output/RSEM/DEG/SQC.STAR.Dysplasia-CIS.tsv Output/RSEM/DEG/SQC.STAR.Dysplasia-Primary.tsv Output/RSEM/DEG/SQC.STAR.CIS-Primary.tsv Output/RSEM/DEG/ADC.STAR.Normal-AAH.tsv Output/RSEM/DEG/ADC.STAR.Normal-AIS.tsv Output/RSEM/DEG/ADC.STAR.Normal-MIA.tsv Output/RSEM/DEG/ADC.STAR.Normal-Primary.tsv Output/RSEM/DEG/ADC.STAR.AAH-AIS.tsv Output/RSEM/DEG/ADC.STAR.AAH-MIA.tsv Output/RSEM/DEG/ADC.STAR.AAH-Primary.tsv Output/RSEM/DEG/ADC.STAR.AIS-MIA.tsv Output/RSEM/DEG/ADC.STAR.AIS-Primary.tsv Output/RSEM/DEG/ADC.STAR.MIA-Primary.tsv
step18-STAR-DEG: Output/RSEM/DEG/SQC.STAR.Normal-Dysplasia.DEG.tsv Output/RSEM/DEG/SQC.STAR.Normal-CIS.DEG.tsv Output/RSEM/DEG/SQC.STAR.Normal-Primary.DEG.tsv Output/RSEM/DEG/SQC.STAR.Dysplasia-CIS.DEG.tsv Output/RSEM/DEG/SQC.STAR.Dysplasia-Primary.DEG.tsv Output/RSEM/DEG/SQC.STAR.CIS-Primary.DEG.tsv Output/RSEM/DEG/ADC.STAR.Normal-AAH.DEG.tsv Output/RSEM/DEG/ADC.STAR.Normal-AIS.DEG.tsv Output/RSEM/DEG/ADC.STAR.Normal-Primary.DEG.tsv Output/RSEM/DEG/ADC.STAR.AAH-AIS.DEG.tsv Output/RSEM/DEG/ADC.STAR.AAH-Primary.DEG.tsv Output/RSEM/DEG/ADC.STAR.AIS-Primary.DEG.tsv
step18: step18-preparation step18-Bowtie2 step18-Bowtie-DEG step18-STAR step18-STAR-DEG
.PHONY: step18 step18-preparation step18-Bowtie2 step18-Bowtie-DEG step18-STAR step18-STAR-DEG

# Step 19 (DEG: Draw plots)
Output/RSEM/Plots:
	mkdir -p $@

Output/RSEM/Plots/ADC.Up.%.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/ADC.%.Normal-AAH.DEG.tsv Output/RSEM/DEG/ADC.%.Normal-AIS.DEG.tsv Output/RSEM/DEG/ADC.%.Normal-Primary.DEG.tsv Output/RSEM/DEG/ADC.%.AAH-AIS.DEG.tsv Output/RSEM/DEG/ADC.%.AAH-Primary.DEG.tsv Output/RSEM/DEG/ADC.%.AIS-Primary.DEG.tsv | Output/RSEM/Plots build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/ADC.Down.%.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/ADC.%.Normal-AAH.DEG.tsv Output/RSEM/DEG/ADC.%.Normal-AIS.DEG.tsv Output/RSEM/DEG/ADC.%.Normal-Primary.DEG.tsv Output/RSEM/DEG/ADC.%.AAH-AIS.DEG.tsv Output/RSEM/DEG/ADC.%.AAH-Primary.DEG.tsv Output/RSEM/DEG/ADC.%.AIS-Primary.DEG.tsv | Output/RSEM/Plots build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/SQC.Up.%.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/SQC.%.Normal-Dysplasia.DEG.tsv Output/RSEM/DEG/SQC.%.Normal-CIS.DEG.tsv Output/RSEM/DEG/SQC.%.Normal-Primary.DEG.tsv Output/RSEM/DEG/SQC.%.Dysplasia-CIS.DEG.tsv Output/RSEM/DEG/SQC.%.Dysplasia-Primary.DEG.tsv Output/RSEM/DEG/SQC.%.CIS-Primary.DEG.tsv | Output/RSEM/Plots build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/SQC.Down.%.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/SQC.%.Normal-Dysplasia.DEG.tsv Output/RSEM/DEG/SQC.%.Normal-CIS.DEG.tsv Output/RSEM/DEG/SQC.%.Normal-Primary.DEG.tsv Output/RSEM/DEG/SQC.%.Dysplasia-CIS.DEG.tsv Output/RSEM/DEG/SQC.%.Dysplasia-Primary.DEG.tsv Output/RSEM/DEG/SQC.%.CIS-Primary.DEG.tsv | Output/RSEM/Plots build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/ADC.Up.%.pair-venn.pdf: Python/draw_venn_plot_2.py Output/RSEM/DEG/ADC.%.Normal-AAH.DEG.tsv Output/RSEM/DEG/ADC.%.Normal-AIS.DEG.tsv Output/RSEM/DEG/ADC.%.Normal-Primary.DEG.tsv | Output/RSEM/Plots build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/ADC.Down.%.pair-venn.pdf: Python/draw_venn_plot_2.py Output/RSEM/DEG/ADC.%.Normal-AAH.DEG.tsv Output/RSEM/DEG/ADC.%.Normal-AIS.DEG.tsv Output/RSEM/DEG/ADC.%.Normal-Primary.DEG.tsv | Output/RSEM/Plots build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/SQC.Up.%.pair-venn.pdf: Python/draw_venn_plot_2.py Output/RSEM/DEG/SQC.%.Normal-Dysplasia.DEG.tsv Output/RSEM/DEG/SQC.%.Normal-CIS.DEG.tsv Output/RSEM/DEG/SQC.%.Normal-Primary.DEG.tsv | Output/RSEM/Plots build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/SQC.Down.%.pair-venn.pdf: Python/draw_venn_plot_2.py Output/RSEM/DEG/SQC.%.Normal-Dysplasia.DEG.tsv Output/RSEM/DEG/SQC.%.Normal-CIS.DEG.tsv Output/RSEM/DEG/SQC.%.Normal-Primary.DEG.tsv | Output/RSEM/Plots build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/ADC.%.volcano.pdf: Python/draw_volcano_plot.py Output/RSEM/DEG/ADC.%.DEG.tsv | Output/RSEM/Plots build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/SQC.%.volcano.pdf: Python/draw_volcano_plot.py Output/RSEM/DEG/SQC.%.DEG.tsv | Output/RSEM/Plots build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/SQC.%.violin.tar: Python/draw_violin_plots.py Output/RSEM/DEG/SQC.%.Normal-Dysplasia.tsv Output/RSEM/DEG/SQC.%.Normal-CIS.tsv Output/RSEM/DEG/SQC.%.Normal-Primary.tsv Output/RSEM/DEG/SQC.%.Dysplasia-CIS.tsv Output/RSEM/DEG/SQC.%.Dysplasia-Primary.tsv Output/RSEM/DEG/SQC.%.CIS-Primary.tsv | Output/RSEM/Plots build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/ADC.%.violin.tar: Python/draw_violin_plots.py Output/RSEM/DEG/ADC.%.Normal-AAH.tsv Output/RSEM/DEG/ADC.%.Normal-AIS.tsv Output/RSEM/DEG/ADC.%.Normal-MIA.tsv Output/RSEM/DEG/ADC.%.Normal-Primary.tsv Output/RSEM/DEG/ADC.%.AAH-AIS.tsv Output/RSEM/DEG/ADC.%.AAH-MIA.tsv Output/RSEM/DEG/ADC.%.AAH-Primary.tsv Output/RSEM/DEG/ADC.%.AIS-MIA.tsv Output/RSEM/DEG/ADC.%.AIS-Primary.tsv Output/RSEM/DEG/ADC.%.MIA-Primary.tsv | Output/RSEM/Plots build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/SQC.%.pair.pdf: Python/draw_pair_plot.py Output/RSEM/DEG/SQC.%.tsv Output/RSEM/DEG/SQC.%.DEG.tsv | Output/RSEM/Plots build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC

Output/RSEM/Plots/ADC.%.pair.pdf: Python/draw_pair_plot.py Output/RSEM/DEG/ADC.%.tsv Output/RSEM/DEG/ADC.%.DEG.tsv | Output/RSEM/Plots build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC

step19-venn: Output/RSEM/Plots/ADC.Up.Bowtie2.venn.pdf Output/RSEM/Plots/ADC.Down.Bowtie2.venn.pdf Output/RSEM/Plots/SQC.Up.Bowtie2.venn.pdf Output/RSEM/Plots/SQC.Down.Bowtie2.venn.pdf Output/RSEM/Plots/ADC.Up.STAR.venn.pdf Output/RSEM/Plots/ADC.Down.STAR.venn.pdf Output/RSEM/Plots/SQC.Up.STAR.venn.pdf Output/RSEM/Plots/SQC.Down.STAR.venn.pdf
step19-pair-venn: Output/RSEM/Plots/ADC.Up.Bowtie2.pair-venn.pdf Output/RSEM/Plots/ADC.Down.Bowtie2.pair-venn.pdf Output/RSEM/Plots/SQC.Up.Bowtie2.pair-venn.pdf Output/RSEM/Plots/SQC.Down.Bowtie2.pair-venn.pdf Output/RSEM/Plots/ADC.Up.STAR.pair-venn.pdf Output/RSEM/Plots/ADC.Down.STAR.pair-venn.pdf Output/RSEM/Plots/SQC.Up.STAR.pair-venn.pdf Output/RSEM/Plots/SQC.Down.STAR.pair-venn.pdf
step19-ADC-Bowtie2-volcano: Output/RSEM/Plots/ADC.Bowtie2.Normal-AAH.volcano.pdf Output/RSEM/Plots/ADC.Bowtie2.Normal-AIS.volcano.pdf Output/RSEM/Plots/ADC.Bowtie2.Normal-Primary.volcano.pdf Output/RSEM/Plots/ADC.Bowtie2.AAH-AIS.volcano.pdf Output/RSEM/Plots/ADC.Bowtie2.AAH-Primary.volcano.pdf Output/RSEM/Plots/ADC.Bowtie2.AIS-Primary.volcano.pdf Output/RSEM/Plots/ADC.STAR.Normal-AAH.volcano.pdf Output/RSEM/Plots/ADC.STAR.Normal-AIS.volcano.pdf Output/RSEM/Plots/ADC.STAR.Normal-Primary.volcano.pdf Output/RSEM/Plots/ADC.STAR.AAH-AIS.volcano.pdf Output/RSEM/Plots/ADC.STAR.AAH-Primary.volcano.pdf Output/RSEM/Plots/ADC.STAR.AIS-Primary.volcano.pdf
step19-SQC-Bowtie2-volcano: Output/RSEM/Plots/SQC.Bowtie2.Normal-Dysplasia.volcano.pdf Output/RSEM/Plots/SQC.Bowtie2.Normal-CIS.volcano.pdf Output/RSEM/Plots/SQC.Bowtie2.Normal-Primary.volcano.pdf Output/RSEM/Plots/SQC.Bowtie2.Dysplasia-CIS.volcano.pdf Output/RSEM/Plots/SQC.Bowtie2.Dysplasia-Primary.volcano.pdf Output/RSEM/Plots/SQC.Bowtie2.CIS-Primary.volcano.pdf Output/RSEM/Plots/SQC.STAR.Normal-Dysplasia.volcano.pdf Output/RSEM/Plots/SQC.STAR.Normal-CIS.volcano.pdf Output/RSEM/Plots/SQC.STAR.Normal-Primary.volcano.pdf Output/RSEM/Plots/SQC.STAR.Dysplasia-CIS.volcano.pdf Output/RSEM/Plots/SQC.STAR.Dysplasia-Primary.volcano.pdf Output/RSEM/Plots/SQC.STAR.CIS-Primary.volcano.pdf
step19-violin: Output/RSEM/Plots/SQC.Bowtie2.violin.tar Output/RSEM/Plots/SQC.STAR.violin.tar Output/RSEM/Plots/ADC.Bowtie2.violin.tar Output/RSEM/Plots/ADC.STAR.violin.tar Output/RSEM/Plots/SQC.STAR.violin.tar Output/RSEM/Plots/SQC.STAR.violin.tar Output/RSEM/Plots/ADC.STAR.violin.tar Output/RSEM/Plots/ADC.STAR.violin.tar
step19: step19-venn step19-pair-venn step19-ADC-Bowtie2-volcano step19-SQC-Bowtie2-volcano step19-violin
.PHONY: step19 step19-venn step19-pair-venn step19-ADC-Bowtie2-volcano step19-SQC-Bowtie2-volcano step19-violin

# Step 20 (Get ready to run MesKit)
Output/MesKit:
	mkdir -p $@

Output/MesKit/SQC_clinical.tsv: Python/make_clinical_MesKit.py | Output/MesKit build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^) $(SQC_WES_SAMPLES) ".BWA.sam" $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/MesKit/ADC_clinical.tsv: Python/make_clinical_MesKit.py | Output/MesKit build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^) $(ADC_WES_SAMPLES) ".BWA.sam" $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/MesKit/%_gene.tsv: Python/select_gene_MesKit.py Output/Witer/%.txt | Output/MesKit build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@)

step20: Output/MesKit/SQC_clinical.tsv Output/MesKit/ADC_clinical.tsv Output/MesKit/SQC_gene.tsv Output/MesKit/ADC_gene.tsv
.PHONY: step20

# Step 21 (Run PyClone)
Output/PyClone/BWA Output/PyClone/Bowtie2 Output/PyClone/BWA/Plots Output/PyClone/Bowtie2/Plots:
	mkdir -p $@

Output/PyClone/BWA/%.tsv Output/PyClone/BWA/%.Cellularity.txt: Python/make_input_pyclone.py Output/Mutect2/BWA/%.BWA.PASS.maf Output/Sequenza/BWA/%.tar.gz | Output/PyClone/BWA build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(@D)/$*) --cpus $(CPUS) 1> $(@D)/$*.preprocess.stdout 2> $(@D)/$*.preprocess.stderr

Output/PyClone/Bowtie2/%.tsv Output/PyClone/Bowtie2/%.Cellularity.txt: Python/make_input_pyclone.py Output/Mutect2/Bowtie2/%.Bowtie2.PASS.maf Output/Sequenza/Bowtie2/%.tar.gz | Output/PyClone/Bowtie2 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(@D)/$*) --cpus $(CPUS) 1> $(@D)/$*.preprocess.stdout 2> $(@D)/$*.preprocess.stderr

# I know that this recipe contains unnecessary dependencies; but I cannot see the solution.
Output/PyClone/BWA/%/loci.tsv: Bash/run_pyclone.bash $(addprefix Output/PyClone/BWA/,$(addsuffix .tsv,$(SQC_WES_SAMPLES) $(ADC_WES_SAMPLES))) $(addprefix Output/PyClone/BWA/,$(addsuffix .Cellularity.txt,$(SQC_WES_SAMPLES) $(ADC_WES_SAMPLES))) | build.log
	$(DOCKER) /bin/bash $(addprefix /,$< Output/PyClone/BWA/$*.tsv Output/PyClone/BWA/$(PAIRED_PRIMARY).tsv Output/PyClone/BWA/$*.Cellularity.txt Output/PyClone/BWA/$(PAIRED_PRIMARY).Cellularity.txt $(@D)) 1> $(@D).stdout 2> $(@D).stderr || true

# I know that this recipe contains unnecessary dependencies; but I cannot see the solution.
Output/PyClone/Bowtie2/%/loci.tsv: Bash/run_pyclone.bash $(addprefix Output/PyClone/Bowtie2/,$(addsuffix .tsv,$(SQC_WES_SAMPLES) $(ADC_WES_SAMPLES))) $(addprefix Output/PyClone/Bowtie2/,$(addsuffix .Cellularity.txt,$(SQC_WES_SAMPLES) $(ADC_WES_SAMPLES))) | build.log
	$(DOCKER) /bin/bash $(addprefix /,$< Output/PyClone/Bowtie2/$*.tsv Output/PyClone/Bowtie2/$(PAIRED_PRIMARY).tsv Output/PyClone/Bowtie2/$*.Cellularity.txt Output/PyClone/Bowtie2/$(PAIRED_PRIMARY).Cellularity.txt $(@D)) 1> $(@D).stdout 2> $(@D).stderr || true

Output/PyClone/BWA/Plots/%.ADC.CCF.pdf: Python/draw_output_pyclone.py Output/PyClone/BWA/%/loci.tsv Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/PyClone/BWA/Plots
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --CCF 1> $@.stdout 2> $@.stderr

Output/PyClone/BWA/Plots/%.SQC.CCF.pdf: Python/draw_output_pyclone.py Output/PyClone/BWA/%/loci.tsv Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/PyClone/BWA/Plots
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --CCF 1> $@.stdout 2> $@.stderr

Output/PyClone/BWA/Plots/%.ADC.VAF.pdf: Python/draw_output_pyclone.py Output/PyClone/BWA/%/loci.tsv Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/PyClone/BWA/Plots
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --VAF 1> $@.stdout 2> $@.stderr

Output/PyClone/BWA/Plots/%.SQC.VAF.pdf: Python/draw_output_pyclone.py Output/PyClone/BWA/%/loci.tsv Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/PyClone/BWA/Plots
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --VAF 1> $@.stdout 2> $@.stderr

Output/PyClone/Bowtie2/Plots/%.ADC.CCF.pdf: Python/draw_output_pyclone.py Output/PyClone/Bowtie2/%/loci.tsv Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/PyClone/Bowtie2/Plots
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --CCF 1> $@.stdout 2> $@.stderr

Output/PyClone/Bowtie2/Plots/%.SQC.CCF.pdf: Python/draw_output_pyclone.py Output/PyClone/Bowtie2/%/loci.tsv Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/PyClone/Bowtie2/Plots
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --CCF 1> $@.stdout 2> $@.stderr

Output/PyClone/Bowtie2/Plots/%.ADC.VAF.pdf: Python/draw_output_pyclone.py Output/PyClone/Bowtie2/%/loci.tsv Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/PyClone/Bowtie2/Plots
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --VAF 1> $@.stdout 2> $@.stderr

Output/PyClone/Bowtie2/Plots/%.SQC.VAF.pdf: Python/draw_output_pyclone.py Output/PyClone/Bowtie2/%/loci.tsv Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/PyClone/Bowtie2/Plots
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --VAF 1> $@.stdout 2> $@.stderr

step21-preprocess: $(addprefix Output/PyClone/BWA/,$(addsuffix .tsv,$(WES_TUMOR_ID))) $(addprefix Output/PyClone/Bowtie2/,$(addsuffix .tsv,$(WES_TUMOR_ID)))
step21-pyclone: $(addprefix Output/PyClone/BWA/,$(addsuffix /loci.tsv,$(filter-out %P,$(filter-out %P1,$(SQC_WES_SAMPLES) $(ADC_WES_SAMPLES))))) $(addprefix Output/PyClone/Bowtie2/,$(addsuffix /loci.tsv,$(filter-out %P,$(filter-out %P1,$(SQC_WES_SAMPLES) $(ADC_WES_SAMPLES)))))
step21-plot: $(addprefix Output/PyClone/BWA/Plots/,$(addsuffix .ADC.CCF.pdf,$(filter-out %P,$(filter-out %P1,$(ADC_WES_SAMPLES))))) $(addprefix Output/PyClone/BWA/Plots/,$(addsuffix .SQC.CCF.pdf,$(filter-out %P,$(filter-out %P1,$(SQC_WES_SAMPLES))))) $(addprefix Output/PyClone/BWA/Plots/,$(addsuffix .ADC.VAF.pdf,$(filter-out %P,$(filter-out %P1,$(ADC_WES_SAMPLES))))) $(addprefix Output/PyClone/BWA/Plots/,$(addsuffix .SQC.VAF.pdf,$(filter-out %P,$(filter-out %P1,$(SQC_WES_SAMPLES))))) $(addprefix Output/PyClone/Bowtie2/Plots/,$(addsuffix .ADC.CCF.pdf,$(filter-out %P,$(filter-out %P1,$(ADC_WES_SAMPLES))))) $(addprefix Output/PyClone/Bowtie2/Plots/,$(addsuffix .SQC.CCF.pdf,$(filter-out %P,$(filter-out %P1,$(SQC_WES_SAMPLES))))) $(addprefix Output/PyClone/Bowtie2/Plots/,$(addsuffix .ADC.VAF.pdf,$(filter-out %P,$(filter-out %P1,$(ADC_WES_SAMPLES))))) $(addprefix Output/PyClone/Bowtie2/Plots/,$(addsuffix .SQC.VAF.pdf,$(filter-out %P,$(filter-out %P1,$(SQC_WES_SAMPLES)))))
step21: step21-preprocess step21-pyclone step21-plot
.PHONY: step21-preprocess step21-pyclone step21-plot step21

# Step 22 (Run CIBERSORTx)
Output/CIBERSORTx:
	mkdir -p $@

Output/CIBERSORTx/SQC.Bowtie2.tsv: Python/aggregate_RSEM_cibersort.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(SQC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/CIBERSORTx build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/ADC.Bowtie2.tsv: Python/aggregate_RSEM_cibersort.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/CIBERSORTx build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/SQC.STAR.tsv: Python/aggregate_RSEM_cibersort.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(SQC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/CIBERSORTx build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/ADC.STAR.tsv: Python/aggregate_RSEM_cibersort.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(ADC_WTS_SAMPLES))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/CIBERSORTx build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/GSE131907_signature.tsv: Python/make_signature_cibersort.py Output/CIBERSORTx/GSE131907_Lung_Cancer_cell_annotation.txt $(LUNG_SINGLECELL) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/%/CIBERSORTxGEP_NA_Weights.txt: Output/CIBERSORTx/%.tsv Output/CIBERSORTx/GSE131907_signature.tsv | build.log
	$(DOCKER) rm -rfv $(@D)
	docker run $(RUN_OPTS) --volume $(abspath $(dir $<)):/src/data --volume $(abspath $(@D)):/src/outdir cibersortx/hires $(CIBERSORTx_TOKEN) --mixture $(addprefix /,$(word 1,$^)) --sigmatrix $(addprefix /,$(word 2,$^)) --threads $(CPUS) --QN TRUE --heatmap TRUE --cluster TRUE 1> $(@D).stdout 2> $(@D).stderr || true

Output/CIBERSORTx/SQC.%.heatmap.pdf: Python/draw_heatmap_cibersort.py Output/CIBERSORTx/SQC.%/CIBERSORTxGEP_NA_Weights.txt | build
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/ADC.%.heatmap.pdf: Python/draw_heatmap_cibersort.py Output/CIBERSORTx/ADC.%/CIBERSORTxGEP_NA_Weights.txt | build
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/SQC.%.violin.tar: Python/draw_violin_plot_cibersort.py Output/CIBERSORTx/SQC.%/CIBERSORTxGEP_NA_Weights.txt | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/ADC.%.violin.tar: Python/draw_violin_plot_cibersort.py Output/CIBERSORTx/ADC.%/CIBERSORTxGEP_NA_Weights.txt | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/SQC.%.cluster.pdf: Python/draw_clustermap_cibersort.py Output/CIBERSORTx/SQC.%/CIBERSORTxGEP_NA_Weights.txt | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/ADC.%.cluster.pdf: Python/draw_clustermap_cibersort.py Output/CIBERSORTx/ADC.%/CIBERSORTxGEP_NA_Weights.txt | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

step22-RSEM: Output/CIBERSORTx/SQC.Bowtie2.tsv Output/CIBERSORTx/ADC.Bowtie2.tsv Output/CIBERSORTx/SQC.STAR.tsv Output/CIBERSORTx/ADC.STAR.tsv
step22-CIBERSORTx: Output/CIBERSORTx/GSE131907_signature.tsv Output/CIBERSORTx/SQC.Bowtie2/CIBERSORTxGEP_NA_Weights.txt Output/CIBERSORTx/ADC.Bowtie2/CIBERSORTxGEP_NA_Weights.txt Output/CIBERSORTx/SQC.STAR/CIBERSORTxGEP_NA_Weights.txt Output/CIBERSORTx/ADC.STAR/CIBERSORTxGEP_NA_Weights.txt
step22-plot: Output/CIBERSORTx/SQC.Bowtie2.heatmap.pdf Output/CIBERSORTx/ADC.Bowtie2.heatmap.pdf Output/CIBERSORTx/SQC.STAR.heatmap.pdf Output/CIBERSORTx/ADC.STAR.heatmap.pdf Output/CIBERSORTx/SQC.Bowtie2.violin.tar Output/CIBERSORTx/ADC.Bowtie2.violin.tar Output/CIBERSORTx/SQC.STAR.violin.tar Output/CIBERSORTx/ADC.STAR.violin.tar Output/CIBERSORTx/SQC.Bowtie2.cluster.pdf Output/CIBERSORTx/SQC.STAR.cluster.pdf Output/CIBERSORTx/ADC.Bowtie2.cluster.pdf Output/CIBERSORTx/ADC.STAR.cluster.pdf
step22: step22-RSEM step22-CIBERSORTx step22-plot
.PHONY: step22-RSEM step22-CIBERSORTx step22-plot step22

# Step 23 (VAF plot)
Output/VAF:
	mkdir -p $@

Output/VAF/SQC.BWA.tar: Python/draw_vaf_plots.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(SQC_WES_SAMPLES))) $(GENE_CENSUS_CSV) | Output/VAF build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/VAF/ADC.BWA.tar: Python/draw_vaf_plots.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(ADC_WES_SAMPLES))) $(GENE_CENSUS_CSV) | Output/VAF build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/VAF/SQC.Bowtie2.tar: Python/draw_vaf_plots.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(SQC_WES_SAMPLES))) $(GENE_CENSUS_CSV) | Output/VAF build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/VAF/ADC.Bowtie2.tar: Python/draw_vaf_plots.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(ADC_WES_SAMPLES))) $(GENE_CENSUS_CSV) | Output/VAF build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step23: Output/VAF/SQC.BWA.tar Output/VAF/ADC.BWA.tar Output/VAF/SQC.Bowtie2.tar Output/VAF/ADC.Bowtie2.tar
.PHONY: step23

# Step 24 (Revolver)
Output/Revolver/BWA Output/Revolver/Bowtie2:
	mkdir -p $@

Output/Revolver/BWA/SQC.tsv: Python/make_input_revolver.py $(addprefix Output/PyClone/BWA/,$(addsuffix /loci.tsv,$(filter-out %P,$(filter-out %P1,$(SQC_WES_SAMPLES))))) Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/Revolver/BWA
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Revolver/BWA/ADC.tsv: Python/make_input_revolver.py $(addprefix Output/PyClone/BWA/,$(addsuffix /loci.tsv,$(filter-out %P,$(filter-out %P1,$(ADC_WES_SAMPLES))))) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/Revolver/BWA
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Revolver/Bowtie2/SQC.tsv: Python/make_input_revolver.py $(addprefix Output/PyClone/Bowtie2/,$(addsuffix /loci.tsv,$(filter-out %P,$(filter-out %P1,$(SQC_WES_SAMPLES))))) Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/Revolver/Bowtie2
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Revolver/Bowtie2/ADC.tsv: Python/make_input_revolver.py $(addprefix Output/PyClone/Bowtie2/,$(addsuffix /loci.tsv,$(filter-out %P,$(filter-out %P1,$(ADC_WES_SAMPLES))))) Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/Revolver/Bowtie2
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Revolver/BWA/%.cluster.pdf Output/Revolver/BWA/%.driver.pdf Output/Revolver/BWA/%.dendrogram.pdf Output/Revolver/BWA/%.det.pdf Output/Revolver/BWA/%.clonality.pdf: R/run_revolver.R Output/Revolver/BWA/%.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.tsv,$^)) --cluster $(addprefix /,$(@D)/$*.cluster.pdf) --driver $(addprefix /,$(@D)/$*.driver.pdf) --dendrogram $(addprefix /,$(@D)/$*.dendrogram.pdf) --det $(addprefix /,$(@D)/$*.det.pdf) --clonality $(addprefix /,$(@D)/$*.clonality.pdf) 1> $@.stdout 2> $@.stderr

Output/Revolver/Bowtie2/%.cluster.pdf Output/Revolver/Bowtie2/%.driver.pdf Output/Revolver/Bowtie2/%.dendrogram.pdf Output/Revolver/Bowtie2/%.det.pdf Output/Revolver/Bowtie2/%.clonality.pdf: R/run_revolver.R Output/Revolver/Bowtie2/%.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.tsv,$^)) --cluster $(addprefix /,$(@D)/$*.cluster.pdf) --driver $(addprefix /,$(@D)/$*.driver.pdf) --dendrogram $(addprefix /,$(@D)/$*.dendrogram.pdf) --det $(addprefix /,$(@D)/$*.det.pdf) --clonality $(addprefix /,$(@D)/$*.clonality.pdf) 1> $@.stdout 2> $@.stderr

step24-input: Output/Revolver/BWA/SQC.tsv Output/Revolver/BWA/ADC.tsv Output/Revolver/Bowtie2/SQC.tsv Output/Revolver/Bowtie2/ADC.tsv
step24-revolver: Output/Revolver/BWA/SQC.cluster.pdf Output/Revolver/Bowtie2/SQC.cluster.pdf# Output/Revolver/BWA/ADC.cluster.pdf Output/Revolver/Bowtie2/ADC.cluster.pdf
step24: step24-input step24-revolver
.PHONY: step24-input step24-revolver step24

# Step 25 (Arriba)
Output/Arriba/Raw Output/Arriba/Filtered:
	mkdir -p $@

Output/Arriba/Raw/%/fusions.tsv: $(GENCODE_GTF) $(REFERENCE_FASTA) Data/WTS/%_*.fastq.gz | Output/Arriba/Raw build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D)) && mkdir -p $(@D)
	$(ARRIBA_DOCKER) run_arriba.sh $(addprefix /,$(ARRIBA_STAR_DIR) $(filter %.gtf,$^) $(filter %.fasta,$^) $(ARRIBA_BLACKLIST_TSV) $(ARRIBA_KNOWN_FUSIONS_TSV) $(ARRIBA_PROTEIN_DOMAIN_GFF3)) $(CPUS) $(addprefix /,$(sort $(filter %.fastq.gz,$^))) 1> $(@D).stdout 2> $(@D).stderr

# I know that this recipe contains unnecessary dependencies; but I cannot see the solution.
Output/Arriba/Filtered/%.tsv: Python/filter_arriba.py Output/Arriba/Raw/%/fusions.tsv $(addprefix Output/Arriba/Raw/,$(addsuffix /fusions.tsv,$(WTS_NORMAL_ID))) | build.log Output/Arriba/Filtered
	$(DOCKER) $(PYTHON) $(addprefix /,$< Output/Arriba/Raw/$*/fusions.tsv Output/Arriba/Raw/$(PAIRED_NORMAL)/fusions.tsv $@) 1> $@.stdout 2> $@.stderr

step25-Raw: $(addprefix Output/Arriba/Raw/,$(addsuffix /fusions.tsv,$(WTS_ID)))
step25-Filtered: $(addprefix Output/Arriba/Filtered/,$(addsuffix .tsv,$(WTS_TUMOR_ID)))
step25: step25-Raw
.PHONY: step25-Raw step25

# Step 26 (FusionCatcher)
Output/FusionCatcher/STAR Output/FusionCatcher/Bowtie2:
	mkdir -p $@

Output/FusionCatcher/STAR/%/preliminary-candidate-fusions-found-in-normal.log Output/FusionCatcher/STAR/%/partial-preliminary-candidate-fusions-found-in-normal.log: Data/WTS/%_*.fastq.gz | Output/FusionCatcher/STAR build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D)) && mkdir -p $(@D)
	if [[ -e "Data/WTS/$(PAIRED_NORMAL)_*.fastq.gz" ]]; then $(DOCKER) fusioncatcher.py --input=$(addprefix /,$(word 1,$(sort $(filter %.fastq.gz,$^)))),$(addprefix /,$(word 2,$(sort $(filter %.fastq.gz,$^)))) --output=$(addprefix /,$(@D)) --data=$(addprefix /,$(FUSIONCATCHER_DATA_DIR)) --tmp=$(TMPFS) --threads=$(CPUS) --no-update-check --aligners=star --skip-conversion-grch37 --Xmx=$(MEMS)g 1> $(@D).stdout 2> $(@D).stderr; else $(DOCKER) fusioncatcher.py --input=$(addprefix /,$(word 1,$(sort $(filter %.fastq.gz,$^)))),$(addprefix /,$(word 2,$(sort $(filter %.fastq.gz,$^)))) --normal=$(addprefix /,$(word 1,$(sort $(wildcard Data/WTS/$(PAIRED_NORMAL)_*.fastq.gz)))),$(addprefix /,$(word 2,$(sort $(wildcard Data/WTS/$(PAIRED_NORMAL)_*.fastq.gz)))) --output=$(addprefix /,$(@D)) --data=$(addprefix /,$(FUSIONCATCHER_DATA_DIR)) --threads=$(CPUS) --no-update-check --aligners=star --skip-conversion-grch37 --Xmx=$(MEMS)g 1> $(@D).stdout 2> $(@D).stderr; fi;

step26: $(addprefix Output/FusionCatcher/STAR/,$(addsuffix /preliminary-candidate-fusions-found-in-normal.log,$(WTS_TUMOR_ID)))
.PHONY: step26

# Step ## (BisqueRNA)
# https://www.rdocumentation.org/packages/BisqueRNA/versions/1.0.5
# https://cran.r-project.org/web/packages/BisqueRNA/vignettes/bisque.html
Output/BisqueRNA:
	mkdir -p $@

# Step ## (MuSiC)
# https://xuranw.github.io/MuSiC/articles/MuSiC.html
Output/MuSiC:
	mkdir -p $@

# Step ## (SCDC)
# https://meichendong.github.io/SCDC/
# https://meichendong.github.io/SCDC/articles/SCDC.html
Output/SCDC:
	mkdir -p $@

# Step ## (Fusioncatcher)
# Step ## (STAR-fusion)
# Step ## (MesKit)
# Step ## (WTS - samtools depth)
# Step ## (Picard - CollectMultipleMetrics)
# Step ## (PureCN)
# Step ## (Histogram of file size)
