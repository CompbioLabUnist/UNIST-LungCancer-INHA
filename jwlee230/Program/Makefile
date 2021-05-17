# Variables
DATE := $(shell date "+%y%m%d")
RANDOM := $(shell bash -c 'echo $$RANDOM')
IMAGE_NAME = lungcancer:latest
CPUS = 30
MEMS = 100
PARALLEL = 1
PWD := $(shell pwd)
TOOLS = 
DOCKER = docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) $(IMAGE_NAME)
PYTHON = python3 -B
CHROMOSOMES = chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
PAIRED_NORMAL = $(shell echo $* | grep --perl-regexp --only-match "^(cn)?\d+")N
TMPFS = /tmpfs

# Options
VOLUME_OPTS = --volume $(abspath Output):/Output --volume $(abspath Data):/Data:ro --volume $(abspath Python):/Python --volume $(abspath R):/R:ro --volume $(abspath Bash):/Bash:ro --tmpfs $(TMPFS)
RUN_OPTS = --tty --cpus="$(CPUS)" --memory="$(MEMS)G"
JAVA_OPTS = -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Xmx$(MEMS)g
GATK = gatk --java-options "$(JAVA_OPTS)"
PICARD = java $(JAVA_OPTS) -jar /Tools/picard.jar
VCF2MAF = perl /Tools/vcf2maf-1.6.19/vcf2maf.pl
MUTENRICHER = python3 /Tools/MutEnricher-1.3.1/mutEnricher.py

# Files
REFERENCE_DIR = Data/hg38
REFERENCE_FASTA = $(REFERENCE_DIR)/Homo_sapiens_assembly38.fasta
REFERENCE_KNOWN_VCF = $(wildcard $(REFERENCE_DIR)/*.vcf)
REFERENCE_DBSNP_VCF = $(REFERENCE_DIR)/dbsnp_146.hg38.vcf
REFERENCE_1000G_VCF = $(REFERENCE_DIR)/1000G_phase1.snps.high_confidence.hg38.vcf
REFERENCE_SIZE_TSV = $(REFERENCE_DIR)/hg38.chrom.sizes
REFERENCE_CENTROMERE = Data/UCSC_hg38_centromere.txt
LIBRARY_BED = Data/S04380110_Regions.bed
GENCODE_GTF = Data/GENCODE/gencode.v38.annotation.gtf.gz
TREMBL_TSV = Data/Homo_sapiens.GRCh38.104.uniprot.tsv.gz

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
ADC_PATIENTS = 39 51 52 54 55 61 62 65 cn80 cn90 cn93 cn95 cn96 cn97 cn99 cn101 cn102 cn103 cn105 cn106 cn107 cn108 cn109 cn110 cn111 cn112 cn113 cn115 cn116
SQC_WES_SAMPLES = 8P1 8C1 12P1 12C1 14P1 14C1 15P1 15C1 33D1 34P1 34C1 38P1 38C1 48P1 48A1 48A2 48A3 49P1 49C1 53P1 53A1 57P1 57C1 57D1 57D1 58P1 58A1 58A2 58A3 58A4 63P1 63C1 63C2 64P1 64D1 cn1P cn1C cn1D cn2P cn2C cn3P cn3C cn4P cn4C cn5P cn5C cn7P cn7C cn10P cn10C cn11P cn11C cn13P cn13C cn15P cn15C cn16P cn16C cn17P cn17C cn18P cn18C cn21P cn21C cn23P cn23C cn26P cn26C cn27P cn27C cn30P cn30C cn31P cn31C cn32P cn32C cn33P cn33C cn34P cn34C cn35P cn35C cn36P cn36C cn37P cn37C cn38P cn38C cn39P cn39C cn40P cn40C cn41P cn41C cn42P cn42C cn43P cn43C cn44P cn44C cn45P cn45C cn47P cn47C cn48P cn48C cn49P cn49C cn50P cn50C cn51P cn51C cn52P cn52C cn53P cn53C cn54P cn54C cn55P cn55C cn56P cn56C cn57P cn57C cn60P cn60C cn61P cn61C cn62P cn62C cn63P cn63C cn64P cn64C cn65P cn65C cn66P cn66C cn67P cn67C cn68P cn68C cn69P cn69C cn70P cn70C cn71P cn71C cn72P cn72C cn73P cn73C cn74P cn74C cn75P cn75C cn76P cn76C cn78P cn78C cn79P cn79C
SQC_WTS_SAMPLES = cn1P cn1C cn1D cn2P cn2C cn3P cn3C cn4P cn4C cn5P cn5C cn6P cn6C cn7P cn7C cn8P cn8C cn10P cn10C cn16P cn16C cn17P cn17C cn23P cn23C cn30P cn30C cn32P cn32C cn38P cn38C cn39P cn39C cn40P cn40C cn44P cn44C cn49P cn49C cn52P cn52C cn53P cn53C cn54P cn54C cn55P cn55C cn57P cn57C cn60P cn60C cn63P cn63C cn64P cn64C cn65P cn65C cn66P cn66C cn67P cn67C cn68P cn68C cn69P cn69C cn70P cn70C cn74P cn74C cn75P cn75C cn76P cn76C cn77P cn77C cn77D cn78P cn78C cn79P cn79C
ADC_WES_SAMPLES = 39P1 39A1 39A2 51P1 51A1 52P1 52M1 54P1 54A1 55P1 55A1 55A2 61P1 61A1 61A2 62P1 62A1 62A2 62A3 62A4 65P1 65A1 cn93P cn93C cn95P cn95A cn97P cn97C cn101P cn101C cn105P cn105C cn106P cn106C cn107P cn107C cn107A1 cn108P cn108C cn110P cn110C
ADC_WTS_SAMPLES = cn90P cn90C cn95P cn95A cn96P cn96C cn99P cn99C cn101P cn101C cn102P cn102C cn103C cn107P cn107C cn107A1 cn109P cn110P cn110C cn111P cn111C cn112P cn112C cn113P cn113C cn115P cn115C cn115A cn116P cn116C

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
	wget "https://github.com/alexdobin/STAR/archive/2.7.7a.tar.gz" -O $@
TOOLS += Docker/Tools/STAR.tar.gz

Docker/Tools/Bowtie2.zip: | Docker/Tools
	wget "https://github.com/BenLangmead/bowtie2/releases/download/v2.4.2/bowtie2-2.4.2-linux-x86_64.zip" -O $@
TOOLS += Docker/Tools/Bowtie2.zip

Docker/Tools/vcf2maf.tar.gz: | Docker/Tools
	wget "https://github.com/mskcc/vcf2maf/archive/v1.6.19.tar.gz" -O $@
TOOLS += Docker/Tools/vcf2maf.tar.gz

Docker/Tools/VEP.zip: | Docker/Tools
	wget "https://github.com/Ensembl/ensembl-vep/archive/release/103.zip" -O $@
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
latest: step18 # step15
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
	$(DOCKER) bowtie2 --threads $(CPUS) --rg-id "$(@F)" --time --qc-filter -x $(basename $(notdir $(addprefix /,$(filter %.fasta,$^)))) -1 $(addprefix /,$(word 1,$(sort $(filter %.fastq.gz,$^)))) -2 $(addprefix /,$(word 2,$(sort $(filter %.fastq.gz,$^)))) -S $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

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

step05: $(addprefix Output/WES/,$(addsuffix .BWA.Sort.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .BWA.Sort.bai,$(WES_ID)))
.PHONY: step05

# Step 06 (WES: Sort.bam >> Sort.MarkDuplicate.bam)
Output/WES/%.Sort.dm-tmp.bam: $(REFERENCE_FASTA) Output/WES/%.Sort.bam Output/WES/%.Sort.bai | build.log
	$(DOCKER) $(PICARD) MarkDuplicates --INPUT $(addprefix /,$(filter %.bam,$^)) --REFERENCE_SEQUENCE $(addprefix /,$(filter %.fasta,$^)) --OUTPUT $(addprefix /,$@) --METRICS_FILE $(addprefix /,$@.metrics) --ASSUME_SORTED true --VALIDATION_STRINGENCY SILENT --COMPRESSION_LEVEL 9 1> $@.stdout 2> $@.stderr

Output/WES/%.Sort.MarkDuplicate.bam: Bash/remove_CR.bash Output/WES/%.Sort.dm-tmp.bam | build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/WES/%.Sort.MarkDuplicate.bai: Output/WES/%.Sort.MarkDuplicate.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step06: $(addprefix Output/WES/,$(addsuffix .BWA.Sort.dm-tmp.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.bai,$(WES_ID)))
.PHONY: step06

# Step 07 (WES: Sort.MarkDuplicate.bam >> Sort.MarkDuplicate.BQSR.bam)
Output/WES/%.Sort.MarkDuplicate.BQSR.table: $(REFERENCE_FASTA) $(REFERENCE_KNOWN_VCF) Output/WES/%.Sort.MarkDuplicate.bam Output/WES/%.Sort.MarkDuplicate.bai | build.log
	$(DOCKER) $(GATK) BaseRecalibrator --input $(addprefix /,$(filter %.bam,$^)) $(foreach site,$(addprefix /,$(filter %.vcf,$^)),--known-sites $(site)) --output $(addprefix /,$@) --reference $(addprefix /,$(filter %.fasta,$^)) 1> $@.stdout 2> $@.stderr

Output/WES/%.Sort.MarkDuplicate.BQSR.bam: Output/WES/%.Sort.MarkDuplicate.bam Output/WES/%.Sort.MarkDuplicate.BQSR.table $(REFERENCE_FASTA) | build.log
	$(DOCKER) $(GATK) ApplyBQSR --bqsr-recal-file $(addprefix /,$(filter %.table,$^)) --input $(addprefix /,$(filter %.bam,$^)) --output $(addprefix /,$@) --reference $(addprefix /,$(filter %.fasta,$^)) 1> $@.stdout 2> $@.stderr

Output/WES/%.Sort.MarkDuplicate.BQSR.bai: Output/WES/%.Sort.MarkDuplicate.BQSR.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step07: $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.bai,$(WES_ID)))
.PHONY: step07

clean-step07:
	rm -fv $(filter-out %.Sort.MarkDuplicate.BQSR.bam,$(shell find Output/WES -size +1G))
.PHONY: clean-step07

# Step 08 (Panel of Normal: Run Mutect2 in tumor-only mode for each normal sample)
Output/PanelOfNormal/BWA:
	mkdir -p $@

Output/PanelOfNormal/BWA/%.vcf: $(REFERENCE_FASTA) Output/WES/%.BWA.Sort.MarkDuplicate.BQSR.bam Output/WES/%.BWA.Sort.MarkDuplicate.BQSR.bai | Output/PanelOfNormal/BWA build.log
	$(DOCKER) $(GATK) Mutect2 --reference $(addprefix /,$(filter %.fasta,$^)) --input $(addprefix /,$(filter %.bam,$^)) --output $(addprefix /,$@) --max-mnp-distance 0 1> $@.stdout 2> $@.stderr

step08: $(addprefix Output/PanelOfNormal/BWA/,$(addsuffix .vcf,$(WES_NORMAL_ID)))
.PHONY: step08

# Step 09 (Panel of Normal: Aggregate VCF files)
Output/PanelOfNormal/BWA.genomicsdb/callset.json: $(REFERENCE_FASTA) $(addprefix Output/PanelOfNormal/BWA/,$(addsuffix .vcf,$(WES_NORMAL_ID))) | build.log
	$(DOCKER) $(GATK) GenomicsDBImport --reference $(addprefix /,$(filter %.fasta,$^)) --genomicsdb-workspace-path $(addprefix /,$(@D)) $(foreach sample,$(addprefix /,$(filter %.vcf,$^)),--variant $(sample)) --reader-threads $(CPUS) --max-num-intervals-to-import-in-parallel $(CPUS) --overwrite-existing-genomicsdb-workspace true --merge-input-intervals true $(foreach chr,$(CHROMOSOMES),--intervals $(chr)) 1> $(@D).stdout 2> $(@D).stderr
	$(DOCKER) chmod --recursive 777 $(addprefix /,$(@D))

Output/PanelOfNormal/PON.BWA.vcf.gz: $(REFERENCE_FASTA) Output/PanelOfNormal/BWA.genomicsdb/callset.json | build.log
	$(DOCKER) $(GATK) CreateSomaticPanelOfNormals --reference $(addprefix /,$(filter %.fasta,$^)) --variant "gendb://$(addprefix /,$(dir $(filter %.json,$^)))" --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step09: Output/PanelOfNormal/BWA.genomicsdb/callset.json Output/PanelOfNormal/PON.BWA.vcf.gz
.PHONY: step09

# Step 10 (WES - samtools depth)
Output/WESDepth:
	mkdir -p $@

Output/WESDepth/%.tsv: $(REFERENCE_FASTA) $(LIBRARY_BED) Output/WES/%.Sort.MarkDuplicate.BQSR.bam Output/WES/%.Sort.MarkDuplicate.BQSR.bai | Output/WESDepth build.log
	$(DOCKER) samtools depth -H -b $(addprefix /,$(filter %.bed,$^)) --reference $(addprefix /,$(filter %.fasta,$^)) -o $(addprefix /,$@) $(addprefix /,$(filter %.bam,$^)) 1> $@.stdout 2> $@.stderr

Output/WESDepth/BWA.pdf: Python/analysis_depth.py $(addprefix Output/WESDepth/,$(addsuffix .BWA.tsv,$(WES_ID))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step10: $(addprefix Output/WESDepth/,$(addsuffix .BWA.tsv,$(WES_ID))) Output/WESDepth/BWA.pdf
.PHONY: step10

# Step 11 (Sequenza: Make reference)
Output/Sequenza:
	mkdir -p $@

Output/Sequenza/hg38.wig.gz: $(REFERENCE_FASTA) | Output/Sequenza build.log
	$(DOCKER) sequenza-utils gc_wiggle --fasta $(addprefix /,$^) -o $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Sequenza/hg37.wig.gz: Data/human_g1k_v37_decoy.fasta | Output/Sequenza build.log
	$(DOCKER) sequenza-utils gc_wiggle --fasta $(addprefix /,$^) -o $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step11: Output/Sequenza/hg38.wig.gz Output/Sequenza/hg37.wig.gz
.PHONY: step11

# Step 12 (Sequenza: Normal - Tumor pairs)
Output/Sequenza/BWA:
	mkdir -p $@

# I know that this recipe contains unnecessary dependencies; but I cannot see the solution.
Output/Sequenza/BWA/%.seqz.gz: $(REFERENCE_FASTA) Output/Sequenza/hg38.wig.gz Output/WES/%.BWA.Sort.MarkDuplicate.BQSR.bam $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.bam,$(WES_NORMAL_ID))) | Output/Sequenza/BWA build.log
	@echo "This step requires multiple CPUs!!"
	$(DOCKER) sequenza-utils bam2seqz --normal "/Output/WES/$(PAIRED_NORMAL).BWA.Sort.MarkDuplicate.BQSR.bam" --tumor "/Output/WES/$*.BWA.Sort.MarkDuplicate.BQSR.bam" --fasta $(addprefix /,$(filter %.fasta,$^)) -gc $(addprefix /,$(filter %.wig.gz,$^)) --output $(addprefix /,$@) --chromosome $(CHROMOSOMES) --parallel $(CPUS) 1> $@.stdout 2> $@.stderr
	cd $(@D) && gzip --decompress --stdout --keep $(foreach chr,$(CHROMOSOMES),$*_$(chr).seqz.gz) | gawk '{if (NR!=1 && $$1 != "chromosome") {print $$0}}' | bgzip --threads $(CPUS) > $(@F)
	$(DOCKER) rm $(addprefix /$(@D)/,$(foreach chr,$(CHROMOSOMES),$*_$(chr).seqz.gz))

Output/Sequenza/BWA/%.seqz.gz.tbi: Output/Sequenza/BWA/%.seqz.gz | build.log
	$(DOCKER) tabix --force --sequence 1 --begin 2 --end 2 --skip-lines 1 $(addprefix /,$<) 1> $@.stdout 2> $@.stderr

step12: $(addprefix Output/Sequenza/BWA/,$(addsuffix .seqz.gz,$(WES_TUMOR_ID))) $(addprefix Output/Sequenza/BWA/,$(addsuffix .seqz.gz.tbi,$(WES_TUMOR_ID)))
.PHONY: step12

# Step 13 (Sequenza: Draw the plots)
Output/Sequenza/BWA/%.tar.gz: R/run_sequenza.R Output/Sequenza/BWA/%.seqz.gz Output/Sequenza/BWA/%.seqz.gz.tbi | build.log
	$(DOCKER) Rscript --vanilla $(addprefix /,$(filter %.R,$^)) --sampleid $* --input $(addprefix /,$(filter %.seqz.gz,$^)) --outdir $(addprefix /,$(@D)) --parallel $(CPUS) 1> $@.stdout 2> $@.stderr
	cd $(abspath $(@D)) && tar -zcvf $(@F) $*_* && rm -fv $*_*

Output/Sequenza/BWA-sequenza-SQC.pdf: Python/aggregate_sequenza.py $(addprefix Output/Sequenza/BWA/,$(addsuffix .tar.gz,$(SQC_WES_SAMPLES))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-sequenza-ADC.pdf: Python/aggregate_sequenza.py $(addprefix Output/Sequenza/BWA/,$(addsuffix .tar.gz,$(ADC_WES_SAMPLES))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-SQC.pdf: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/BWA/,$(addsuffix .tar.gz,$(SQC_WES_SAMPLES))) Data/hg38/hg38.chrom.sizes | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-ADC.pdf: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/BWA/,$(addsuffix .tar.gz,$(ADC_WES_SAMPLES))) Data/hg38/hg38.chrom.sizes | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step13: $(addprefix Output/Sequenza/BWA/,$(addsuffix .tar.gz,$(WES_TUMOR_ID))) Output/Sequenza/BWA-sequenza-SQC.pdf Output/Sequenza/BWA-sequenza-ADC.pdf Output/Sequenza/BWA-genome-SQC.pdf Output/Sequenza/BWA-genome-ADC.pdf
.PHONY: step13

# Step 14 (Mutect2)
Output/Mutect2/BWA:
	mkdir -p $@

# I know that this recipe contains unnecessary dependencies; but I cannot see the solution.
Output/Mutect2/BWA/%.BWA.vcf: $(REFERENCE_FASTA) Output/PanelOfNormal/PON.BWA.vcf.gz Output/WES/%.BWA.Sort.MarkDuplicate.BQSR.bam $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.bam,$(WES_NORMAL_ID))) | Output/Mutect2/BWA build.log
	$(DOCKER) $(GATK) Mutect2 --reference $(addprefix /,$(filter %.fasta,$^)) --input "/Output/WES/$*.BWA.Sort.MarkDuplicate.BQSR.bam" --input "/Output/WES/$(PAIRED_NORMAL).BWA.Sort.MarkDuplicate.BQSR.bam" --output $(addprefix /,$@) --normal-sample "$(PAIRED_NORMAL).BWA.sam" --panel-of-normals $(addprefix /,$(filter %.vcf.gz,$^)) 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA/%.BWA.filter.vcf: $(REFERENCE_FASTA) Output/Mutect2/BWA/%.BWA.vcf | build.log
	$(DOCKER) $(GATK) FilterMutectCalls --reference $(addprefix /,$(filter %.fasta,$^)) --variant $(addprefix /,$(filter %.vcf,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA/%.BWA.PASS.vcf: Bash/select_pass.bash Output/Mutect2/BWA/%.BWA.filter.vcf | build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA/%.BWA.PASS.maf: $(REFERENCE_FASTA) Output/Mutect2/BWA/%.BWA.PASS.vcf | build.log
	$(DOCKER) $(VCF2MAF) --vep-path "/Tools/ensembl-vep-release-103" --vep-data "/Data" --vep-forks $(CPUS) --ncbi-build "GRCh38" --input-vcf $(addprefix /,$(filter %.vcf,$^)) --output $(addprefix /,$@) --tumor-id "$*.BWA.sam" --normal-id "$(PAIRED_NORMAL).BWA.sam" --ref-fasta $(addprefix /,$(filter %.fasta,$^)) --tmp-dir $(TMPFS) 1> $@.stdout 2> $@.stderr

Output/Mutect2/Mutect2-SQC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(SQC_WES_SAMPLES))) Output/merged.sqc.BWA.PASS.maf.gene.mutationburden.txt | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Mutect2/Mutect2-ADC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(ADC_WES_SAMPLES))) Output/merged.adc.BWA.PASS.maf.gene.mutationburden.txt | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step14: $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.vcf,$(WES_TUMOR_ID))) $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.filter.vcf,$(WES_TUMOR_ID))) $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.vcf,$(WES_TUMOR_ID))) $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) Output/Mutect2/Mutect2-SQC.pdf Output/Mutect2/Mutect2-ADC.pdf
.PHONY: step14

# Step 15 (HaplotypeCaller)
Output/HaplotypeCaller/BWA:
	mkdir -p $@

Output/HaplotypeCaller/BWA/%.BWA.g.vcf.gz: $(REFERENCE_FASTA) Output/WES/%.BWA.Sort.MarkDuplicate.BQSR.bam $(REFERENCE_DBSNP_VCF) | Output/HaplotypeCaller/BWA build.log
	$(DOCKER) $(GATK) HaplotypeCaller --reference $(addprefix /,$(filter %.fasta,$^)) --input $(addprefix /,$(filter %.bam,$^)) --output $(addprefix /,$@) --emit-ref-confidence GVCF --dbsnp $(addprefix /,$(filter %.vcf,$^)) --native-pair-hmm-threads $(CPUS) --native-pair-hmm-use-double-precision true 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/HaplotypeCaller.genomicsdb/callset.json: $(REFERENCE_FASTA) $(addprefix Output/HaplotypeCaller/BWA/,$(addsuffix .BWA.g.vcf.gz,$(WES_TUMOR_ID))) | build.log
	$(DOCKER) $(GATK) GenomicsDBImport --reference $(addprefix /,$(filter %.fasta,$^)) --genomicsdb-workspace-path $(addprefix /,$(@D))  $(foreach sample,$(addprefix /,$(filter %.vcf,$^)),--variant $(sample)) --reader-threads $(CPUS) --max-num-intervals-to-import-in-parallel $(CPUS) --overwrite-existing-genomicsdb-workspace true --merge-input-intervals true $(foreach chr,$(CHROMOSOMES),--intervals $(chr)) 1> $@.stdout 2> $@.stderr
	$(DOCKER) chmod --recursive 777 $(addprefix /,$(@D))

Output/HaplotypeCaller/HaplotypeCaller.vcf.gz: $(REFERENCE_FASTA) Output/HaplotypeCaller/HaplotypeCaller.genomicsdb/callset.json $(REFERENCE_DBSNP_VCF) | build.log
	$(DOCKER) $(GATK) GenotypeGVCFs --reference $(addprefix /,$(filter %.fasta,$^)) --variant "gendb://$(addprefix /,$(dir $(filter %.json,$^)))" --output $(addprefix /,$@) --dbsnp $(addprefix /,$(filter %.vcf,$^)) 1> $@.stdout 2> $@.stderr

step15: $(addprefix Output/HaplotypeCaller/BWA/,$(addsuffix .BWA.g.vcf.gz,$(WES_TUMOR_ID))) Output/HaplotypeCaller/HaplotypeCaller.genomicsdb/callset.json Output/HaplotypeCaller/HaplotypeCaller.vcf.gz
.PHONY: step15
# TODO: wallace >> /BiO/Research/UNIST-Pancreas_Liver_Metastasis-SMC-2018-05/Workspace/seunghoon/hg38/haplotypecaller

# Step 16 (MutEnricher)
Output/MutEnricher Output/MutEnricher/BWA:
	mkdir -p $@

Output/MutEnricher/BWA/%.vcf.gz: Output/Mutect2/BWA/%.BWA.PASS.vcf | build.log Output/MutEnricher/BWA
	$(DOCKER) bcftools view --output-type z --output $(addprefix /,$@) --threads $(CPUS) $(addprefix /,$^) 1> $@.stdout 2> $@.stderr
	$(DOCKER) htsfile $(addprefix /,$@)

Output/MutEnricher/BWA/%.vcf.gz.tbi: Output/MutEnricher/BWA/%.vcf.gz | build.log
	$(DOCKER) tabix --force $(addprefix /,$^) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/SQC.tsv: Python/make_list.py $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz,$(SQC_WES_SAMPLES))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz.tbi,$(SQC_WES_SAMPLES))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$< $(filter %.vcf.gz,$^) $@) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/ADC.tsv: Python/make_list.py $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz,$(ADC_WES_SAMPLES))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz.tbi,$(ADC_WES_SAMPLES))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$< $(filter %.vcf.gz,$^) $@) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/hg38.gtf.gz: | Output/MutEnricher
	wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz" -O $@
	touch $@

Output/MutEnricher/gene_covariates.txt: Output/MutEnricher/hg38.gtf.gz $(REFERENCE_FASTA) | build.log
	$(DOCKER) $(PYTHON) /Tools/MutEnricher-1.3.1/utilities/get_gene_covariates.py --outname $(addprefix /,$@) --processors $(CPUS) $(addprefix /,$^) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/%_gene_enrichments.txt Output/MutEnricher/%_hotspot.txt Output/MutEnricher/%_gene_hotspot_Fisher_enrichments.txt Output/MutEnricher/%_gene_data.pkl Output/MutEnricher/%.log: Output/MutEnricher/gene_covariates.txt Output/MutEnricher/hg38.gtf.gz Output/MutEnricher/%.tsv | build.log
	$(DOCKER) $(MUTENRICHER) coding --outdir $(addprefix /,$(@D)) --prefix $* --processors $(CPUS) --use-local --covariates-file $(addprefix /,$(filter %.txt,$^)) $(addprefix /,$(filter %.gtf.gz,$^)) $(addprefix /,$(filter %.tsv,$^)) 1> $@.stdout 2> $@.stderr

step16: Output/MutEnricher/ADC.log Output/MutEnricher/SQC.log
.PHONY: step16

# Step 17 (RSEM: Prepare reference genome)
Output/RSEM:
	mkdir -p $@

Output/RSEM/hg38.gtf.gz: | Output/RSEM
	wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz" -O $@
	touch $@

Output/RSEM/hg38.gtf: Output/RSEM/hg38.gtf.gz | build.log
	$(DOCKER) gzip --decompress --keep --force $(addprefix /,$<)

Output/RSEM/hg38.seq Output/RSEM/hg38.grp: $(REFERENCE_FASTA) Output/RSEM/hg38.gtf | build.log
	$(DOCKER) rsem-prepare-reference --gtf $(addprefix /,$(filter %.gtf,$^)) --bowtie2 --star --num-threads $(CPUS) $(addprefix /,$(filter %.fasta,$^)) /$(@D)/"hg38" 1> $@.stdout 2> $@.stderr

step17: Output/RSEM/hg38.gtf Output/RSEM/hg38.seq
.PHONY: step17

# Step 18 (RSEM: Calculate expression)
Output/RSEM/Bowtie2 Output/RSEM/STAR:
	mkdir -p $@

Output/RSEM/Bowtie2/%.genome.bam Output/RSEM/Bowtie2/%.genes.results: Data/WTS/%_*.fastq.gz Output/RSEM/hg38.grp | Output/RSEM/Bowtie2 build.log
	$(DOCKER) rsem-calculate-expression --num-threads $(CPUS) --estimate-rspd --ci-memory $(shell echo "$(MEMS) * 1024" | bc) --bowtie2 --output-genome-bam --paired-end $(addprefix /,$(sort $(filter %.fastq.gz,$^)) $(basename $(filter %.grp,$^)) $(@D)/$*) 1> $@.stdout 2> $@.stderr

Output/RSEM/STAR/%.genome.bam Output/RSEM/STAR/%.genes.results: Data/WTS/%_*.fastq.gz Output/RSEM/hg38.grp | Output/RSEM/STAR build.log
	$(DOCKER) rsem-calculate-expression --num-threads $(CPUS) --estimate-rspd --ci-memory $(shell echo "$(MEMS) * 1024" | bc) --star --star-gzipped-read-file --output-genome-bam --paired-end $(addprefix /,$(sort $(filter %.fastq.gz,$^)) $(basename $(filter %.grp,$^)) $(@D)/$*) 1> $@.stdout 2> $@.stderr

Output/RSEM/SQC.Bowtie2.tsv: Python/aggregate_RSEM.py Output/RSEM/Bowtie2/cn102N.genes.results Output/RSEM/Bowtie2/cn102C.genes.results Output/RSEM/Bowtie2/cn102P.genes.results $(GENCODE_GTF) $(TREMBL_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS)

step18: Output/RSEM/SQC.Bowtie2.tsv # $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(WTS_ID))) # $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(WTS_ID)))
.PHONY: step18

# Step ## (WTS - fq.gz >> bam) - STAR
# Step ## (WTS - samtools depth)
# Step ## (Picard - CollectMultipleMetrics)
# step ## (PureCN)
# Step ## (Histogram of file size)
