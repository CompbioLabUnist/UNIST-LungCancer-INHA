# Variables
DATE := $(shell date "+%y%m%d")
RANDOM := $(shell bash -c 'echo $$RANDOM')
IMAGE_NAME = lungcancer:latest
CPUS = 10
MEMS = 50G
PWD := $(shell pwd)
TOOLS = 
DOCKER = docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) $(IMAGE_NAME)
CHROMOSOMES = chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM

# Options
VOLUME_OPTS = --volume $(abspath Output):/Output --volume $(abspath Data):/Data:ro --volume $(abspath Python):/Python
RUN_OPTS = --tty --cpus="$(CPUS)" --memory="$(MEMS)"
JAVA_OPTS = -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Xmx$(MEMS)
GATK = gatk --java-options "$(JAVA_OPTS)"
PICARD = java $(JAVA_OPTS) -jar /Tools/picard.jar

# Files
REFERENCE_DIR = Data/hg38
REFERENCE_FASTA = $(REFERENCE_DIR)/Homo_sapiens_assembly38.fasta
REFERENCE_KNOWN_VCF = $(wildcard $(REFERENCE_DIR)/*.vcf)
WES_ID = $(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | cut --delimiter _ --fields 1))
WES_NORMAL_ID = $(sort $(shell ls -1 Data/RawData/*.fastq.gz | grep "N" | cut --delimiter / --fields 3 | cut --delimiter _ --fields 1))
WTS_ID = $(sort $(shell ls -1 Data/WTS/*.fastq.gz | cut --delimiter / --fields 3 | cut --delimiter _ --fields 1))
VCF_ID = $(sort $(shell ls -1 Data/VCF/*.vcf | cut --delimiter / --fields 3 | cut --delimiter . --fields 1))
BAM_ID = $(sort $(shell ls -1 Data/BAM/*.bam | cut --delimiter / --fields 3 | cut --delimiter . --fields 1))

# General
all:
.PHONY += all

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

# Docker
build.log: Docker/Dockerfile $(TOOLS) | log Output
	rm -f $@
	docker images | grep $(IMAGE_NAME) && docker rmi $(IMAGE_NAME) || true
	docker build --rm --tag $(IMAGE_NAME) $(dir $(word 1,$^)) | tee $@

build: build.log
.PHONY += build

interactive: build.log
	docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) --interactive $(IMAGE_NAME) /bin/bash || true
.PHONY += interactive

delete: build.log
	docker rmi $(IMAGE_NAME)
	rm -rfv build.log
.PHONY += delete

stop:
	docker rm $(CONTAINER_NAME)
.PHONY += stop

# SGE
run: | log Output build.log
	echo "make -j 15 -C $(PWD) latest" > tmp.sh
	qsub -cwd -l h_vmem=750G -m abe -M "230@fumire.moe" -N Lung_$(DATE) -pe smp 150 -o $(abspath log) -e $(abspath log) tmp.sh
.PHONY += run

# Actual
latest: step06
.PHONY += latest

# Step 01 (FastQC for WES)
Output/WES/FastQC:
	mkdir -p $@

Output/WES/FastQC/%_fastqc.zip: Data/RawData/%.fastq.gz | Output/WES/FastQC build.log
	$(DOCKER) fastqc --outdir $(addprefix /,$(@D)) --format fastq --noextract $(addprefix /,$^) 1> $@.stdout 2> $@.stderr

Output/WES/FastQC/FastQC_WES.png: Python/aggregate_fastqc.py $(addprefix Output/WES/FastQC/,$(addsuffix _fastqc.zip,$(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | cut --delimiter . --fields 1))))
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step01: Output/WES/FastQC/FastQC_WES.png
.PHONY += step01

# Step 02 (FastQC for WTS)
Output/WTS/FastQC:
	mkdir -p $@

Output/WTS/FastQC/%_fastqc.zip: Data/WTS/%.fastq.gz | Output/WTS/FastQC build.log
	$(DOCKER) fastqc --outdir $(addprefix /,$(@D)) --format fastq --noextract $(addprefix /,$^) 1> $@.stdout 2> $@.stderr

Output/WTS/FastQC/FastQC_WTS.png: Python/aggregate_fastqc.py $(addprefix Output/WTS/FastQC/,$(addsuffix _fastqc.zip,$(sort $(shell ls -1 Data/WTS/*.fastq.gz | cut --delimiter / --fields 3 | cut --delimiter . --fields 1))))
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step02: Output/WTS/FastQC/FastQC_WTS.png
.PHONY += step02

# Step 03 (WES: fq.gz >> sam)
Output/WES/BAM:
	mkdir -p $@

Output/WES/BAM/%.BWA.sam: $(REFERENCE_FASTA) Data/RawData/%_*.fastq.gz | Output/WES/BAM build.log
	$(DOCKER) bwa mem -M -t $(CPUS) -R "@RG\tID:$(@F)\tPL:ILLUMINA\tLB:$(@F)\tSM:$(@F)\tCN:UNIST" -v 3 $(addprefix /,$(filter %.fasta,$^)) $(sort $(addprefix /,$(filter %.fastq.gz,$^))) 1> $@ 2> $@.stderr
	$(DOCKER) sed --in-place "/^\[.*/d" $(addprefix /,$@)

Output/WES/BAM/%.Bowtie2.sam: $(REFERENCE_FASTA) Data/RawData/%_*.fastq.gz | Output/WES/BAM build.log
	$(DOCKER) bowtie2 -x $(basename $(addprefix /,$(filter %.fasta,$^))) -1 $(addprefix /,$(word 1,$(sort $(filter %.fastq.gz,$^)))) -2 $(addprefix /,$(word 2,$(sort $(filter %.fastq.gz,$^)))) -S $(addprefix /,$@) --threads $(CPUS) --rg-id $(@F) --mm 1> $@.stdout 2> $@.stderr

step03: $(sort $(addprefix Output/WES/BAM/,$(addsuffix .BWA.sam,$(WES_ID))))
.PHONY += step03

# Step 04 (WES: sam >> bam)
Output/WES/BAM/%.bam: $(REFERENCE_FASTA) Output/WES/BAM/%.sam | build.log
	$(DOCKER) samtools view -b -h --threads $(CPUS) --reference $(addprefix /,$(filter %.fasta,$^)) $(addprefix /,$(filter %.sam,$^)) -o $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step04: $(sort $(addprefix Output/WES/BAM/,$(addsuffix .BWA.bam,$(WES_ID))))
.PHONY += step04

# Step 05 (WES: bam >> st.bam)
Output/WES/BAM/%.st.bam: $(REFERENCE_FASTA) Output/WES/BAM/%.bam | build.log
	$(DOCKER) samtools sort -l 9 -m 5G --threads $(CPUS) --reference $(addprefix /,$(filter %.fasta,$^)) -o $(addprefix /,$@) $(addprefix /,$(filter %.bam,$^)) 1> $@.stdout 2> $@.stderr

Output/WES/BAM/%.st.bai: Output/WES/BAM/%.st.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step05: $(sort $(addprefix Output/WES/BAM/,$(addsuffix .BWA.st.bam,$(WES_ID))) $(addprefix Output/WES/BAM/,$(addsuffix .BWA.st.bai,$(WES_ID))))
.PHONY += step05

# Step 06 (WES: st.bam >> st.dm.bam)
Output/WES/BAM/%.st.dm-tmp.bam: $(REFERENCE_FASTA) Output/WES/BAM/%.st.bam Output/WES/BAM/%.st.bai | build.log
	$(DOCKER) $(PICARD) MarkDuplicates --INPUT $(addprefix /,$(filter %.bam,$^)) --REFERENCE_SEQUENCE $(addprefix /,$(filter %.fasta,$^)) --OUTPUT $(addprefix /,$@) --METRICS_FILE $(addprefix /,$@.metrics) --ASSUME_SORTED true --VALIDATION_STRINGENCY SILENT --COMPRESSION_LEVEL 9 1> $@.stdout 2> $@.stderr

Output/WES/BAM/%.st.dm-tmp.sam: Output/WES/BAM/%.st.dm-tmp.bam | build.log
	$(DOCKER) samtools view -h -o $(addprefix /,$@) $(addprefix /,$<) 1> $@.stdout 2> $@.stderr

Output/WES/BAM/%.st.dm.sam: Output/WES/BAM/%.st.dm-tmp.sam | build.log
	$(DOCKER) sed "s/\r//" $(addprefix /,$<) > $@

Output/WES/BAM/%.st.dm.bam: Output/WES/BAM/%.st.dm.sam | build.log
	$(DOCKER) samtools view -b -o $(addprefix /,$@) $(addprefix /,$<)

Output/WES/BAM/%.st.dm.bai: Output/WES/BAM/%.st.dm.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step06: $(sort $(addprefix Output/WES/BAM/,$(addsuffix .BWA.st.dm.bam,$(WES_ID))) $(addprefix Output/WES/BAM/,$(addsuffix .BWA.st.dm.bai,$(WES_ID))))
.PHONY += step06

# Step 07 (WES: st.dm.bam >> st.dm.rc.bam)
Output/WES/BAM/%.st.dm.rc.table: $(REFERENCE_FASTA) $(REFERENCE_KNOWN_VCF) Output/WES/BAM/%.st.dm.bam Output/WES/BAM/%.st.dm.bai | build.log
	$(DOCKER) $(GATK) BaseRecalibrator --input $(addprefix /,$(filter %.bam,$^)) $(foreach site,$(addprefix /,$(filter %.vcf,$^)),--known-sites $(site)) --output $(addprefix /,$@) --reference $(addprefix /,$(filter %.fasta,$^)) 1> $@.stdout 2> $@.stderr

Output/WES/BAM/%.st.dm.rc.bam: Output/WES/BAM/%.st.dm.bam Output/WES/BAM/%.st.dm.rc.table $(REFERENCE_FASTA) | build.log
	$(DOCKER) $(GATK) ApplyBQSR --bqsr-recal-file $(addprefix /,$(filter %.table,$^)) --input $(addprefix /,$(filter %.bam,$^)) --output $(addprefix /,$@) --reference $(addprefix /,$(filter %.fasta,$^)) 1> $@.stdout 2> $@.stderr

Output/WES/BAM/%.st.dm.rc.bai: Output/WES/BAM/%.st.dm.rc.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step07: $(sort $(addprefix Output/WES/BAM/,$(addsuffix .BWA.st.dm.rc.bam,$(WES_ID))) $(addprefix Output/WES/BAM/,$(addsuffix .BWA.st.dm.rc.bai,$(WES_ID))))
.PHONY += step07

clean-step07:
	rm -fv $(filter-out %.st.dm.rc.bam,$(shell find /Output/WES/BAM -size +1G))
.PHONY += clean-step07

# Step 08 (Panel of Normal: Run Mutect2 in tumor-only mode for each normal sample)
Output/PanelOfNormal:
	mkdir -p $@

Output/PanelOfNormal/%.vcf: $(REFERENCE_FASTA) Output/WES/BAM/%.st.dm.rc.bam Output/WES/BAM/%.st.dm.rc.bai | Output/PanelOfNormal build.log
	$(DOCKER) gatk Mutect2 --reference $(addprefix /,$(filter %.fasta,$^)) --input $(addprefix /,$(filter %.bam,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step08: $(sort $(addprefix Output/WES/PanelOfNormal/,$(addsuffix .vcf,$(BAM_NORMAL_ID))))
.PHONY += step08

# Step 09 (Panel of Normal: Aggregate VCF files)
Output/PanelOfNormal/PON.vcf: $(REFERENCE_FASTA) $(sort $(addprefix Output/PanelOfNormal/,$(addsuffix .vcf,$(BAM_NORMAL_ID)))) | build.log
	$(DOCKER) gatk CreateSomaticPanelOfNormals --reference $(addprefix /,$(filter %.fasta,$^)) --output $(addprefix /,$@) $(foreach sample,$(addprefix /,$(filter %.vcf,$^)),--variant $(sample)) 1> $@.stdout 2> $@.stderr

step09: Output/PanelOfNormal/PON.vcf
.PHONY += step09

# Step 10 (WES - samtools depth)
Output/WES/Depth:
	mkdir -p $@

Output/WES/Depth/%.tsv: $(REFERENCE_FASTA) Output/WES/BAM/%.st.dm.rc.bam Output/WES/BAM/%.st.dm.rc.bai | Output/WES/Depth build.log
	$(DOCKER) samtools depth -H --reference $(addprefix /,$(filter %.fasta,$^)) -o $(addprefix /,$@) $(addprefix /,$(filter %.bam,$^)) 1> $@.stdout 2> $@.stderr

step10: $(sort $(addprefix Output/WES/Depth/,$(addsuffix .tsv,$(WES_ID))))
.PHONY += step10

# Step 11 (Sequenza: Make reference)
Output/Sequenza:
	mkdir -p $@

Output/Sequenza/hg38.wig.gz: $(REFERENCE_FASTA) | Output/Sequenza build.log
	$(DOCKER) sequenza-utils gc_wiggle --fasta $(addprefix /,$^) -o $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

step11: Output/Sequenza/hg38.wig.gz
.PHONY += step11

# Step 12 (Sequenza: Normal - Tumor pairs)

# Step ## (WTS - fq.gz >> bam) - STAR
# Step ## (R-sem)
# Step ## (WTS - samtools depth)
# Step ## (Picard - CollectMultipleMetrics)
# Step ## (PureCN)
