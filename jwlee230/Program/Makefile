# Variables
DATE := $(shell date "+%y%m%d")
RANDOM := $(shell bash -c 'echo $$RANDOM')
IMAGE_NAME = lungcancer:latest
CPUS = 20
MEMS = 100G
PWD := $(shell pwd)
TOOLS = 
DOCKER = docker run --rm $(VOLUME_OPTS) $(RUN_OPTS) $(IMAGE_NAME)

# Options
VOLUME_OPTS = --volume $(abspath Output):/Output --volume $(abspath Data):/Data:ro --volume $(abspath Python):/Python
RUN_OPTS = --tty --cpus="$(CPUS)" --memory="$(MEMS)"
JAVA_OPTS = -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Xmx$(MEMS)
GATK_OPTS = --java-options "$(JAVA_OPTS)"
PICARD = java $(JAVA_OPTS) -jar /Tools/picard.jar

# Files
REFERENCE_DIR = Data/hg38
REFERENCE_FASTA = $(REFERENCE_DIR)/Homo_sapiens_assembly38.fasta
REFERENCE_KNOWN_VCF = $(wildcard $(REFERENCE_DIR)/*.vcf)
WGS_ID = $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | cut --delimiter _ --fields 1)

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
	wget "https://github.com/broadinstitute/picard/releases/download/2.24.0/picard.jar" -O $@
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
	echo "make -C $(PWD) latest" > tmp.sh
	qsub -cwd -l h_vmem=$(MEMS) -m abe -M "230@fumire.moe" -N Lung_$(DATE) -pe smp $(CPUS) -o $(abspath log) -e $(abspath log) tmp.sh
.PHONY += run

parallel_run: | log Output build.log
	echo "make -j $(CPUS) -C $(PWD) CPUS=1 MEMS=5G latest" > tmp.sh
	qsub -cwd -l h_vmem=$(MEMS) -m abe -M "230@fumire.moe" -N Lung_$(DATE) -pe smp $(CPUS) -o $(abspath log) -e $(abspath log) tmp.sh
.PHONY += parallel_run

# Actual
latest: step03
.PHONY += latest

# Step 01 (FastQC for WGS)
Output/WGS/FastQC:
	mkdir -p $@

Output/WGS/FastQC/%_fastqc.zip: Data/RawData/%.fastq.gz | Output/WGS/FastQC build.log
	$(DOCKER) fastqc --outdir $(addprefix /,$(@D)) --format fastq --noextract $(addprefix /,$^) 1> $@.stdout 2> $@.stderr

Output/WGS/FastQC/summary.png: Python/aggregate_fastqc.py $(addprefix Output/WGS/FastQC/,$(addsuffix _fastqc.zip,$(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | cut --delimiter . --fields 1))))
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step01: Output/WGS/FastQC/summary.png
.PHONY += step01

# Step 02 (FastQC for WTS)
Output/WTS/FastQC:
	mkdir -p $@

Output/WTS/FastQC/%_fastqc.zip: Data/WTS/%.fastq.gz | Output/WTS/FastQC build.log
	$(DOCKER) fastqc --outdir $(addprefix /,$(@D)) --format fastq --noextract $(addprefix /,$^) 1> $@.stdout 2> $@.stderr

Output/WTS/FastQC/summary.png: Python/aggregate_fastqc.py $(addprefix Output/WTS/FastQC/,$(addsuffix _fastqc.zip,$(sort $(shell ls -1 Data/WTS/*.fastq.gz | cut --delimiter / --fields 3 | cut --delimiter . --fields 1))))
	$(DOCKER) python3 $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step02: Output/WTS/FastQC/summary.png
.PHONY += step02

# Step 03 (WGS - fq.gz >> bam; BAM with multi-threading)
Output/WGS/BAM:
	mkdir -p $@

Output/WGS/BAM/%.BWA.sam: $(REFERENCE_FASTA) Data/RawData/%_*.fastq.gz | Output/WGS/BAM build.log
	$(DOCKER) bwa mem -M -t $(CPUS) -R "@RG\tID:$(@F)\tPL:ILLUMINA\tLB:$(@F)\tSM:$(@F)\tCN:UNIST" -v 3 $(addprefix /,$(filter %.fasta,$^)) $(sort $(addprefix /,$(filter %.fastq.gz,$^))) 1> $@ 2> $@.stderr
	sed --in-place "/^\[.*/d" $@

Output/WGS/BAM/%.BWA.bam: $(REFERENCE_FASTA) Output/WGS/BAM/%.BWA.sam | build.log
	$(DOCKER) samtools view -b -h -S --threads $(CPUS) --reference $(addprefix /,$(filter %.fasta,$^)) $(addprefix /,$(filter %.sam,$^)) -o $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/WGS/BAM/%.Bowtie2.bam: $(REFERENCE_FASTA) Data/RawData/%_*.fastq.gz | Output/WGS/BAM build.log
	$(DOCKER) bowtie2 -x $(addprefix /,$(filter %.fasta,$^)) -1 $(addprefix /,$(word 1,$(sort $(filter %.fastq.gz,$^)))) -2 $(addprefix /,$(word 2,$(sort $(filter %.fastq.gz,$^)))) -b $(addprefix /,$@) --threads $(CPUS) --rg-id $(@F) --preserve-tags

Output/WGS/BAM/%.st.bam: $(REFERENCE_FASTA) Output/WGS/BAM/%.bam | build.log
	$(DOCKER) samtools sort -l 9 -m $(MEMS) --threads $(CPUS) --reference $(addprefix /,$(filter %.fasta,$^)) -o $(addprefix /,$@) $(addprefix /,$(filter %.bam,$^)) 1> $@.stdout 2> $@.stderr

Output/WGS/BAM/%.st.bam.bai: Output/WGS/BAM/%.st.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/WGS/BAM/%.st.dm.bam: $(REFERENCE_FASTA) Output/WGS/BAM/%.st.bam Output/WGS/BAM/%.st.bam.bai | build.log
	$(DOCKER) $(PICARD) MarkDuplicates --INPUT=$(addprefix /,$(filter %.bam,$^)) --REFERENCE_SEQUENCE=$(addprefix /,$(filter %.fasta,$^)) --OUTPUT=$(addprefix /,$@) --METRICS_FILE $(addprefix /,$@.metrics) --ASSUME_SORTED=true --COMPRESSION_LEVEL=9 --CREATE_INDEX=true --VALIDATION_STRINGENCY=SILENT --VERBOSITY=INFO 1> $@.stdout 2> $@.stderr

Output/WGS/BAM/%.st.dm.bam.bai: Output/WGS/BAM/%.st.dm.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step03: $(addprefix Output/WGS/BAM/,$(addsuffix .BWA.st.dm.bam,$(WGS_ID))) $(addprefix Output/WGS/BAM/,$(addsuffix .BWA.st.dm.bam.bai,$(WGS_ID))) $(addprefix Output/WGS/BAM/,$(addsuffix .Bowtie2.st.dm.bam,$(WGS_ID))) $(addprefix Output/WGS/BAM/,$(addsuffix .Bowtie2.st.dm.bam.bai,$(WGS_ID)))
.PHONY += step03

clean-step03:
	$(DOCKER) rm -fv $(filter-out %.st.dm.bam,$(shell find /Output/WGS/BAM -size +1G))
.PHONY += clean-step03

# Step 04 (WGS - fq.gz >> bam; single-threading)
Output/WGS/BAM/%.st.dm.rc.table: $(REFERENCE_FASTA) $(REFERENCE_KNOWN_VCF) Output/WGS/BAM/%.st.dm.bam Output/WGS/BAM/%.st.dm.bam.bai | build.log
	$(DOCKER) gatk $(GATK_OPTS) BaseRecalibrator --input $(addprefix /,$(filter %.bam,$^)) $(foreach site,$(REFERENCE_KNOWN_VCF),--known-sites $(site)) --output $(addprefix /,$@) --reference $(addprefix /,$(filter %.fasta,$^)) 1> $@.stdout 2> $@.stderr

Output/WGS/BAM/%.st.dm.rc.bam: Output/WGS/BAM/%.st.dm.bam Output/WGS/BAM/%.st.dm.rc.table $(REFERENCE_FASTA) | build.log
	$(DOCKER) gatk $(GATK_OPTS) ApplyBQSR --bqsr-recal-file $(addprefix /,$(filter %.table,$^)) --input $(addprefix /,$(filter %.bam,$^)) --output $(addprefix /,$@) --reference $(addprefix /,$(filter %.fasta,$^)) 1> $@.stdout 2> $@.stderr

step04: $(addprefix Output/WGS/BAM/,$(addsuffix .BWA.st.dm.rc.bam,$(WGS_ID))) $(addprefix Output/WGS/BAM/,$(addsuffix .Bowtie2.st.dm.rc.bam,$(WGS_ID)))
.PHONY += step04
