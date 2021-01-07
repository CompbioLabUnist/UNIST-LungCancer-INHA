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
VOLUME_OPTS = --volume $(abspath Output):/Output --volume $(abspath Data):/Data --volume $(abspath Python):/Python
RUN_OPTS = --tty --cpus="$(CPUS)" --memory="$(MEMS)"
GATK_OPTS = --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=1 -Xmx$(MEMS)"
PICARD = /Tools/picard.jar

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
	qsub -cwd -l h_vmem=$(MEMS) -m abe -M "230@fumire.moe" -N Colon_$(DATE) -pe smp $(CPUS) -o $(abspath log) -e $(abspath log) tmp.sh
.PHONY += run

parallel_run: | log Output build.log
	echo "make -j $(CPUS) -C $(PWD) latest" > tmp.sh
	qsub -cwd -l h_vmem=$(MEMS) -m abe -M "230@fumire.moe" -N Colon_$(DATE) -pe smp $(CPUS) -o $(abspath log) -e $(abspath log) tmp.sh
.PHONY += parallel_run

# Actual2021-01-07T14:33:40+0900
latest: step00
.PHONY += latest
