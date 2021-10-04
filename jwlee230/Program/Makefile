# Variables
DATE := $(shell date "+%y%m%d")
RANDOM := $(shell bash -c 'echo $$RANDOM')
IMAGE_NAME = lungcancer:latest
SHELL := /bin/bash
CPUS = 10
MEMS = 100
PARALLEL = 10
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
CLINICAL_DATA_CSV = Data/Clinical.csv
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
GENE_CENSUS_CSV = Data/CGC_tier1.csv
LUNG_SINGLECELL_ANNOTATION_TXT = Data/Lung_SingleCell/GSE131907_Lung_Cancer_cell_annotation.txt
LUNG_SINGLECELL_RAW_MATRIX_TXT = Data/Lung_SingleCell/GSE131907_Lung_Cancer_raw_UMI_matrix.txt
LUNG_SINGLECELL_RAW_MATRIX_RDS = Data/Lung_SingleCell/GSE131907_Lung_Cancer_raw_UMI_matrix.rds
LUNG_SINGLECELL_NORMALIZED_MATRIX_TXT = Data/Lung_SingleCell/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt
LUNG_SINGLECELL_NORMALIZED_MATRIX_RDS = Data/Lung_SingleCell/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds
KNOWN_GENE_GTF = Data/hg38.gtf.gz
ARRIBA_STAR_DIR = Tools/arriba_v2.1.0/STAR_index_hg38_RefSeq_hg38
ARRIBA_DATABASE_DIR = Tools/arriba_v2.1.0/database
ARRIBA_BLACKLIST_TSV = $(ARRIBA_DATABASE_DIR)/blacklist_hg38_GRCh38_v2.1.0.tsv.gz
ARRIBA_KNOWN_FUSIONS_TSV = $(ARRIBA_DATABASE_DIR)/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz
ARRIBA_PROTEIN_DOMAIN_GFF3 = $(ARRIBA_DATABASE_DIR)/protein_domains_hg38_GRCh38_v2.1.0.gff3
CTAT_RESOURCE_LIB_DIR = Data/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir

# Sample IDs
PATIENT_ID := $(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+"))
WES_ID := $(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+[PCADMN]\d?"))
WES_NORMAL_ID := $(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+N"))
WES_TUMOR_ID := $(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+[PCADM]\d?"))
WTS_ID := $(sort $(shell ls -1 Data/WTS/*.fastq.gz | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+[PCADMN]"))
WTS_NORMAL_ID := $(filter %N,$(WTS_ID))
WTS_TUMOR_ID := $(filter-out %N,$(WTS_ID))

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

Docker/Tools/Spruce.zip: | Docker/Tools
	wget "https://github.com/raphael-group/spruce/archive/refs/heads/master.zip" -O $@
TOOLS += Docker/Tools/Spruce.zip

Docker/Tools/QAPA.zip: | Docker/Tools
	wget "https://codeload.github.com/morrislab/qapa/zip/refs/heads/master" -O $@
TOOLS += Docker/Tools/QAPA.zip

Docker/Tools/annovar.latest.tar.gz: | Docker/Tools
	@echo "You must manullay download <annovar> on https://www.openbioinformatics.org/annovar/annovar_download_form.php"
	touch $@
TOOLS += Docker/Tools/annovar.latest.tar.gz

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
latest: step03 step04 step05 step06 step07 step08 step09 step10 step11 step12 step13 step14 step15 step16 step17 step18 step19 step20 step21 step22 step23 step24 step25 step26 step27 step28 step29 step30
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
Output/Sequenza/BWA/%/sample_alternative_solutions.txt Output/Sequenza/BWA/%/sample_segments.txt: R/run_sequenza.R Output/Sequenza/BWA/%.seqz.gz Output/Sequenza/BWA/%.seqz.gz.tbi | build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D)) && mkdir -p $(@D)
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --sampleid "sample" --input $(addprefix /,$(filter %.seqz.gz,$^)) --outdir $(addprefix /,$(@D)) --parallel $(CPUS) 1> $(@D).stdout 2> $(@D).stderr

Output/Sequenza/Bowtie2/%/sample_alternative_solutions.txt Output/Sequenza/BWA/%/sample_segments.txt: R/run_sequenza.R Output/Sequenza/Bowtie2/%.seqz.gz Output/Sequenza/Bowtie2/%.seqz.gz.tbi | build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D)) && mkdir -p $(@D)
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --sampleid "sample" --input $(addprefix /,$(filter %.seqz.gz,$^)) --outdir $(addprefix /,$(@D)) --parallel $(CPUS) 1> $(@D).stdout 2> $(@D).stderr

Output/Sequenza/BWA-sequenza-SQC.pdf: Python/aggregate_sequenza.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-sequenza-ADC.pdf: Python/aggregate_sequenza.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-sequenza-SQC.pdf: Python/aggregate_sequenza.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-sequenza-ADC.pdf: Python/aggregate_sequenza.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-SQC.pdf: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-ADC.pdf: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-genome-SQC.pdf: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-genome-ADC.pdf: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-simple-SQC.pdf: Python/aggregate_sequenza_simple.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-simple-ADC.pdf: Python/aggregate_sequenza_simple.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-simple-SQC.pdf: Python/aggregate_sequenza_simple.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-simple-ADC.pdf: Python/aggregate_sequenza_simple.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequnza/BWA-correlation-SQC.pdf: Python/aggregate_sequenza_correlation.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(REFERENCE_CENTROMERE) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC

step13-Sequenza: $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID)))
step13-default: Output/Sequenza/BWA-sequenza-SQC.pdf Output/Sequenza/BWA-sequenza-ADC.pdf Output/Sequenza/Bowtie2-sequenza-SQC.pdf Output/Sequenza/Bowtie2-sequenza-ADC.pdf
step13-genome: Output/Sequenza/BWA-genome-SQC.pdf Output/Sequenza/BWA-genome-ADC.pdf Output/Sequenza/Bowtie2-genome-SQC.pdf Output/Sequenza/Bowtie2-genome-ADC.pdf
step13-simple: Output/Sequenza/BWA-simple-SQC.pdf Output/Sequenza/BWA-simple-ADC.pdf Output/Sequenza/Bowtie2-simple-SQC.pdf Output/Sequenza/Bowtie2-simple-ADC.pdf
step13: step13-Sequenza step13-default step13-genome step13-simple
.PHONY: step13-Sequenza step13-default step13-genome step13-simple step13

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
Output/MutEnricher Output/MutEnricher/Reference Output/MutEnricher/BWA Output/MutEnricher/Bowtie2:
	mkdir -p $@

Output/MutEnricher/Reference/annovar_downdb.log Output/MutEnricher/Reference/hg38_refGene.txt Output/MutEnricher/Reference/hg38_refGeneMrna.fa Output/MutEnricher/Reference/hg38_refGeneVersion.txt: | Output/MutEnricher/Reference build.log
	$(DOCKER) annotate_variation.pl --verbose --buildver hg38 --downdb --webfrom annovar refGene $(addprefix /,$(@D)) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/BWA/%.hg38_multianno.vcf: Output/Mutect2/BWA/%.BWA.PASS.vcf Output/MutEnricher/Reference/hg38_refGeneMrna.fa | build.log Output/MutEnricher/BWA
	$(DOCKER) table_annovar.pl $(addprefix /,$(filter %.vcf,$^)) $(addprefix /,$(dir $(filter %.fa,$^))) --verbose --buildver "hg38" --remove --vcfinput --protocol "refGene" --operation g --nastring "." --outfile $(addprefix /,$(@D)/$*) --thread $(CPUS) --maxgenethread $(CPUS) --tempdir $(TMPFS) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/BWA/%.vcf.gz: Bash/prepare_VCF_annotation.bash Output/MutEnricher/BWA/%.hg38_multianno.vcf | build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@)

Output/MutEnricher/BWA/%.vcf.gz.tbi: Output/MutEnricher/BWA/%.vcf.gz | build.log
	$(DOCKER) tabix --preset "vcf" --force $(addprefix /,$^) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/BWA/SQC.tsv: Python/make_list_MutEnricher.py $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz,$(WES_TUMOR_ID))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz.tbi,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^) $(filter %.vcf.gz,$^) $(filter %.csv,$^) $@) --SQC 1> $@.stdout 2> $@.stderr

Output/MutEnricher/BWA/ADC.tsv: Python/make_list_MutEnricher.py $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz,$(WES_TUMOR_ID))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz.tbi,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^) $(filter %.vcf.gz,$^) $(filter %.csv,$^) $@) --ADC 1> $@.stdout 2> $@.stderr

Output/MutEnricher/Bowtie2/%.hg38_multianno.vcf: Output/Mutect2/Bowtie2/%.Bowtie2.PASS.vcf Output/MutEnricher/Reference/hg38_refGeneMrna.fa | build.log Output/MutEnricher/Bowtie2
	$(DOCKER) table_annovar.pl $(addprefix /,$(filter %.vcf,$^)) $(addprefix /,$(dir $(filter %.fa,$^))) --verbose --buildver "hg38" --remove --vcfinput --protocol "refGene" --operation g --nastring "." --outfile $(addprefix /,$(@D)/$*) --thread $(CPUS) --maxgenethread $(CPUS) --tempdir $(TMPFS) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/Bowtie2/%.vcf.gz: Bash/prepare_VCF_annotation.bash Output/MutEnricher/Bowtie2/%.hg38_multianno.vcf | build.log Output/MutEnricher/Bowtie2
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/Bowtie2/%.vcf.gz.tbi: Output/MutEnricher/Bowtie2/%.vcf.gz | build.log
	$(DOCKER) tabix --preset "vcf" --force $(addprefix /,$^) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/Bowtie2/SQC.tsv: Python/make_list_MutEnricher.py $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .vcf.gz,$(WES_TUMOR_ID))) $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .vcf.gz.tbi,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^) $(filter %.vcf.gz,$^) $(filter %.csv,$^) $@) --SQC 1> $@.stdout 2> $@.stderr

Output/MutEnricher/Bowtie2/ADC.tsv: Python/make_list_MutEnricher.py $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .vcf.gz,$(WES_TUMOR_ID))) $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .vcf.gz.tbi,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^) $(filter %.vcf.gz,$^) $(filter %.csv,$^) $@) --ADC 1> $@.stdout 2> $@.stderr

Output/MutEnricher/gene_covariates.txt: $(KNOWN_GENE_GTF) $(REFERENCE_FASTA) | build.log
	$(DOCKER) $(PYTHON) /Tools/MutEnricher-1.3.1/utilities/get_gene_covariates.py --outname $(addprefix /,$@) --processors $(CPUS) $(addprefix /,$^) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/BWA/%/MutEnricher_gene_enrichments.txt Output/MutEnricher/BWA/%/MutEnricher_hotspot.txt Output/MutEnricher/BWA/%/MutEnricher_gene_hotspot_Fisher_enrichments.txt Output/MutEnricher/BWA/%/MutEnricher_gene_data.pkl Output/MutEnricher/BWA/%/MutEnricher.log: Output/MutEnricher/gene_covariates.txt $(KNOWN_GENE_GTF) Output/MutEnricher/BWA/%.tsv | build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D)) && mkdir -p $(@D)
	$(DOCKER) $(MUTENRICHER) coding --outdir $(addprefix /,$(@D)) --prefix "MutEnricher" --anno-type "annovar-refGene" --exome-only --by-contig --processors $(CPUS) --covariates-file $(addprefix /,$(filter %.txt,$^)) $(addprefix /,$(filter %.gtf.gz,$^)) $(addprefix /,$(filter %.tsv,$^)) 1> $(@D).stdout 2> $(@D).stderr

Output/MutEnricher/Bowtie2/%/MutEnricher_gene_enrichments.txt Output/MutEnricher/Bowtie2/%/MutEnricher_hotspot.txt Output/MutEnricher/Bowtie2/%/MutEnricher_gene_hotspot_Fisher_enrichments.txt Output/MutEnricher/Bowtie2/%/MutEnricher_gene_data.pkl Output/MutEnricher/Bowtie2/%/MutEnricher.log: Output/MutEnricher/gene_covariates.txt $(KNOWN_GENE_GTF) Output/MutEnricher/Bowtie2/%.tsv | build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D)) && mkdir -p $(@D)
	$(DOCKER) $(MUTENRICHER) coding --outdir $(addprefix /,$(@D)) --prefix "MutEnricher" --anno-type "annovar-refGene" --exome-only --by-contig --processors $(CPUS) --covariates-file $(addprefix /,$(filter %.txt,$^)) $(addprefix /,$(filter %.gtf.gz,$^)) $(addprefix /,$(filter %.tsv,$^)) 1> $(@D).stdout 2> $(@D).stderr

Output/Mutect2/BWA-SQC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-SQC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-ADC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-ADC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

step16-preprocess: Output/MutEnricher/gene_covariates.txt $(addprefix Output/MutEnricher/BWA/,$(addsuffix .hg38_multianno.vcf,$(WES_TUMOR_ID))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz,$(WES_TUMOR_ID))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz.tbi,$(WES_TUMOR_ID))) Output/MutEnricher/BWA/SQC.tsv Output/MutEnricher/BWA/ADC.tsv $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .hg38_multianno.vcf,$(WES_TUMOR_ID))) $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .vcf.gz,$(WES_TUMOR_ID))) $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .vcf.gz.tbi,$(WES_TUMOR_ID))) Output/MutEnricher/Bowtie2/SQC.tsv Output/MutEnricher/Bowtie2/ADC.tsv
step16-MutEnricher: Output/MutEnricher/BWA/SQC/MutEnricher.log Output/MutEnricher/BWA/ADC/MutEnricher.log Output/MutEnricher/Bowtie2/ADC/MutEnricher.log Output/MutEnricher/Bowtie2/SQC/MutEnricher.log
step16-plot: Output/Mutect2/BWA-SQC.pdf Output/Mutect2/BWA-ADC.pdf Output/Mutect2/Bowtie2-SQC.pdf Output/Mutect2/Bowtie2-ADC.pdf
step16: step16-preprocess step16-MutEnricher step16-plot
.PHONY: step16-preprocess step16-MutEnricher step16-plot step16

# Step 17 (RSEM: Prepare reference genome)
Output/RSEM/Reference:
	mkdir -p $@

Output/RSEM/Reference/hg38.gtf.gz: | Output/RSEM/Reference
	wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz" -O $@ 1> $@.stdout 2> $@.stderr
	touch $@

Output/RSEM/Reference/hg38.gtf: Output/RSEM/Reference/hg38.gtf.gz | build.log
	$(DOCKER) gzip --decompress --keep --force $(addprefix /,$<) 1> $@.stdout 2> $@.stderr

Output/RSEM/Reference/hg38.seq Output/RSEM/Reference/hg38.grp: $(REFERENCE_FASTA) Output/RSEM/Reference/hg38.gtf | build.log
	$(DOCKER) rsem-prepare-reference --gtf $(addprefix /,$(filter %.gtf,$^)) --bowtie2 --star --num-threads $(CPUS) $(addprefix /,$(filter %.fasta,$^)) $(addprefix /,$(basename $@)) 1> $@.stdout 2> $@.stderr

step17: Output/RSEM/Reference/hg38.gtf.gz Output/RSEM/Reference/hg38.gtf Output/RSEM/Reference/hg38.seq
.PHONY: step17

# Step 18 (RSEM: Calculate expression)
Output/RSEM/Bowtie2 Output/RSEM/STAR Output/RSEM/DEG:
	mkdir -p $@

Output/RSEM/Bowtie2/%.genome.bam Output/RSEM/Bowtie2/%.genes.results: Data/WTS/%_*.fastq.gz Output/RSEM/Reference/hg38.grp | Output/RSEM/Bowtie2 build.log
	$(DOCKER) rsem-calculate-expression --num-threads $(CPUS) --estimate-rspd --ci-memory $(shell echo "$(MEMS) * 1024" | bc) --bowtie2 --output-genome-bam --paired-end $(addprefix /,$(sort $(filter %.fastq.gz,$^)) $(basename $(filter %.grp,$^)) $(@D)/$*) 1> $@.stdout 2> $@.stderr

Output/RSEM/STAR/%.genome.bam Output/RSEM/STAR/%.genes.results: Data/WTS/%_*.fastq.gz Output/RSEM/Reference/hg38.grp | Output/RSEM/STAR build.log
	$(DOCKER) rsem-calculate-expression --num-threads $(CPUS) --estimate-rspd --ci-memory $(shell echo "$(MEMS) * 1024" | bc) --star --star-gzipped-read-file --output-genome-bam --paired-end $(addprefix /,$(sort $(filter %.fastq.gz,$^)) $(basename $(filter %.grp,$^)) $(@D)/$*) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/Bowtie2.tsv: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(WTS_ID))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/STAR.tsv: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(WTS_ID))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Normal-Dysplasia.tsv Output/RSEM/DEG/%.SQC.Normal-Dysplasia.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --compare "stage" "Normal" "Dysplasia" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Normal-CIS.tsv Output/RSEM/DEG/%.SQC.Normal-CIS.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --compare "stage" "Normal" "CIS+AIS" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Normal-Primary.tsv Output/RSEM/DEG/%.SQC.Normal-Primary.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --compare "stage" "Normal" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Dysplasia-CIS.tsv Output/RSEM/DEG/%.SQC.Dysplasia-CIS.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --compare "stage" "Dysplasia" "CIS+AIS" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Dysplasia-Primary.tsv Output/RSEM/DEG/%.SQC.Dysplasia-Primary.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --compare "stage" "Dysplasia" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.CIS-Primary.tsv Output/RSEM/DEG/%.SQC.CIS-Primary.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --compare "stage" "CIS+AIS" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.ADC.Normal-AAH.tsv Output/RSEM/DEG/%.ADC.Normal-AAH.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --ADC --compare "stage" "Normal" "AAH" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.ADC.Normal-AIS.tsv Output/RSEM/DEG/%.ADC.Normal-AIS.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --ADC --compare "stage" "Normal" "CIS+AIS" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.ADC.Normal-Primary.tsv Output/RSEM/DEG/%.ADC.Normal-Primary.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --ADC --compare "stage" "Normal" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.ADC.AAH-AIS.tsv Output/RSEM/DEG/%.ADC.AAH-AIS.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --ADC --compare "stage" "AAH" "CIS+AIS" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.ADC.AAH-Primary.tsv Output/RSEM/DEG/%.ADC.AAH-Primary.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --ADC --compare "stage" "AAH" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.ADC.AIS-Primary.tsv Output/RSEM/DEG/%.ADC.AIS-Primary.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --ADC --compare "stage" "CIS+AIS" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.coldata: Python/aggregate_RSEM_3.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --exact "Recurrence" "YES" --compare "stage" "Normal" "CIS+AIS" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.coldata: Python/aggregate_RSEM_3.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --exact "Recurrence" "YES" --compare "stage" "Normal" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Recur.CIS-Primary.tsv Output/RSEM/DEG/%.SQC.Recur.CIS-Primary.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --exact "Recurrence" "YES" --compare "stage" "CIS+AIS" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.coldata: Python/aggregate_RSEM_3.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --exact "Recurrence" "NO" --compare "stage" "Normal" "CIS+AIS" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.coldata: Python/aggregate_RSEM_3.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --exact "Recurrence" "NO" --compare "stage" "Normal" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Nonrecur.CIS-Primary.tsv Output/RSEM/DEG/%.SQC.NonrecurCIS-Primary.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --exact "Recurrence" "NO" --compare "stage" "CIS+AIS" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.ADC.Nonrecur.Normal-AIS.tsv Output/RSEM/DEG/%.ADC.Nonrecur.Normal-AIS.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --ADC --exact "Recurrence" "NO" --compare "stage" "Normal" "CIS+AIS" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.ADC.Nonrecur.Normal-Primary.tsv Output/RSEM/DEG/%.ADC.Nonrecur.Normal-Primary.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --ADC --exact "Recurrence" "NO" --compare "stage" "Normal" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.ADC.Nonrecur.AIS-Primary.tsv Output/RSEM/DEG/%.ADC.Nonrecur.AIS-Primary.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --ADC --exact "Recurrence" "NO" --compare "stage" "CIS+AIS" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.CIS.Nonrecur-Recur.tsv Output/RSEM/DEG/%.SQC.CIS.Nonrecur-Recur.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --exact "stage" "CIS+AIS" --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Primary.Nonrecur-Recur.tsv Output/RSEM/DEG/%.SQC.Primary.Nonrecur-Recur.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --exact "stage" "Primary" --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.DEG.tsv: R/run_DESeq.R Output/RSEM/DEG/%.tsv Output/RSEM/DEG/%.coldata | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --readcount $(addprefix /,$(filter %.tsv,$^)) --coldata $(addprefix /,$(filter %.coldata,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.Up.tex: Python/select_DEGs.py Output/RSEM/DEG/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --up 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.Down.tex: Python/select_DEGs.py Output/RSEM/DEG/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --down 1> $@.stdout 2> $@.stderr

step18-preparation: $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(WTS_ID))) Output/RSEM/DEG/Bowtie2.tsv $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(WTS_ID))) Output/RSEM/DEG/STAR.tsv
step18-SQC: Output/RSEM/DEG/Bowtie2.SQC.Normal-Dysplasia.tsv Output/RSEM/DEG/Bowtie2.SQC.Normal-CIS.tsv Output/RSEM/DEG/Bowtie2.SQC.Normal-Primary.tsv Output/RSEM/DEG/Bowtie2.SQC.Dysplasia-CIS.tsv Output/RSEM/DEG/Bowtie2.SQC.Dysplasia-Primary.tsv Output/RSEM/DEG/Bowtie2.SQC.CIS-Primary.tsv Output/RSEM/DEG/STAR.SQC.Normal-Dysplasia.tsv Output/RSEM/DEG/STAR.SQC.Normal-CIS.tsv Output/RSEM/DEG/STAR.SQC.Normal-Primary.tsv Output/RSEM/DEG/STAR.SQC.Dysplasia-CIS.tsv Output/RSEM/DEG/STAR.SQC.Dysplasia-Primary.tsv Output/RSEM/DEG/STAR.SQC.CIS-Primary.tsv
step18-SQC-DEG: Output/RSEM/DEG/Bowtie2.SQC.Normal-Dysplasia.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Normal-CIS.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Normal-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Dysplasia-CIS.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Dysplasia-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.CIS-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.Normal-Dysplasia.DEG.tsv Output/RSEM/DEG/STAR.SQC.Normal-CIS.DEG.tsv Output/RSEM/DEG/STAR.SQC.Normal-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.Dysplasia-CIS.DEG.tsv Output/RSEM/DEG/STAR.SQC.Dysplasia-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.CIS-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Normal-Dysplasia.DEG.Up.tex Output/RSEM/DEG/Bowtie2.SQC.Normal-CIS.DEG.Up.tex Output/RSEM/DEG/Bowtie2.SQC.Normal-Primary.DEG.Up.tex Output/RSEM/DEG/Bowtie2.SQC.Dysplasia-CIS.DEG.Up.tex Output/RSEM/DEG/Bowtie2.SQC.Dysplasia-Primary.DEG.Up.tex Output/RSEM/DEG/Bowtie2.SQC.CIS-Primary.DEG.Up.tex Output/RSEM/DEG/STAR.SQC.Normal-Dysplasia.DEG.Up.tex Output/RSEM/DEG/STAR.SQC.Normal-CIS.DEG.Up.tex Output/RSEM/DEG/STAR.SQC.Normal-Primary.DEG.Up.tex Output/RSEM/DEG/STAR.SQC.Dysplasia-CIS.DEG.Up.tex Output/RSEM/DEG/STAR.SQC.Dysplasia-Primary.DEG.Up.tex Output/RSEM/DEG/STAR.SQC.CIS-Primary.DEG.Up.tex Output/RSEM/DEG/Bowtie2.SQC.Normal-Dysplasia.DEG.Down.tex Output/RSEM/DEG/Bowtie2.SQC.Normal-CIS.DEG.Down.tex Output/RSEM/DEG/Bowtie2.SQC.Normal-Primary.DEG.Down.tex Output/RSEM/DEG/Bowtie2.SQC.Dysplasia-CIS.DEG.Down.tex Output/RSEM/DEG/Bowtie2.SQC.Dysplasia-Primary.DEG.Down.tex Output/RSEM/DEG/Bowtie2.SQC.CIS-Primary.DEG.Down.tex Output/RSEM/DEG/STAR.SQC.Normal-Dysplasia.DEG.Down.tex Output/RSEM/DEG/STAR.SQC.Normal-CIS.DEG.Down.tex Output/RSEM/DEG/STAR.SQC.Normal-Primary.DEG.Down.tex Output/RSEM/DEG/STAR.SQC.Dysplasia-CIS.DEG.Down.tex Output/RSEM/DEG/STAR.SQC.Dysplasia-Primary.DEG.Down.tex Output/RSEM/DEG/STAR.SQC.CIS-Primary.DEG.Down.tex
step18-SQC-Recur: Output/RSEM/DEG/Bowtie2.SQC.Recur.Normal-CIS.tsv Output/RSEM/DEG/Bowtie2.SQC.Recur.Normal-Primary.tsv Output/RSEM/DEG/Bowtie2.SQC.Recur.CIS-Primary.tsv Output/RSEM/DEG/STAR.SQC.Recur.Normal-CIS.tsv Output/RSEM/DEG/STAR.SQC.Recur.Normal-Primary.tsv Output/RSEM/DEG/STAR.SQC.Recur.CIS-Primary.tsv Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.Normal-CIS.tsv Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.Normal-Primary.tsv Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.CIS-Primary.tsv Output/RSEM/DEG/STAR.SQC.Nonrecur.Normal-CIS.tsv Output/RSEM/DEG/STAR.SQC.Nonrecur.Normal-Primary.tsv Output/RSEM/DEG/STAR.SQC.Nonrecur.CIS-Primary.tsv
step18-SQC-Recur-DEG: Output/RSEM/DEG/Bowtie2.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Recur.CIS-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/STAR.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.Recur.CIS-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.CIS-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/STAR.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.Nonrecur.CIS-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Recur.Normal-CIS.DEG.Up.tex Output/RSEM/DEG/Bowtie2.SQC.Recur.Normal-Primary.DEG.Up.tex Output/RSEM/DEG/Bowtie2.SQC.Recur.CIS-Primary.DEG.Up.tex Output/RSEM/DEG/STAR.SQC.Recur.Normal-CIS.DEG.Up.tex Output/RSEM/DEG/STAR.SQC.Recur.Normal-Primary.DEG.Up.tex Output/RSEM/DEG/STAR.SQC.Recur.CIS-Primary.DEG.Up.tex Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.Normal-CIS.DEG.Up.tex Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.Normal-Primary.DEG.Up.tex Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.CIS-Primary.DEG.Up.tex Output/RSEM/DEG/STAR.SQC.Nonrecur.Normal-CIS.DEG.Up.tex Output/RSEM/DEG/STAR.SQC.Nonrecur.Normal-Primary.DEG.Up.tex Output/RSEM/DEG/STAR.SQC.Nonrecur.CIS-Primary.DEG.Up.tex Output/RSEM/DEG/Bowtie2.SQC.Recur.Normal-CIS.DEG.Down.tex Output/RSEM/DEG/Bowtie2.SQC.Recur.Normal-Primary.DEG.Down.tex Output/RSEM/DEG/Bowtie2.SQC.Recur.CIS-Primary.DEG.Down.tex Output/RSEM/DEG/STAR.SQC.Recur.Normal-CIS.DEG.Down.tex Output/RSEM/DEG/STAR.SQC.Recur.Normal-Primary.DEG.Down.tex Output/RSEM/DEG/STAR.SQC.Recur.CIS-Primary.DEG.Down.tex Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.Normal-CIS.DEG.Down.tex Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.Normal-Primary.DEG.Down.tex Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.CIS-Primary.DEG.Down.tex Output/RSEM/DEG/STAR.SQC.Nonrecur.Normal-CIS.DEG.Down.tex Output/RSEM/DEG/STAR.SQC.Nonrecur.Normal-Primary.DEG.Down.tex Output/RSEM/DEG/STAR.SQC.Nonrecur.CIS-Primary.DEG.Down.tex
step18-SQC-Stage: Output/RSEM/DEG/Bowtie2.SQC.CIS.Nonrecur-Recur.tsv Output/RSEM/DEG/Bowtie2.SQC.Primary.Nonrecur-Recur.tsv Output/RSEM/DEG/STAR.SQC.CIS.Nonrecur-Recur.tsv Output/RSEM/DEG/STAR.SQC.Primary.Nonrecur-Recur.tsv
step18-SQC-Stage-DEG: Output/RSEM/DEG/Bowtie2.SQC.CIS.Nonrecur-Recur.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Primary.Nonrecur-Recur.DEG.tsv Output/RSEM/DEG/STAR.SQC.CIS.Nonrecur-Recur.DEG.tsv Output/RSEM/DEG/STAR.SQC.Primary.Nonrecur-Recur.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.CIS.Nonrecur-Recur.DEG.Up.tex Output/RSEM/DEG/Bowtie2.SQC.Primary.Nonrecur-Recur.DEG.Up.tex Output/RSEM/DEG/STAR.SQC.CIS.Nonrecur-Recur.DEG.Up.tex Output/RSEM/DEG/STAR.SQC.Primary.Nonrecur-Recur.DEG.Up.tex Output/RSEM/DEG/Bowtie2.SQC.CIS.Nonrecur-Recur.DEG.Down.tex Output/RSEM/DEG/Bowtie2.SQC.Primary.Nonrecur-Recur.DEG.Down.tex Output/RSEM/DEG/STAR.SQC.CIS.Nonrecur-Recur.DEG.Down.tex Output/RSEM/DEG/STAR.SQC.Primary.Nonrecur-Recur.DEG.Down.tex
step18-ADC: Output/RSEM/DEG/Bowtie2.ADC.Normal-AAH.tsv Output/RSEM/DEG/Bowtie2.ADC.Normal-AIS.tsv Output/RSEM/DEG/Bowtie2.ADC.Normal-Primary.tsv Output/RSEM/DEG/Bowtie2.ADC.AAH-AIS.tsv Output/RSEM/DEG/Bowtie2.ADC.AAH-Primary.tsv Output/RSEM/DEG/Bowtie2.ADC.AIS-Primary.tsv Output/RSEM/DEG/STAR.ADC.Normal-AAH.tsv Output/RSEM/DEG/STAR.ADC.Normal-AIS.tsv Output/RSEM/DEG/STAR.ADC.Normal-Primary.tsv Output/RSEM/DEG/STAR.ADC.AAH-AIS.tsv Output/RSEM/DEG/STAR.ADC.AAH-Primary.tsv Output/RSEM/DEG/STAR.ADC.AIS-Primary.tsv
step18-ADC-DEG: Output/RSEM/DEG/Bowtie2.ADC.Normal-AAH.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.Normal-AIS.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.Normal-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.AAH-AIS.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.AAH-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.AIS-Primary.DEG.tsv Output/RSEM/DEG/STAR.ADC.Normal-AAH.DEG.tsv Output/RSEM/DEG/STAR.ADC.Normal-AIS.DEG.tsv Output/RSEM/DEG/STAR.ADC.Normal-Primary.DEG.tsv Output/RSEM/DEG/STAR.ADC.AAH-AIS.DEG.tsv Output/RSEM/DEG/STAR.ADC.AAH-Primary.DEG.tsv Output/RSEM/DEG/STAR.ADC.AIS-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.Normal-AAH.DEG.Up.tex Output/RSEM/DEG/Bowtie2.ADC.Normal-AIS.DEG.Up.tex Output/RSEM/DEG/Bowtie2.ADC.Normal-Primary.DEG.Up.tex Output/RSEM/DEG/Bowtie2.ADC.AAH-AIS.DEG.Up.tex Output/RSEM/DEG/Bowtie2.ADC.AAH-Primary.DEG.Up.tex Output/RSEM/DEG/Bowtie2.ADC.AIS-Primary.DEG.Up.tex Output/RSEM/DEG/STAR.ADC.Normal-AAH.DEG.Up.tex Output/RSEM/DEG/STAR.ADC.Normal-AIS.DEG.Up.tex Output/RSEM/DEG/STAR.ADC.Normal-Primary.DEG.Up.tex Output/RSEM/DEG/STAR.ADC.AAH-AIS.DEG.Up.tex Output/RSEM/DEG/STAR.ADC.AAH-Primary.DEG.Up.tex Output/RSEM/DEG/STAR.ADC.AIS-Primary.DEG.Up.tex Output/RSEM/DEG/Bowtie2.ADC.Normal-AAH.DEG.Down.tex Output/RSEM/DEG/Bowtie2.ADC.Normal-AIS.DEG.Down.tex Output/RSEM/DEG/Bowtie2.ADC.Normal-Primary.DEG.Down.tex Output/RSEM/DEG/Bowtie2.ADC.AAH-AIS.DEG.Down.tex Output/RSEM/DEG/Bowtie2.ADC.AAH-Primary.DEG.Down.tex Output/RSEM/DEG/Bowtie2.ADC.AIS-Primary.DEG.Down.tex Output/RSEM/DEG/STAR.ADC.Normal-AAH.DEG.Down.tex Output/RSEM/DEG/STAR.ADC.Normal-AIS.DEG.Down.tex Output/RSEM/DEG/STAR.ADC.Normal-Primary.DEG.Down.tex Output/RSEM/DEG/STAR.ADC.AAH-AIS.DEG.Down.tex Output/RSEM/DEG/STAR.ADC.AAH-Primary.DEG.Down.tex Output/RSEM/DEG/STAR.ADC.AIS-Primary.DEG.Down.tex
step18-ADC-Recur: Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.Normal-AIS.tsv Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.Normal-Primary.tsv Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.AIS-Primary.tsv Output/RSEM/DEG/STAR.ADC.Nonrecur.Normal-AIS.tsv Output/RSEM/DEG/STAR.ADC.Nonrecur.Normal-Primary.tsv Output/RSEM/DEG/STAR.ADC.Nonrecur.AIS-Primary.tsv
step18-ADC-Recur-DEG: Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.Normal-AIS.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.AIS-Primary.DEG.tsv Output/RSEM/DEG/STAR.ADC.Nonrecur.Normal-AIS.DEG.tsv Output/RSEM/DEG/STAR.ADC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/STAR.ADC.Nonrecur.AIS-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.Normal-AIS.DEG.Up.tex Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.Normal-Primary.DEG.Up.tex Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.AIS-Primary.DEG.Up.tex Output/RSEM/DEG/STAR.ADC.Nonrecur.Normal-AIS.DEG.Up.tex Output/RSEM/DEG/STAR.ADC.Nonrecur.Normal-Primary.DEG.Up.tex Output/RSEM/DEG/STAR.ADC.Nonrecur.AIS-Primary.DEG.Up.tex Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.Normal-AIS.DEG.Down.tex Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.Normal-Primary.DEG.Down.tex Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.AIS-Primary.DEG.Down.tex Output/RSEM/DEG/STAR.ADC.Nonrecur.Normal-AIS.DEG.Down.tex Output/RSEM/DEG/STAR.ADC.Nonrecur.Normal-Primary.DEG.Down.tex Output/RSEM/DEG/STAR.ADC.Nonrecur.AIS-Primary.DEG.Down.tex
step18: step18-preparation step18-SQC step18-SQC-DEG step18-SQC-Recur step18-SQC-Recur-DEG step18-SQC-Stage step18-SQC-Stage-DEG step18-ADC step18-ADC-DEG step18-ADC-Recur step18-ADC-Recur-DEG
.PHONY: step18-preparation step18-SQC step18-SQC-DEG step18-SQC-Recur step18-SQC-Recur-DEG step18-SQC-Stage step18-SQC-Stage-DEG step18-ADC step18-ADC-DEG step18-ADC-Recur step18-ADC-Recur-DEG step18

# Step 19 (DEG: further analysis)
Output/RSEM/Plots/Cluster Output/RSEM/Plots/Venn Output/RSEM/Plots/Volcano Output/RSEM/Plots/Violin Output/RSEM/Enrichment:
	mkdir -p $@

Output/RSEM/Plots/Cluster/%.SQC.pdf: Python/draw_clustermap_RSEM.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | Output/RSEM/Plots/Cluster build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Cluster/%.ADC.pdf: Python/draw_clustermap_RSEM.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | Output/RSEM/Plots/Cluster build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Cluster/%.ADC-SQC.pdf: Python/draw_clustermap_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | Output/RSEM/Plots/Cluster build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.ADC.Up.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.ADC.Normal-AAH.DEG.tsv Output/RSEM/DEG/%.ADC.Normal-AIS.DEG.tsv Output/RSEM/DEG/%.ADC.Normal-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-AAH" "Normal-AIS" "Normal-Primary" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.ADC.Down.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.ADC.Normal-AAH.DEG.tsv Output/RSEM/DEG/%.ADC.Normal-AIS.DEG.tsv Output/RSEM/DEG/%.ADC.Normal-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-AAH" "Normal-AIS" "Normal-Primary" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.ADC.Nonrecur.Up.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.ADC.Nonrecur.Normal-AIS.DEG.tsv Output/RSEM/DEG/%.ADC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.ADC.Nonrecur.AIS-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-AIS" "Normal-Primary" "AIS-Primary" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.ADC.Nonrecur.Down.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.ADC.Nonrecur.Normal-AIS.DEG.tsv Output/RSEM/DEG/%.ADC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.ADC.Nonrecur.AIS-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-AIS" "Normal-Primary" "AIS-Primary" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.Up.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Normal-Dysplasia.DEG.tsv Output/RSEM/DEG/%.SQC.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Normal-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-Dysplasia" "Normal-CIS" "Normal-Primary" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.Down.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Normal-Dysplasia.DEG.tsv Output/RSEM/DEG/%.SQC.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Normal-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-Dysplasia" "Normal-CIS" "Normal-Primary" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.Recur.Up.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.CIS-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-CIS" "Normal-Primary" "CIS-Primary" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.Recur.Down.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.CIS-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-CIS" "Normal-Primary" "CIS-Primary" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.Nonrecur.Up.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.CIS-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-CIS" "Normal-Primary" "CIS-Primary" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.Nonrecur.Down.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.CIS-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-CIS" "Normal-Primary" "CIS-Primary" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC-CIS.Recur-Nonrecur.Up.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Recur" "Non-recur" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC-CIS.Recur-Nonrecur.Down.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Recur" "Non-recur" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC-Primary.Recur-Nonrecur.Up.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Recur" "Non-recur" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC-Primary.Recur-Nonrecur.Down.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Recur" "Non-recur" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Volcano/%.volcano.pdf: Python/draw_volcano_plot.py Output/RSEM/DEG/%.DEG.tsv | Output/RSEM/Plots/Volcano build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Violin/%.ADC.violin.tar: Python/draw_violin_plots.py Output/RSEM/DEG/%.ADC.Normal-AAH.tsv Output/RSEM/DEG/%.ADC.Normal-AIS.tsv Output/RSEM/DEG/%.ADC.Normal-Primary.tsv Output/RSEM/DEG/%.ADC.AAH-AIS.tsv Output/RSEM/DEG/%.ADC.AAH-Primary.tsv Output/RSEM/DEG/%.ADC.AIS-Primary.tsv | Output/RSEM/Plots/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Violin/%.SQC.violin.tar: Python/draw_violin_plots.py Output/RSEM/DEG/%.SQC.Normal-Dysplasia.tsv Output/RSEM/DEG/%.SQC.Normal-CIS.tsv Output/RSEM/DEG/%.SQC.Normal-Primary.tsv Output/RSEM/DEG/%.SQC.Dysplasia-CIS.tsv Output/RSEM/DEG/%.SQC.Dysplasia-Primary.tsv Output/RSEM/DEG/%.SQC.CIS-Primary.tsv | Output/RSEM/Plots/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.Up.KEGG.tsv Output/RSEM/Enrichment/%.Up.KEGG.tex: Python/get_enrichment_enrichr.py Output/RSEM/DEG/%.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.Down.KEGG.tsv Output/RSEM/Enrichment/%.Down.KEGG.tex: Python/get_enrichment_enrichr.py Output/RSEM/DEG/%.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.Up.MSigDB.tsv Output/RSEM/Enrichment/%.Up.MSigDB.tex: Python/get_enrichment_enrichr.py Output/RSEM/DEG/%.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.Down.MSigDB.tsv Output/RSEM/Enrichment/%.Down.MSigDB.tex: Python/get_enrichment_enrichr.py Output/RSEM/DEG/%.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Recur.Up.KEGG.tsv Output/RSEM/Enrichment/%.SQC.Recur.Up.Subtract.KEGG.tex: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --up --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Recur.Down.KEGG.tsv Output/RSEM/Enrichment/%.SQC.Recur.Down.Subtract.KEGG.tex: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --down --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Recur.Up.MSigDB.tsv Output/RSEM/Enrichment/%.SQC.Recur.Up.Subtract.MSigDB.tex: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --up --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Recur.Down.MSigDB.tsv Output/RSEM/Enrichment/%.SQC.Recur.Down.Subtract.MSigDB.tex: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --down --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Nonrecur.Up.KEGG.tsv Output/RSEM/Enrichment/%.SQC.Nonrecur.Up.Subtract.KEGG.tex: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --up --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Nonrecur.Down.KEGG.tsv Output/RSEM/Enrichment/%.SQC.Nonrecur.Down.Subtract.KEGG.tex: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --down --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Nonrecur.Up.MSigDB.tsv Output/RSEM/Enrichment/%.SQC.Nonrecur.Up.Subtract.MSigDB.tex: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --up --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Nonrecur.Down.MSigDB.tsv Output/RSEM/Enrichment/%.SQC.Nonrecur.Down.Subtract.MSigDB.tex: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --down --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Recur.Up.KEGG.tsv Output/RSEM/Enrichment/%.SQC.Recur.Up.Subtract.KEGG.tex: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --up --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Recur.Down.KEGG.tsv Output/RSEM/Enrichment/%.SQC.Recur.Down.Subtract.KEGG.tex: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --down --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Recur.Up.MSigDB.tsv Output/RSEM/Enrichment/%.SQC.Recur.Up.Subtract.MSigDB.tex: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --up --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Recur.Down.MSigDB.tsv Output/RSEM/Enrichment/%.SQC.Recur.Down.Subtract.MSigDB.tex: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --down --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Nonrecur.Up.KEGG.tsv Output/RSEM/Enrichment/%.SQC.Nonrecur.Up.Subtract.KEGG.tex: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --up --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Nonrecur.Down.KEGG.tsv Output/RSEM/Enrichment/%.SQC.Nonrecur.Down.Subtract.KEGG.tex: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --down --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Nonrecur.Up.MSigDB.tsv Output/RSEM/Enrichment/%.SQC.Nonrecur.Up.Subtract.MSigDB.tex: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --up --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Nonrecur.Down.MSigDB.tsv Output/RSEM/Enrichment/%.SQC.Nonrecur.Down.Subtract.MSigDB.tex: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --down --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/SQC.%.pair.pdf: Python/draw_pair_plot.py Output/RSEM/DEG/SQC.%.tsv Output/RSEM/DEG/SQC.%.DEG.tsv | Output/RSEM/Plots build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^  $@) --SQC

step19-cluster: Output/RSEM/Plots/Cluster/Bowtie2.SQC.pdf Output/RSEM/Plots/Cluster/STAR.SQC.pdf Output/RSEM/Plots/Cluster/Bowtie2.ADC.pdf Output/RSEM/Plots/Cluster/STAR.ADC.pdf Output/RSEM/Plots/Cluster/Bowtie2.ADC-SQC.pdf Output/RSEM/Plots/Cluster/STAR.ADC-SQC.pdf
step19-venn-ADC: Output/RSEM/Plots/Venn/Bowtie2.ADC.Up.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.ADC.Down.venn.pdf Output/RSEM/Plots/Venn/STAR.ADC.Up.venn.pdf Output/RSEM/Plots/Venn/STAR.ADC.Down.venn.pdf
step19-venn-ADC-Recur: Output/RSEM/Plots/Venn/Bowtie2.ADC.Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.ADC.Nonrecur.Down.venn.pdf Output/RSEM/Plots/Venn/STAR.ADC.Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/STAR.ADC.Nonrecur.Down.venn.pdf
step19-venn-SQC: Output/RSEM/Plots/Venn/Bowtie2.SQC.Up.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC.Down.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.Up.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.Down.venn.pdf
step19-venn-SQC-Recur: Output/RSEM/Plots/Venn/Bowtie2.SQC.Recur.Up.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC.Recur.Down.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.Recur.Up.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.Recur.Down.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC.Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC.Nonrecur.Down.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.Nonrecur.Down.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC-CIS.Recur-Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC-CIS.Recur-Nonrecur.Down.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC-Primary.Recur-Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC-Primary.Recur-Nonrecur.Down.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC-CIS.Recur-Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC-CIS.Recur-Nonrecur.Down.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC-Primary.Recur-Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC-Primary.Recur-Nonrecur.Down.venn.pdf
step19-volcano-ADC: Output/RSEM/Plots/Volcano/Bowtie2.ADC.Normal-AAH.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.ADC.Normal-AIS.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.ADC.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.ADC.AAH-AIS.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.ADC.AAH-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.ADC.AIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.Normal-AAH.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.Normal-AIS.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.AAH-AIS.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.AAH-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.AIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.ADC.Nonrecur.Normal-AIS.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.ADC.Nonrecur.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.ADC.Nonrecur.AIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.Nonrecur.Normal-AIS.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.Nonrecur.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.Nonrecur.AIS-Primary.volcano.pdf
step19-volcano-SQC: Output/RSEM/Plots/Volcano/Bowtie2.SQC.Normal-Dysplasia.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Normal-CIS.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Dysplasia-CIS.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Dysplasia-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.CIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Normal-Dysplasia.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Normal-CIS.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Dysplasia-CIS.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Dysplasia-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.CIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Recur.Normal-CIS.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Recur.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Recur.Normal-CIS.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Recur.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Recur.Normal-CIS.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Recur.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Recur.CIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Recur.Normal-CIS.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Recur.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Recur.CIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Nonrecur.Normal-CIS.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Nonrecur.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Nonrecur.CIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Nonrecur.Normal-CIS.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Nonrecur.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Nonrecur.CIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.CIS.Nonrecur-Recur.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Primary.Nonrecur-Recur.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.CIS.Nonrecur-Recur.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Primary.Nonrecur-Recur.volcano.pdf
step19-violin-ADC: Output/RSEM/Plots/Violin/Bowtie2.ADC.violin.tar Output/RSEM/Plots/Violin/STAR.ADC.violin.tar
step19-violin-SQC: Output/RSEM/Plots/Violin/Bowtie2.SQC.violin.tar Output/RSEM/Plots/Violin/STAR.SQC.violin.tar
step19-KEGG-ADC: Output/RSEM/Enrichment/Bowtie2.ADC.Normal-AAH.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-AIS.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AAH-AIS.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AAH-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AIS-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-AAH.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-AIS.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AAH-AIS.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AAH-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AIS-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-AAH.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-AIS.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.AAH-AIS.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.AAH-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.AIS-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-AAH.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-AIS.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.AAH-AIS.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.AAH-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.AIS-Primary.Down.KEGG.tsv
step19-MSigDB-ADC: Output/RSEM/Enrichment/Bowtie2.ADC.Normal-AAH.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-AIS.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AAH-AIS.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AAH-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AIS-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-AAH.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-AIS.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AAH-AIS.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AAH-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AIS-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-AAH.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-AIS.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.AAH-AIS.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.AAH-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.AIS-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-AAH.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-AIS.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.AAH-AIS.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.AAH-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.AIS-Primary.Down.MSigDB.tsv
step19-KEGG-ADC-Recur: Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.Normal-AIS.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.AIS-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.Normal-AIS.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.AIS-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.Normal-AIS.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.AIS-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.Normal-AIS.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.AIS-Primary.Down.KEGG.tsv
step19-MSigDB-ADC-Recur: Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.Normal-AIS.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.AIS-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.Normal-AIS.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.AIS-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.Normal-AIS.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.AIS-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.Normal-AIS.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.AIS-Primary.Down.MSigDB.tsv
step19-KEGG-SQC: Output/RSEM/Enrichment/Bowtie2.SQC.Normal-Dysplasia.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-CIS.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Dysplasia-CIS.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Dysplasia-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.CIS-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-Dysplasia.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-CIS.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Dysplasia-CIS.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Dysplasia-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.CIS-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-Dysplasia.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-CIS.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Dysplasia-CIS.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Dysplasia-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.CIS-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-Dysplasia.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-CIS.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Dysplasia-CIS.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Dysplasia-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.CIS-Primary.Down.KEGG.tsv
step19-MSigDB-SQC: Output/RSEM/Enrichment/Bowtie2.SQC.Normal-Dysplasia.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-CIS.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Dysplasia-CIS.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Dysplasia-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.CIS-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-Dysplasia.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-CIS.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Dysplasia-CIS.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Dysplasia-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.CIS-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-Dysplasia.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-CIS.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Dysplasia-CIS.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Dysplasia-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.CIS-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-Dysplasia.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-CIS.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Dysplasia-CIS.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Dysplasia-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.CIS-Primary.Down.MSigDB.tsv
step19-KEGG-SQC-Recur: Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Recur.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Recur.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Nonrecur.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Nonrecur.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Recur.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Recur.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Nonrecur.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Nonrecur.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Recur.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Recur.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Nonrecur.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Nonrecur.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Recur.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Recur.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Nonrecur.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Nonrecur.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Recur.Normal-CIS.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Recur.Normal-CIS.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Recur.Normal-CIS.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Recur.Normal-CIS.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Nonrecur.Normal-CIS.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Nonrecur.Normal-CIS.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Nonrecur.Normal-CIS.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Nonrecur.Normal-CIS.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Recur.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Recur.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Recur.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Recur.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Nonrecur.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Nonrecur.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Nonrecur.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Nonrecur.Normal-Primary.Down.KEGG.tsv
step19-MSigDB-SQC-Recur: Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Recur.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Recur.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Nonrecur.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Nonrecur.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Recur.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Recur.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Nonrecur.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Nonrecur.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Recur.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Recur.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Nonrecur.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Nonrecur.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Recur.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Recur.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Nonrecur.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Nonrecur.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Recur.Normal-CIS.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Recur.Normal-CIS.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Recur.Normal-CIS.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Recur.Normal-CIS.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Nonrecur.Normal-CIS.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Nonrecur.Normal-CIS.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Nonrecur.Normal-CIS.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Nonrecur.Normal-CIS.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Recur.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Recur.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Recur.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Recur.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Nonrecur.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Nonrecur.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Nonrecur.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Nonrecur.Normal-Primary.Down.MSigDB.tsv
step19: step19-cluster step19-venn-ADC step19-venn-ADC-Recur step19-venn-SQC step19-venn-SQC-Recur step19-volcano-ADC step19-volcano-SQC step19-violin-ADC step19-violin-SQC step19-KEGG-ADC step19-MSigDB-ADC step19-KEGG-ADC-Recur step19-MSigDB-ADC-Recur step19-KEGG-SQC step19-MSigDB-SQC step19-KEGG-SQC-Recur step19-MSigDB-SQC-Recur
.PHONY: step19-cluster step19-venn-ADC step19-venn-ADC-Recur step19-venn-SQC step19-venn-SQC-Recur step19-volcano-ADC step19-volcano-SQC step19-violin-ADC step19-violin-SQC step19-KEGG-ADC step19-MSigDB-ADC step19-KEGG-ADC-Recur step19-MSigDB-ADC-Recur step19-KEGG-SQC step19-MSigDB-SQC step19-KEGG-SQC-Recur step19-MSigDB-SQC-Recur step19

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
Output/PyClone/BWA Output/PyClone/Bowtie2 Output/PyClone/BWA-Plots Output/PyClone/Bowtie2-Plots:
	mkdir -p $@

Output/PyClone/BWA/%.tsv Output/PyClone/BWA/%.Cellularity.txt: Python/make_input_pyclone.py Output/Mutect2/BWA/%.BWA.PASS.maf Output/Sequenza/BWA/%.tar.gz | Output/PyClone/BWA build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(@D)/$*) --cpus $(CPUS) 1> $(@D)/$*.preprocess.stdout 2> $(@D)/$*.preprocess.stderr

Output/PyClone/Bowtie2/%.tsv Output/PyClone/Bowtie2/%.Cellularity.txt: Python/make_input_pyclone.py Output/Mutect2/Bowtie2/%.Bowtie2.PASS.maf Output/Sequenza/Bowtie2/%.tar.gz | Output/PyClone/Bowtie2 build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(@D)/$*) --cpus $(CPUS) 1> $(@D)/$*.preprocess.stdout 2> $(@D)/$*.preprocess.stderr

# I know that this recipe contains unnecessary dependencies; but I cannot see the solution.
Output/PyClone/BWA/%/loci.tsv: Bash/run_pyclone.bash $(addprefix Output/PyClone/BWA/,$(addsuffix .tsv,$(WES_ID))) $(addprefix Output/PyClone/BWA/,$(addsuffix .Cellularity.txt,$(WES_ID))) | build.log
	$(DOCKER) /bin/bash $(addprefix /,$< Output/PyClone/BWA/$*.tsv Output/PyClone/BWA/$(PAIRED_PRIMARY).tsv Output/PyClone/BWA/$*.Cellularity.txt Output/PyClone/BWA/$(PAIRED_PRIMARY).Cellularity.txt $(@D)) 1> $(@D).stdout 2> $(@D).stderr || true

# I know that this recipe contains unnecessary dependencies; but I cannot see the solution.
Output/PyClone/Bowtie2/%/loci.tsv: Bash/run_pyclone.bash $(addprefix Output/PyClone/Bowtie2/,$(addsuffix .tsv,$(WES_ID))) $(addprefix Output/PyClone/Bowtie2/,$(addsuffix .Cellularity.txt,$(WES_ID))) | build.log
	$(DOCKER) /bin/bash $(addprefix /,$< Output/PyClone/Bowtie2/$*.tsv Output/PyClone/Bowtie2/$(PAIRED_PRIMARY).tsv Output/PyClone/Bowtie2/$*.Cellularity.txt Output/PyClone/Bowtie2/$(PAIRED_PRIMARY).Cellularity.txt $(@D)) 1> $(@D).stdout 2> $(@D).stderr || true

Output/PyClone/BWA-Plots/%.CCF.pdf: Python/draw_output_pyclone.py Output/PyClone/BWA/%/loci.tsv Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/PyClone/BWA-Plots
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --CCF 1> $@.stdout 2> $@.stderr

Output/PyClone/BWA-Plots/%.VAF.pdf: Python/draw_output_pyclone.py Output/PyClone/BWA/%/loci.tsv Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/PyClone/BWA-Plots
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --VAF 1> $@.stdout 2> $@.stderr

Output/PyClone/Bowtie2-Plots/%.CCF.pdf: Python/draw_output_pyclone.py Output/PyClone/Bowtie2/%/loci.tsv Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/PyClone/Bowtie2-Plots
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --CCF 1> $@.stdout 2> $@.stderr

Output/PyClone/Bowtie2-Plots/%.VAF.pdf: Python/draw_output_pyclone.py Output/PyClone/Bowtie2/%/loci.tsv Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/PyClone/Bowtie2-Plots
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --VAF 1> $@.stdout 2> $@.stderr

step21-preprocess: $(addprefix Output/PyClone/BWA/,$(addsuffix .tsv,$(WES_TUMOR_ID))) $(addprefix Output/PyClone/Bowtie2/,$(addsuffix .tsv,$(WES_TUMOR_ID)))
step21-pyclone: $(addprefix Output/PyClone/BWA/,$(addsuffix /loci.tsv,$(filter-out %P,$(filter-out %P1,$(WES_TUMOR_ID))))) $(addprefix Output/PyClone/Bowtie2/,$(addsuffix /loci.tsv,$(filter-out %P,$(filter-out %P1,$(WES_TUMOR_ID)))))
step21-plot: $(addprefix Output/PyClone/BWA-Plots/,$(addsuffix .CCF.pdf,$(filter-out %P,$(filter-out %P1,$(WES_TUMOR_ID))))) $(addprefix Output/PyClone/BWA-Plots/,$(addsuffix .VAF.pdf,$(filter-out %P,$(filter-out %P1,$(WES_TUMOR_ID))))) $(addprefix Output/PyClone/Bowtie2-Plots/,$(addsuffix .CCF.pdf,$(filter-out %P,$(filter-out %P1,$(WES_TUMOR_ID))))) $(addprefix Output/PyClone/Bowtie2-Plots/,$(addsuffix .VAF.pdf,$(filter-out %P,$(filter-out %P1,$(WES_TUMOR_ID)))))
step21: step21-preprocess step21-pyclone step21-plot
.PHONY: step21-preprocess step21-pyclone step21-plot step21

# Step 22 (Run CIBERSORTx)
Output/CIBERSORTx:
	mkdir -p $@

Output/CIBERSORTx/Bowtie2.SQC.tsv: Python/aggregate_RSEM_cibersort.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(WTS_ID))) $(CLINICAL_DATA_CSV) $(GENCODE_GTF) $(TREMBL_TSV) | Output/CIBERSORTx build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/Bowtie2.ADC.tsv: Python/aggregate_RSEM_cibersort.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(WTS_ID))) $(CLINICAL_DATA_CSV) $(GENCODE_GTF) $(TREMBL_TSV) | Output/CIBERSORTx build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/STAR.SQC.tsv: Python/aggregate_RSEM_cibersort.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(WTS_ID))) $(CLINICAL_DATA_CSV) $(GENCODE_GTF) $(TREMBL_TSV) | Output/CIBERSORTx build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/STAR.ADC.tsv: Python/aggregate_RSEM_cibersort.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(WTS_ID))) $(CLINICAL_DATA_CSV) $(GENCODE_GTF) $(TREMBL_TSV) | Output/CIBERSORTx build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/GSE131907_signature.tsv: Python/make_signature_cibersort.py $(LUNG_SINGLECELL_ANNOTATION_TXT) $(LUNG_SINGLECELL_NORMALIZED_MATRIX_TXT) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/%/CIBERSORTxGEP_NA_Weights.txt: Output/CIBERSORTx/%.tsv Output/CIBERSORTx/GSE131907_signature.tsv | build.log
	$(DOCKER) rm -rfv $(@D)
	docker run $(RUN_OPTS) --volume $(abspath $(dir $<)):/src/data --volume $(abspath $(@D)):/src/outdir cibersortx/hires $(CIBERSORTx_TOKEN) --mixture $(addprefix /,$(word 1,$^)) --sigmatrix $(addprefix /,$(word 2,$^)) --threads $(CPUS) --QN TRUE --heatmap TRUE --cluster TRUE 1> $(@D).stdout 2> $(@D).stderr || true

Output/CIBERSORTx/%.heatmap.pdf: Python/draw_heatmap_cibersort.py Output/CIBERSORTx/%/CIBERSORTxGEP_NA_Weights.txt | build
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/%.violin.tar: Python/draw_violin_plot_cibersort.py Output/CIBERSORTx/%/CIBERSORTxGEP_NA_Weights.txt | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/%.cluster.pdf: Python/draw_clustermap_cibersort.py Output/CIBERSORTx/%/CIBERSORTxGEP_NA_Weights.txt | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step22-RSEM: Output/CIBERSORTx/Bowtie2.SQC.tsv Output/CIBERSORTx/Bowtie2.ADC.tsv Output/CIBERSORTx/STAR.SQC.tsv Output/CIBERSORTx/STAR.ADC.tsv
step22-CIBERSORTx: Output/CIBERSORTx/GSE131907_signature.tsv Output/CIBERSORTx/Bowtie2.SQC/CIBERSORTxGEP_NA_Weights.txt Output/CIBERSORTx/Bowtie2.ADC/CIBERSORTxGEP_NA_Weights.txt Output/CIBERSORTx/STAR.SQC/CIBERSORTxGEP_NA_Weights.txt Output/CIBERSORTx/STAR.ADC/CIBERSORTxGEP_NA_Weights.txt
step22-heatmap: Output/CIBERSORTx/Bowtie2.SQC.heatmap.pdf Output/CIBERSORTx/Bowtie2.ADC.heatmap.pdf Output/CIBERSORTx/STAR.SQC.heatmap.pdf Output/CIBERSORTx/STAR.ADC.heatmap.pdf
step22-violin: Output/CIBERSORTx/Bowtie2.SQC.violin.tar Output/CIBERSORTx/Bowtie2.ADC.violin.tar Output/CIBERSORTx/STAR.SQC.violin.tar Output/CIBERSORTx/STAR.ADC.violin.tar
step22-cluster: Output/CIBERSORTx/Bowtie2.SQC.cluster.pdf Output/CIBERSORTx/STAR.SQC.cluster.pdf Output/CIBERSORTx/Bowtie2.ADC.cluster.pdf Output/CIBERSORTx/STAR.ADC.cluster.pdf
step22: step22-RSEM step22-CIBERSORTx step22-heatmap step22-violin step22-cluster
.PHONY: step22-RSEM step22-CIBERSORTx step22-heatmap step22-violin step22-cluster step22

# Step 23 (VAF plot)
Output/VAF Output/VAF-tSNE:
	mkdir -p $@

Output/VAF/BWA.SQC.tar: Python/draw_vaf_plots.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENSUS_CSV) | Output/VAF build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/VAF/BWA.ADC.tar: Python/draw_vaf_plots.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENSUS_CSV) | Output/VAF build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/VAF/Bowtie2.SQC.tar: Python/draw_vaf_plots.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENSUS_CSV) | Output/VAF build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/VAF/Bowtie2.ADC.tar: Python/draw_vaf_plots.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENSUS_CSV) | Output/VAF build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/VAF-tSNE/BWA.SQC.pdf: Python/draw_vaf_tsne.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/VAF-tSNE build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/VAF-tSNE/BWA.ADC.pdf: Python/draw_vaf_tsne.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENSUS_CSV) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt | Output/VAF-tSNE build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC

Output/VAF-tSNE/Bowtie2.SQC.pdf: Python/draw_vaf_tsne.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/VAF-tSNE build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/VAF-tSNE/Bowtie2.ADC.pdf: Python/draw_vaf_tsne.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/VAF-tSNE build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

step23-VAF: Output/VAF/BWA.SQC.tar Output/VAF/BWA.ADC.tar Output/VAF/Bowtie2.SQC.tar Output/VAF/Bowtie2.ADC.tar
step23-tSNE: Output/VAF-tSNE/BWA.SQC.pdf Output/VAF-tSNE/BWA.ADC.pdf Output/VAF-tSNE/Bowtie2.SQC.pdf Output/VAF-tSNE/Bowtie2.ADC.pdf
step23: step23-VAF
.PHONY: step23-VAF step23

# Step 24 (Revolver)
Output/Revolver/BWA Output/Revolver/Bowtie2:
	mkdir -p $@

Output/Revolver/BWA/SQC.tsv: Python/make_input_revolver.py $(addprefix Output/PyClone/BWA/,$(addsuffix /loci.tsv,$(filter-out %P,$(filter-out %P1,$(WES_TUMOR_ID))))) $(CLINICAL_DATA_CSV) Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/Revolver/BWA
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Revolver/BWA/ADC.tsv: Python/make_input_revolver.py $(addprefix Output/PyClone/BWA/,$(addsuffix /loci.tsv,$(filter-out %P,$(filter-out %P1,$(WES_TUMOR_ID))))) $(CLINICAL_DATA_CSV) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/Revolver/BWA
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Revolver/Bowtie2/SQC.tsv: Python/make_input_revolver.py $(addprefix Output/PyClone/Bowtie2/,$(addsuffix /loci.tsv,$(filter-out %P,$(filter-out %P1,$(WES_TUMOR_ID))))) $(CLINICAL_DATA_CSV) Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/Revolver/Bowtie2
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Revolver/Bowtie2/ADC.tsv: Python/make_input_revolver.py $(addprefix Output/PyClone/Bowtie2/,$(addsuffix /loci.tsv,$(filter-out %P,$(filter-out %P1,$(WES_TUMOR_ID))))) $(CLINICAL_DATA_CSV) Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENSUS_CSV) | build.log Output/Revolver/Bowtie2
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Revolver/BWA/%.cluster.pdf Output/Revolver/BWA/%.driver.pdf Output/Revolver/BWA/%.dendrogram.pdf Output/Revolver/BWA/%.det.pdf Output/Revolver/BWA/%.clonality.pdf: R/run_revolver.R Output/Revolver/BWA/%.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.tsv,$^)) --cluster $(addprefix /,$(@D)/$*.cluster.pdf) --driver $(addprefix /,$(@D)/$*.driver.pdf) --dendrogram $(addprefix /,$(@D)/$*.dendrogram.pdf) --det $(addprefix /,$(@D)/$*.det.pdf) --clonality $(addprefix /,$(@D)/$*.clonality.pdf) 1> $@.stdout 2> $@.stderr

Output/Revolver/Bowtie2/%.cluster.pdf Output/Revolver/Bowtie2/%.driver.pdf Output/Revolver/Bowtie2/%.dendrogram.pdf Output/Revolver/Bowtie2/%.det.pdf Output/Revolver/Bowtie2/%.clonality.pdf: R/run_revolver.R Output/Revolver/Bowtie2/%.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.tsv,$^)) --cluster $(addprefix /,$(@D)/$*.cluster.pdf) --driver $(addprefix /,$(@D)/$*.driver.pdf) --dendrogram $(addprefix /,$(@D)/$*.dendrogram.pdf) --det $(addprefix /,$(@D)/$*.det.pdf) --clonality $(addprefix /,$(@D)/$*.clonality.pdf) 1> $@.stdout 2> $@.stderr

step24-input: Output/Revolver/BWA/SQC.tsv Output/Revolver/BWA/ADC.tsv Output/Revolver/Bowtie2/SQC.tsv Output/Revolver/Bowtie2/ADC.tsv
step24-revolver: Output/Revolver/BWA/SQC.cluster.pdf Output/Revolver/Bowtie2/SQC.cluster.pdf# Output/Revolver/BWA/ADC.cluster.pdf Output/Revolver/Bowtie2/ADC.cluster.pdf
step24: step24-input step24-revolver
.PHONY: step24-input step24-revolver step24

# Step 25 (Arriba)
Output/Arriba/Reference Output/Arriba/Raw Output/Arriba/Filtered:
	mkdir -p $@

Output/Arriba/Reference/hg38.fa: | Output/Arriba/Reference build.log
	$(ARRIBA_DOCKER) /bin/bash /Tools/arriba_v2.1.0/download_references.sh "hg38+RefSeq" 1> $@.stdout 2> $@.stderr

Output/Arriba/Raw/%/fusions.tsv: $(GENCODE_GTF) $(REFERENCE_FASTA) Data/WTS/%_*.fastq.gz | Output/Arriba/Raw build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D)) && mkdir -p $(@D)
	$(ARRIBA_DOCKER) run_arriba.sh $(addprefix /,$(ARRIBA_STAR_DIR) $(filter %.gtf,$^) $(filter %.fasta,$^) $(ARRIBA_BLACKLIST_TSV) $(ARRIBA_KNOWN_FUSIONS_TSV) $(ARRIBA_PROTEIN_DOMAIN_GFF3)) $(CPUS) $(addprefix /,$(sort $(filter %.fastq.gz,$^))) 1> $(@D).stdout 2> $(@D).stderr

# I know that this recipe contains unnecessary dependencies; but I cannot see the solution.
# It is working even if there is no matched normal sample.
Output/Arriba/Filtered/%.tsv: Python/filter_arriba.py Output/Arriba/Raw/%/fusions.tsv $(addprefix Output/Arriba/Raw/,$(addsuffix /fusions.tsv,$(WTS_NORMAL_ID))) | build.log Output/Arriba/Filtered
	$(DOCKER) $(PYTHON) $(addprefix /,$< Output/Arriba/Raw/$*/fusions.tsv Output/Arriba/Raw/$(PAIRED_NORMAL)/fusions.tsv $@) 1> $@.stdout 2> $@.stderr

Output/Arriba/ADC.pdf: Python/aggregate_arriba.py $(addprefix Output/Arriba/Filtered/,$(addsuffix .tsv,$(WTS_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENSUS_CSV) | Output/Arriba build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Arriba/SQC.pdf: Python/aggregate_arriba.py $(addprefix Output/Arriba/Filtered/,$(addsuffix .tsv,$(WTS_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENSUS_CSV) | Output/Arriba build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

step25-Raw: $(addprefix Output/Arriba/Raw/,$(addsuffix /fusions.tsv,$(WTS_ID)))
step25-Filtered: $(addprefix Output/Arriba/Filtered/,$(addsuffix .tsv,$(WTS_TUMOR_ID)))
step25-Plot: Output/Arriba/ADC.pdf Output/Arriba/SQC.pdf
step25: step25-Raw step25-Filtered step25-Plot
.PHONY: step25-Raw step25-Filtered step25-Plot step25

# Step 26 (FusionCatcher)
Output/FusionCatcher/Reference Output/FusionCatcher/STAR Output/FusionCatcher/Bowtie2 Output/FusionCatcher/Plot:
	mkdir -p $@

Output/FusionCatcher/Reference/human_v102/genes.fa: Bash/download_human-db.bash | build.log Output/FusionCatcher/Reference
	$(DOCKER) /bin/bash $(addprefix /,$^ $(@D)) 1> $(@D).stdout 2> $(@D).stderr

Output/FusionCatcher/STAR/%/normal.tsv Output/FusionCatcher/STAR/%/tumor.tsv: Output/FusionCatcher/Reference/human_v102/genes.fa Data/WTS/%_*.fastq.gz | Output/FusionCatcher/STAR build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D)) && mkdir -p $(@D)
	if [[ -e "Data/WTS/$(PAIRED_NORMAL)_*.fastq.gz" ]]; then $(DOCKER) fusioncatcher.py --input=$(addprefix /,$(word 1,$(sort $(filter %.fastq.gz,$^)))),$(addprefix /,$(word 2,$(sort $(filter %.fastq.gz,$^)))) --output=$(addprefix /,$(@D)) --data=$(addprefix /,$(dir $(filter %.fa,$^))) --tmp=$(TMPFS) --threads=$(CPUS) --no-update-check --aligners=star --skip-conversion-grch37 --Xmx=$(MEMS)g 1> $(@D).stdout 2> $(@D).stderr; else $(DOCKER) fusioncatcher.py --input=$(addprefix /,$(word 1,$(sort $(filter %.fastq.gz,$^)))),$(addprefix /,$(word 2,$(sort $(filter %.fastq.gz,$^)))) --normal=$(addprefix /,$(word 1,$(sort $(wildcard Data/WTS/$(PAIRED_NORMAL)_*.fastq.gz)))),$(addprefix /,$(word 2,$(sort $(wildcard Data/WTS/$(PAIRED_NORMAL)_*.fastq.gz)))) --output=$(addprefix /,$(@D)) --data=$(addprefix /,$(dir $(filter %.fa,$^))) --threads=$(CPUS) --no-update-check --aligners=star --skip-conversion-grch37 --Xmx=$(MEMS)g 1> $(@D).stdout 2> $(@D).stderr; fi;
	cp $(@D)/1_*/final-list_candidate-fusion-genes.txt $(@D)/normal.tsv || touch $(@D)/normal.tsv
	cp $(@D)/2_*/final-list_candidate-fusion-genes.txt $(@D)/tumor.tsv

Output/FusionCatcher/Bowtie2/%/normal.tsv Output/FusionCatcher/Bowtie2/%/tumor.tsv: Output/FusionCatcher/Reference/human_v102/genes.fa Data/WTS/%_*.fastq.gz | Output/FusionCatcher/Bowtie2 build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D)) && mkdir -p $(@D)
	if [[ -e "Data/WTS/$(PAIRED_NORMAL)_*.fastq.gz" ]]; then $(DOCKER) fusioncatcher.py --input=$(addprefix /,$(word 1,$(sort $(filter %.fastq.gz,$^)))),$(addprefix /,$(word 2,$(sort $(filter %.fastq.gz,$^)))) --output=$(addprefix /,$(@D)) --data=$(addprefix /,$(dir $(filter %.fa,$^))) --tmp=$(TMPFS) --threads=$(CPUS) --no-update-check --aligners=bowtie2 --skip-conversion-grch37 --Xmx=$(MEMS)g 1> $(@D).stdout 2> $(@D).stderr; else $(DOCKER) fusioncatcher.py --input=$(addprefix /,$(word 1,$(sort $(filter %.fastq.gz,$^)))),$(addprefix /,$(word 2,$(sort $(filter %.fastq.gz,$^)))) --normal=$(addprefix /,$(word 1,$(sort $(wildcard Data/WTS/$(PAIRED_NORMAL)_*.fastq.gz)))),$(addprefix /,$(word 2,$(sort $(wildcard Data/WTS/$(PAIRED_NORMAL)_*.fastq.gz)))) --output=$(addprefix /,$(@D)) --data=$(addprefix /,$(dir $(filter %.fa,$^))) --threads=$(CPUS) --no-update-check --aligners=bowtie2 --skip-conversion-grch37 --Xmx=$(MEMS)g 1> $(@D).stdout 2> $(@D).stderr; fi;
	cp $(@D)/1_*/final-list_candidate-fusion-genes.txt $(@D)/normal.tsv || touch $(@D)/normal.tsv
	cp $(@D)/2_*/final-list_candidate-fusion-genes.txt $(@D)/tumor.tsv

Output/FusionCatcher/STAR/%.tsv: Python/filter_fusioncatcher.py Output/FusionCatcher/STAR/%/normal.tsv Output/FusionCatcher/STAR/%/tumor.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/FusionCatcher/Bowtie2/%.tsv: Python/filter_fusioncatcher.py Output/FusionCatcher/Bowtie2/%/normal.tsv Output/FusionCatcher/Bowtie2/%/tumor.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/FusionCatcher/Plot/ADC.STAR.pdf: Python/aggregate_fusioncatcher.py $(addprefix Output/FusionCatcher/STAR/,$(addsuffix .tsv,$(filter-out %N,$(ADC_WTS_SAMPLES)))) | Output/FusionCatcher/Plot build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS)

step26-prepar: Output/FusionCatcher/Reference/human_v102/genes.fa
step26-STAR: $(addprefix Output/FusionCatcher/STAR/,$(addsuffix /normal.tsv,$(WTS_TUMOR_ID))) $(addprefix Output/FusionCatcher/STAR/,$(addsuffix .tsv,$(WTS_TUMOR_ID)))
step26-Bowtie: $(addprefix Output/FusionCatcher/Bowtie2/,$(addsuffix /normal.tsv,$(WTS_TUMOR_ID))) $(addprefix Output/FusionCatcher/Bowtie2/,$(addsuffix .tsv,$(WTS_TUMOR_ID)))
step26-Plot:
step26: step26-prepare step26-STAR step26-Bowtie
.PHONY: step26-prepare step26-STAR step26-Bowtie step26

# Step 27 (STAR-Fusion)
Output/STAR-Fusion:
	mkdir -p $@

Output/STAR-Fusion/%/star-fusion.fusion_predictions.tsv Output/STAR-Fusion/%/star-fusion.fusion_predictions.abridged.tsv: Data/WTS/%_*.fastq.gz | Output/STAR-Fusion build.log
	$(DOCKER) rm -rf $(addsuffix /,$(@D)) && mkdir -p $(@D)
	$(DOCKER) STAR-Fusion --genome_lib_dir $(addprefix /,$(CTAT_RESOURCE_LIB_DIR)) --left_fq $(addprefix /,$(word 1,$(sort $(filter %.fastq.gz,$^)))) --right_fq $(addprefix /,$(word 2,$(sort $(filter %.fastq.gz,$^)))) --CPU $(CPUS) --output_dir $(addprefix /,$(@D)) 1> $(@D).stdout 2> $(@D).stderr

step27: $(addprefix Output/STAR-Fusion/,$(addsuffix /star-fusion.fusion_predictions.tsv,$(WTS_TUMOR_ID)))
.PHONY: step27

# Step 28 (MuSiC)
Output/MuSiC:
	mkdir -p $@

Output/MuSiC/Bowtie2.SQC.RSEM.tsv: Python/aggregate_RSEM_deconvolution.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(WTS_ID))) $(CLINICAL_DATA_CSV) $(GENCODE_GTF) $(TREMBL_TSV) | Output/MuSiC build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/MuSiC/Bowtie2.ADC.RSEM.tsv: Python/aggregate_RSEM_deconvolution.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(WTS_ID))) $(CLINICAL_DATA_CSV) $(GENCODE_GTF) $(TREMBL_TSV) | Output/MuSiC build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/MuSiC/STAR.SQC.RSEM.tsv: Python/aggregate_RSEM_deconvolution.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(WTS_ID))) $(CLINICAL_DATA_CSV) $(GENCODE_GTF) $(TREMBL_TSV) | Output/MuSiC build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/MuSiC/STAR.ADC.RSEM.tsv: Python/aggregate_RSEM_deconvolution.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(WTS_ID))) $(CLINICAL_DATA_CSV) $(GENCODE_GTF) $(TREMBL_TSV) | Output/MuSiC build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/MuSiC/%.tsv: R/run_MuSiC.R Output/MuSiC/%.RSEM.tsv $(LUNG_SINGLECELL_ANNOTATION_TXT) $(LUNG_SINGLECELL_RAW_MATRIX_RDS) | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.tsv,$^)) --annotation $(addprefix /,$(filter %.txt,$^)) --reference $(addprefix /,$(filter %.rds,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/MuSiC/%.heatmap.pdf: Python/draw_heatmap_MuSiC.py Output/MuSiC/%.tsv | build
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/MuSiC/%.violin.tar: Python/draw_violin_plot_MuSiC.py Output/MuSiC/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MuSiC/%.cluster.pdf: Python/draw_clustermap_MuSiC.py Output/MuSiC/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step28-prepare: Output/MuSiC/Bowtie2.SQC.RSEM.tsv Output/MuSiC/Bowtie2.ADC.RSEM.tsv Output/MuSiC/STAR.SQC.RSEM.tsv Output/MuSiC/STAR.ADC.RSEM.tsv
step28-MuSiC: Output/MuSiC/Bowtie2.SQC.tsv Output/MuSiC/Bowtie2.ADC.tsv Output/MuSiC/STAR.SQC.tsv Output/MuSiC/STAR.ADC.tsv
step28-heatmap: Output/MuSiC/Bowtie2.SQC.heatmap.pdf Output/MuSiC/Bowtie2.ADC.heatmap.pdf Output/MuSiC/STAR.SQC.heatmap.pdf Output/MuSiC/STAR.ADC.heatmap.pdf
step28-violin: Output/MuSiC/Bowtie2.SQC.violin.tar Output/MuSiC/Bowtie2.ADC.violin.tar Output/MuSiC/STAR.SQC.violin.tar Output/MuSiC/STAR.ADC.violin.tar
step28-cluster: Output/MuSiC/Bowtie2.SQC.cluster.pdf Output/MuSiC/Bowtie2.ADC.cluster.pdf Output/MuSiC/STAR.SQC.cluster.pdf Output/MuSiC/STAR.ADC.cluster.pdf
step28: step28-prepare step28-MuSiC step28-heatmap step28-violin step28-cluster
.PHONY: step28-prepare step28-MuSiC step28-heatmap step28-violin step28-cluster step28

# Step 29 (BisqueRNA)
Output/BisqueRNA:
	mkdir -p $@

Output/BisqueRNA/%.tsv: R/run_BisqueRNA.R Output/MuSiC/%.RSEM.tsv $(LUNG_SINGLECELL_ANNOTATION_TXT) $(LUNG_SINGLECELL_RAW_MATRIX_RDS) | Output/BisqueRNA build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.tsv,$^)) --annotation $(addprefix /,$(filter %.txt,$^)) --reference $(addprefix /,$(filter %.rds,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/BisqueRNA/%.heatmap.pdf: Python/draw_heatmap_BisqueRNA.py Output/BisqueRNA/%.tsv | build
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/BisqueRNA/%.violin.tar: Python/draw_violin_plot_BisqueRNA.py Output/BisqueRNA/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/BisqueRNA/%.cluster.pdf: Python/draw_clustermap_BisqueRNA.py Output/BisqueRNA/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step29-BisqueRNA: Output/BisqueRNA/Bowtie2.SQC.tsv Output/BisqueRNA/Bowtie2.ADC.tsv Output/BisqueRNA/STAR.SQC.tsv Output/BisqueRNA/STAR.ADC.tsv
step29-heatmap: Output/BisqueRNA/Bowtie2.SQC.heatmap.pdf Output/BisqueRNA/Bowtie2.ADC.heatmap.pdf Output/BisqueRNA/STAR.SQC.heatmap.pdf Output/BisqueRNA/STAR.ADC.heatmap.pdf
step29-violin: Output/BisqueRNA/Bowtie2.SQC.violin.tar Output/BisqueRNA/Bowtie2.ADC.violin.tar Output/BisqueRNA/STAR.SQC.violin.tar Output/BisqueRNA/STAR.ADC.violin.tar
step29-cluster: Output/BisqueRNA/Bowtie2.SQC.cluster.pdf Output/BisqueRNA/Bowtie2.ADC.cluster.pdf Output/BisqueRNA/STAR.SQC.cluster.pdf Output/BisqueRNA/STAR.ADC.cluster.pdf
step29: step29-BisqueRNA step29-heatmap step29-violin step29-cluster
.PHONY: step29-BisqueRNA step29-heatmap step29-violin step29-cluster step29

# Step 30 (SCDC)
Output/SCDC:
	mkdir -p $@

Output/SCDC/%.tsv: R/run_SCDC.R Output/MuSiC/%.RSEM.tsv $(LUNG_SINGLECELL_ANNOTATION_TXT) $(LUNG_SINGLECELL_RAW_MATRIX_RDS) | Output/SCDC build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.tsv,$^)) --annotation $(addprefix /,$(filter %.txt,$^)) --reference $(addprefix /,$(filter %.rds,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/SCDC/%.heatmap.pdf: Python/draw_heatmap_MuSiC.py Output/SCDC/%.tsv | build
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/SCDC/%.violin.tar: Python/draw_violin_plot_MuSiC.py Output/SCDC/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/SCDC/%.cluster.pdf: Python/draw_clustermap_MuSiC.py Output/SCDC/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step30-SCDC: Output/SCDC/Bowtie2.SQC.tsv Output/SCDC/Bowtie2.ADC.tsv Output/SCDC/STAR.SQC.tsv Output/SCDC/STAR.ADC.tsv
step30-heatmap: Output/SCDC/Bowtie2.SQC.heatmap.pdf Output/SCDC/Bowtie2.ADC.heatmap.pdf Output/SCDC/STAR.SQC.heatmap.pdf Output/SCDC/STAR.ADC.heatmap.pdf
step30-violin: Output/SCDC/Bowtie2.SQC.violin.tar Output/SCDC/Bowtie2.ADC.violin.tar Output/SCDC/STAR.SQC.violin.tar Output/SCDC/STAR.ADC.violin.tar
step30-cluster: Output/SCDC/Bowtie2.SQC.cluster.pdf Output/SCDC/Bowtie2.ADC.cluster.pdf Output/SCDC/STAR.SQC.cluster.pdf Output/SCDC/STAR.ADC.cluster.pdf
step30: step30-SCDC step30-heatmap step30-violin step30-cluster
.PHONY: step30-SCDC step30-heatmap step30-violin step30-cluster step30

# Step 31 (TIDE)
Output/TIDE:
	mkdir -p $@

Output/TIDE/Bowtie2.SQC.tsv: Python/aggregate_RSEM_TIDE.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(WTS_ID))) $(CLINICAL_DATA_CSV) $(GENCODE_GTF) $(TREMBL_TSV) | Output/TIDE build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/TIDE/Bowtie2.ADC.tsv: Python/aggregate_RSEM_TIDE.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(WTS_ID))) $(CLINICAL_DATA_CSV) $(GENCODE_GTF) $(TREMBL_TSV) | Output/TIDE build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/TIDE/STAR.SQC.tsv: Python/aggregate_RSEM_TIDE.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(WTS_ID))) $(CLINICAL_DATA_CSV) $(GENCODE_GTF) $(TREMBL_TSV) | Output/TIDE build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/TIDE/STAR.ADC.tsv: Python/aggregate_RSEM_TIDE.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(WTS_ID))) $(CLINICAL_DATA_CSV) $(GENCODE_GTF) $(TREMBL_TSV) | Output/TIDE build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/TIDE/%.TIDE.tsv: Bash/run_TIDE.bash Output/TIDE/%.tsv | build.log
	$(DOCKER) /bin/bash $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step31-prepare: Output/TIDE/Bowtie2.SQC.tsv Output/TIDE/Bowtie2.ADC.tsv Output/TIDE/STAR.SQC.tsv Output/TIDE/STAR.ADC.tsv
step31-tide: Output/TIDE/Bowtie2.SQC.TIDE.tsv Output/TIDE/Bowtie2.ADC.TIDE.tsv Output/TIDE/STAR.SQC.TIDE.tsv Output/TIDE/STAR.ADC.TIDE.tsv
step31: step31-prepare step31-tide
.PHONY: step31-prepare step31-tide step31

# Step ## (Gistic2.0)
# http://portals.broadinstitute.org/cgi-bin/cancer/publications/pub_paper.cgi?mode=view&paper_id=216&p=t
# ftp://ftp.broadinstitute.org/pub/GISTIC2.0/INSTALL.txt

# Step ## (QAPA)
# https://github.com/morrislab/qapa

# Step ## (Roar)
# http://bioconductor.org/packages/release/bioc/html/roar.html

# Step ## (MesKit)
# Step ## (Spruce)
# STep ## (CRUST)
# Step ## (WTS - samtools depth)
# Step ## (Picard - CollectMultipleMetrics)
# Step ## (PureCN)
# Step ## (Histogram of file size)
