# Variables
DATE := $(shell date "+%y%m%d")
RANDOM := $(shell bash -c 'echo $$RANDOM')
IMAGE_NAME = lungcancer:latest
SHELL := /bin/bash
CPUS = 10
MEMS = 50
PARALLEL = 5
PWD := $(shell pwd)
TOOLS =
DOCKER = docker run --rm $(ENV_OPTS) $(VOLUME_OPTS) $(RUN_OPTS) $(IMAGE_NAME)
WORKDIR_DOCKER = docker run --rm $(ENV_OPTS) $(VOLUME_OPTS) $(RUN_OPTS) --workdir $(addprefix /,$(@D)) $(IMAGE_NAME)
TOOLS_DOCKER = docker run --rm $(ENV_OPTS) $(VOLUME_OPTS) $(RUN_OPTS) --workdir "/Tools" $(IMAGE_NAME)
PYTHON = python3 -B
R = Rscript --vanilla
BASH = /bin/bash
CHROMOSOMES = chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
TAXONOMY_TYPES = phylum class order family genus species
PAIRED_NORMAL = $(shell echo $* | grep --perl-regexp --only-match "^(cn)?\d+")N
PAIRED_PRIMARY = $(shell echo $* | grep --perl-regexp --only-match "^(cn)?\d+")P$(shell echo $* | grep --perl-regexp --only-match "\d$$" > /dev/null && echo "1")

# Options
VOLUME_OPTS = --volume $(abspath Output):/Output --volume $(abspath Data):/Data:ro --volume $(abspath Python):/Python:ro --volume $(abspath R):/R:ro --volume $(abspath Bash):/Bash:ro --volume $(abspath misc):/misc:ro --tmpfs $(TMPFS)
RUN_OPTS = --tty --cpus="$(CPUS)" --memory="$(MEMS)G"
ENV_OPTS = --env MLM_LICENSE_FILE=/misc/license.lic --env NUMEXPR_MAX_THREADS=$(CPUS)
JAVA_OPTS = -XX:+UseParallelGC -XX:ParallelGCThreads=1 -Xmx$(MEMS)g
GATK = gatk --java-options "$(JAVA_OPTS)"
PICARD = java $(JAVA_OPTS) -jar /Tools/picard.jar
VCF2MAF = perl /Tools/vcf2maf-1.6.19/vcf2maf.pl
MUTENRICHER = $(PYTHON) /Tools/MutEnricher-1.3.1/mutEnricher.py
TMPFS = /tmpfs
CIBERSORTx_TOKEN = --username "jwlee230@unist.ac.kr" --token "23aabc1d937bb0ae7fe0d4bf50c04b0c"
PYCLONE = /root/anaconda2/envs/pyclone/bin/PyClone
HELMSMAN = /root/anaconda2/envs/helmsman/bin/python3 /Tools/helmsman-1.5.3/helmsman.py
PURECN_EXTDATA_DIR = /usr/local/lib/R/site-library/PureCN/extdata

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
REFERENCE_BLACKLIST = Data/UCSC_hg38_blacklist.txt
REFERENCE_BAND = Data/cytoBand.txt
LIBRARY_BED = Data/S04380110_Regions.bed
LIBRARY_BED_GZ = $(LIBRARY_BED).gz
GENCODE_GTF = Data/GENCODE/gencode.v38.annotation.gtf
TREMBL_TSV = Data/Homo_sapiens.GRCh38.104.uniprot.tsv.gz
GENE_CENCUS_ONE_CSV = Data/CGC_tier1.csv
GENE_CENCUS_BOTH_CSV = Data/CGC_both.csv
LUNG_GSE131907_ANNOTATION_TXT = Data/Lung_SingleCell/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt
LUNG_GSE131907_RAW_MATRIX_RDS = Data/Lung_SingleCell/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.rds
LUNG_GSE131907_NORMALIZED_MATRIX_RDS = Data/Lung_SingleCell/GSE131907/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds
LUNG_GSE131907_NORMALIZED_MATRIX_TXT = Data/Lung_SingleCell/GSE131907/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt
LUNG_GSE162498_11_TUMORS_RDS = Data/Lung_SingleCell/GSE162498/GSE162498_NSCLC_CD3_11tumors.Rds
LUNG_GSE179994_ANNOTATION_TXT = Data/Lung_SingleCell/GSE179994/GSE179994_all.scTCR.tsv
LUNG_GSE179994_MATRIX_RDS = Data/Lung_SingleCell/GSE179994/GSE179994_all.Tcell.rawCounts.rds
KNOWN_GENE_GTF = Data/hg38.gtf.gz
ARRIBA_STAR_DIR = Tools/arriba_v2.1.0/STAR_index_hg38_RefSeq_hg38
ARRIBA_DATABASE_DIR = Tools/arriba_v2.1.0/database
ARRIBA_BLACKLIST_TSV = $(ARRIBA_DATABASE_DIR)/blacklist_hg38_GRCh38_v2.1.0.tsv.gz
ARRIBA_KNOWN_FUSIONS_TSV = $(ARRIBA_DATABASE_DIR)/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz
ARRIBA_PROTEIN_DOMAIN_GFF3 = $(ARRIBA_DATABASE_DIR)/protein_domains_hg38_GRCh38_v2.1.0.gff3
CTAT_RESOURCE_LIB_DIR = Data/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir
MAPPABILITY_BW = Data/k24.Bismap.MultiTrackMappability.bw
UCSC_REFFLAT_TXT = Data/refFlat.txt

# Sample IDs
PATIENT_ID := $(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+"))
WES_ID := $(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+[PCADMN]\d?"))
WES_REAL_ID := $(filter-out 33P1,$(WES_ID))
WES_NORMAL_ID := $(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+N"))
WES_TUMOR_ID := $(sort $(shell ls -1 Data/RawData/*.fastq.gz | cut --delimiter / --fields 3 | grep --perl-regexp --only-match "^(cn)?\d+[PCADM]\d?"))
WES_REAL_TUMOR_ID := $(filter-out 33P1,$(WES_TUMOR_ID))
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

Docker/Tools/Gistic.tar.gz: | Docker/Tools
	wget "ftp://ftp.broadinstitute.org/distribution/GISTIC2.0/GISTIC_2_0_23.tar.gz" -O $@
TOOLS += Docker/Tools/Gistic.tar.gz

Docker/Tools/MiXCR.zip: | Docker/Tools
	wget "https://github.com/milaboratory/mixcr/releases/download/v3.0.13/mixcr-3.0.13.zip" -O $@
TOOLS += Docker/Tools/MiXCR.zip

Docker/Tools/RTCR.zip: | Docker/Tools
	wget "https://github.com/uubram/RTCR/zipball/master" -O $@
TOOLS += Docker/Tools/RTCR.zip

Docker/Tools/helmsman.tar.gz: | Docker/Tools
	wget "https://github.com/carjed/helmsman/archive/refs/tags/1.5.3.tar.gz" -O $@
TOOLS += Docker/Tools/helmsman.tar.gz

# Docker
build.log: Docker/Dockerfile Docker/requirements.txt $(TOOLS) | log Output
	rm -f $@
	docker images | grep $(IMAGE_NAME) && docker rmi $(IMAGE_NAME) || true
	docker build --rm --tag $(IMAGE_NAME) $(dir $(word 1,$^)) | tee $@

build: build.log
.PHONY: build

interactive: build.log
	docker run --rm $(ENV_OPTS) $(VOLUME_OPTS) $(RUN_OPTS) --interactive $(IMAGE_NAME) $(BASH) || true
.PHONY: interactive

delete: build.log
	docker rmi $(IMAGE_NAME)
	rm -rfv build.log
.PHONY: delete

stop:
	docker rm $(CONTAINER_NAME)
.PHONY: stop

# Slurm
run: | log Output build.log
	echo -e "#!/bin/bash\nmake -j $(PARALLEL) -C $(PWD) latest" > tmp.sh
	sbatch --chdir=$(abspath .) --cpus-per-task=$(shell echo "$(CPUS) * $(PARALLEL)" | bc) --error="$(abspath log)/%x.e%A" --job-name="Lung_$(DATE)" --mail-type=ALL --mail-user="230@fumire.moe" --mem=$(shell echo "$(MEMS) * $(PARALLEL)" | bc)G --output="$(abspath log)/%x.o%A" --export=ALL tmp.sh
.PHONY: run

# Actual
latest: step07 step08 step09 step10 step11 step12 step13 step14 step15 step16 step17 step18 step19 step20 step21 step22 step23 step24 step25 step26 step27 step28 step29 step30 step31 step32 step33 step34 step35 step36 step37 step38 step39 step40 step42 step43 step44 step45 step46 step47 step48 step49 step50 step51 step52 step53 step54 step55 step56 step57 step58 step59 step60 step61 step62 step63
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
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --threshold 2 1> $@.stdout 2> $@.stderr

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
	$(DOCKER) $(BASH) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/WES/%.Sort.MarkDuplicate.bai: Output/WES/%.Sort.MarkDuplicate.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step06: $(addprefix Output/WES/,$(addsuffix .BWA.Sort.dm-tmp.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.bai,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.dm-tmp.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.bam,$(WES_ID))) $(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.bai,$(WES_ID)))
.PHONY: step06

# Step 07 (WES: Sort.MarkDuplicate.bam >> Sort.MarkDuplicate.BQSR.bam)
#Output/WES/%.Sort.MarkDuplicate.BQSR.table: $(REFERENCE_FASTA) $(REFERENCE_KNOWN_VCF) Output/WES/%.Sort.MarkDuplicate.bam Output/WES/%.Sort.MarkDuplicate.bai | build.log
#	$(DOCKER) $(GATK) BaseRecalibrator --input $(addprefix /,$(filter %.bam,$^)) $(foreach site,$(addprefix /,$(filter %.vcf,$^)),--known-sites $(site)) --output $(addprefix /,$@) --reference $(addprefix /,$(filter %.fasta,$^)) 1> $@.stdout 2> $@.stderr

#Output/WES/%.Sort.MarkDuplicate.BQSR.bam: Output/WES/%.Sort.MarkDuplicate.bam Output/WES/%.Sort.MarkDuplicate.BQSR.table $(REFERENCE_FASTA) | build.log
#	$(DOCKER) $(GATK) ApplyBQSR --bqsr-recal-file $(addprefix /,$(filter %.table,$^)) --input $(addprefix /,$(filter %.bam,$^)) --output $(addprefix /,$@) --reference $(addprefix /,$(filter %.fasta,$^)) 1> $@.stdout 2> $@.stderr

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

Output/Sequenza/Bowtie2/%/sample_alternative_solutions.txt Output/Sequenza/Bowtie2/%/sample_segments.txt: R/run_sequenza.R Output/Sequenza/Bowtie2/%.seqz.gz Output/Sequenza/Bowtie2/%.seqz.gz.tbi | build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D)) && mkdir -p $(@D)
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --sampleid "sample" --input $(addprefix /,$(filter %.seqz.gz,$^)) --outdir $(addprefix /,$(@D)) --parallel $(CPUS) 1> $(@D).stdout 2> $(@D).stderr

Output/Sequenza/BWA-sequenza-SQC.pdf: Python/aggregate_sequenza.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-sequenza-ADC.pdf: Python/aggregate_sequenza.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-sequenza-SQC.pdf: Python/aggregate_sequenza.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-sequenza-ADC.pdf: Python/aggregate_sequenza.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-SQC.pdf Output/Sequenza/BWA-genome-SQC.png: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-ADC.pdf Output/Sequenza/BWA-genome-ADC.png: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-genome-SQC.pdf Output/Sequenza/Bowtie2-genome-SQC.png: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-genome-ADC.pdf Output/Sequenza/Bowtie2-genome-ADC.png: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-SQC.patient.pdf Output/Sequenza/BWA-genome-SQC.patient.png: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-ADC.patient.pdf Output/Sequenza/BWA-genome-ADC.patient.png: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-genome-SQC.patient.pdf Output/Sequenza/Bowtie2-genome-SQC.patient.png: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-genome-ADC.patient.pdf Output/Sequenza/Bowtie2-genome-ADC.patient.png: Python/aggregate_sequenza_genome.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-simple-SQC.pdf: Python/aggregate_sequenza_simple.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-simple-ADC.pdf: Python/aggregate_sequenza_simple.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-simple-SQC.pdf: Python/aggregate_sequenza_simple.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-simple-ADC.pdf: Python/aggregate_sequenza_simple.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequnza/BWA-correlation-SQC.pdf: Python/aggregate_sequenza_correlation.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC

Output/Sequenza/BWA-genome-SQC.Recurrence.pdf Output/Sequenza/BWA-genome-SQC.Recurrence.png: Python/aggregate_sequenza_genome_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-ADC.Recurrence.pdf Output/Sequenza/BWA-genome-ADC.Recurrence.png: Python/aggregate_sequenza_genome_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-genome-SQC.Recurrence.pdf Output/Sequenza/Bowtie2-genome-SQC.Recurrence.png: Python/aggregate_sequenza_genome_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-genome-ADC.Recurrence.pdf Output/Sequenza/Bowtie2-genome-ADC.Recurrence.png: Python/aggregate_sequenza_genome_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-SQC.Recurrence.patient.pdf Output/Sequenza/BWA-genome-SQC.Recurrence.patient.png: Python/aggregate_sequenza_genome_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-ADC.Recurrence.patient.pdf Output/Sequenza/BWA-genome-ADC.Recurrence.patient.png: Python/aggregate_sequenza_genome_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-genome-SQC.Recurrence.patient.pdf Output/Sequenza/Bowtie2-genome-SQC.Recurrence.patient.png: Python/aggregate_sequenza_genome_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-genome-ADC.Recurrence.patient.pdf Output/Sequenza/Bowtie2-genome-ADC.Recurrence.patient.png: Python/aggregate_sequenza_genome_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-SQC.Smoking.pdf Output/Sequenza/BWA-genome-SQC.Smoking.png: Python/aggregate_sequenza_genome_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-ADC.Smoking.pdf Output/Sequenza/BWA-genome-ADC.Smoking.png: Python/aggregate_sequenza_genome_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-genome-SQC.Smoking.pdf Output/Sequenza/Bowtie2-genome-SQC.Smoking.png: Python/aggregate_sequenza_genome_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-genome-ADC.Smoking.pdf Output/Sequenza/Bowtie2-genome-ADC.Smoking.png: Python/aggregate_sequenza_genome_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-SQC.Smoking.patient.pdf Output/Sequenza/BWA-genome-SQC.Smoking.patient.png: Python/aggregate_sequenza_genome_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-genome-ADC.Smoking.patient.pdf Output/Sequenza/BWA-genome-ADC.Smoking.patient.png: Python/aggregate_sequenza_genome_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-genome-SQC.Smoking.patient.pdf Output/Sequenza/Bowtie2-genome-SQC.Smoking.patient.png: Python/aggregate_sequenza_genome_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-genome-ADC.Smoking.patient.pdf Output/Sequenza/Bowtie2-genome-ADC.Smoking.patient.png: Python/aggregate_sequenza_genome_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-simple-SQC.Recurrence.pdf: Python/aggregate_sequenza_simple_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-simple-ADC.Recurrence.pdf: Python/aggregate_sequenza_simple_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-simple-SQC.Recurrence.pdf: Python/aggregate_sequenza_simple_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-simple-ADC.Recurrence.pdf: Python/aggregate_sequenza_simple_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-simple-SQC.Smoking.pdf: Python/aggregate_sequenza_simple_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-simple-ADC.Smoking.pdf: Python/aggregate_sequenza_simple_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-simple-SQC.Smoking.pdf: Python/aggregate_sequenza_simple_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-simple-ADC.Smoking.pdf: Python/aggregate_sequenza_simple_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-violin-SQC.pdf: Python/aggregate_sequenza_violin.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-violin-ADC.pdf: Python/aggregate_sequenza_violin.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-violin-SQC.pdf: Python/aggregate_sequenza_violin.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-violin-ADC.pdf: Python/aggregate_sequenza_violin.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-violin-SQC.Recurrence.pdf: Python/aggregate_sequenza_violin_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-violin-ADC.Recurrence.pdf: Python/aggregate_sequenza_violin_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-violin-SQC.Recurrence.pdf: Python/aggregate_sequenza_violin_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-violin-ADC.Recurrence.pdf: Python/aggregate_sequenza_violin_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-violin-SQC.Smoking.pdf: Python/aggregate_sequenza_violin_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-violin-ADC.Smoking.pdf: Python/aggregate_sequenza_violin_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-violin-SQC.Smoking.pdf: Python/aggregate_sequenza_violin_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-violin-ADC.Smoking.pdf: Python/aggregate_sequenza_violin_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-gene-SQC.tar: Python/aggregate_sequenza_violin_gene.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-gene-ADC.tar: Python/aggregate_sequenza_violin_gene.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-gene-SQC.tar: Python/aggregate_sequenza_violin_gene.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-gene-ADC.tar: Python/aggregate_sequenza_violin_gene.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-gene-SQC.Recurrence.tar: Python/aggregate_sequenza_violin_gene_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-gene-ADC.Recurrence.tar: Python/aggregate_sequenza_violin_gene_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-gene-SQC.Recurrence.tar: Python/aggregate_sequenza_violin_gene_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-gene-ADC.Recurrence.tar: Python/aggregate_sequenza_violin_gene_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-gene-SQC.Smoking.tar: Python/aggregate_sequenza_violin_gene_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/Sequenza/BWA-gene-ADC.Smoking.tar: Python/aggregate_sequenza_violin_gene_clinical.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-gene-SQC.Smoking.tar: Python/aggregate_sequenza_violin_gene_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/Sequenza/Bowtie2-gene-ADC.Smoking.tar: Python/aggregate_sequenza_violin_gene_clinical.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

step13-Sequenza: $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID)))
step13-default: Output/Sequenza/BWA-sequenza-SQC.pdf Output/Sequenza/BWA-sequenza-ADC.pdf Output/Sequenza/Bowtie2-sequenza-SQC.pdf Output/Sequenza/Bowtie2-sequenza-ADC.pdf
step13-genome: Output/Sequenza/BWA-genome-SQC.pdf Output/Sequenza/BWA-genome-ADC.pdf Output/Sequenza/Bowtie2-genome-SQC.pdf Output/Sequenza/Bowtie2-genome-ADC.pdf Output/Sequenza/BWA-genome-SQC.patient.pdf Output/Sequenza/BWA-genome-ADC.patient.pdf Output/Sequenza/Bowtie2-genome-SQC.patient.pdf Output/Sequenza/Bowtie2-genome-ADC.patient.pdf
step13-simple: Output/Sequenza/BWA-simple-SQC.pdf Output/Sequenza/BWA-simple-ADC.pdf Output/Sequenza/Bowtie2-simple-SQC.pdf Output/Sequenza/Bowtie2-simple-ADC.pdf
step13-genome-Clinical: Output/Sequenza/BWA-genome-SQC.Recurrence.pdf Output/Sequenza/BWA-genome-ADC.Recurrence.pdf Output/Sequenza/Bowtie2-genome-SQC.Recurrence.pdf Output/Sequenza/Bowtie2-genome-ADC.Recurrence.pdf Output/Sequenza/BWA-genome-SQC.Smoking.pdf Output/Sequenza/BWA-genome-ADC.Smoking.pdf Output/Sequenza/Bowtie2-genome-SQC.Smoking.pdf Output/Sequenza/Bowtie2-genome-ADC.Smoking.pdf Output/Sequenza/BWA-genome-SQC.Recurrence.patient.pdf Output/Sequenza/BWA-genome-ADC.Recurrence.patient.pdf Output/Sequenza/Bowtie2-genome-SQC.Recurrence.patient.pdf Output/Sequenza/Bowtie2-genome-ADC.Recurrence.patient.pdf Output/Sequenza/BWA-genome-SQC.Smoking.patient.pdf Output/Sequenza/BWA-genome-ADC.Smoking.patient.pdf Output/Sequenza/Bowtie2-genome-SQC.Smoking.patient.pdf Output/Sequenza/Bowtie2-genome-ADC.Smoking.patient.pdf
step13-simple-Clinical: Output/Sequenza/BWA-simple-SQC.Recurrence.pdf Output/Sequenza/BWA-simple-ADC.Recurrence.pdf Output/Sequenza/Bowtie2-simple-SQC.Recurrence.pdf Output/Sequenza/Bowtie2-simple-ADC.Recurrence.pdf Output/Sequenza/BWA-simple-SQC.Smoking.pdf Output/Sequenza/BWA-simple-ADC.Smoking.pdf Output/Sequenza/Bowtie2-simple-SQC.Smoking.pdf Output/Sequenza/Bowtie2-simple-ADC.Smoking.pdf
step13-violin: Output/Sequenza/BWA-violin-SQC.pdf Output/Sequenza/BWA-violin-ADC.pdf Output/Sequenza/Bowtie2-violin-SQC.pdf Output/Sequenza/Bowtie2-violin-ADC.pdf
step13-violin-Clinical: Output/Sequenza/BWA-violin-SQC.Recurrence.pdf Output/Sequenza/BWA-violin-ADC.Recurrence.pdf Output/Sequenza/Bowtie2-violin-SQC.Recurrence.pdf Output/Sequenza/Bowtie2-violin-ADC.Recurrence.pdf Output/Sequenza/BWA-violin-SQC.Smoking.pdf Output/Sequenza/BWA-violin-ADC.Smoking.pdf Output/Sequenza/Bowtie2-violin-SQC.Smoking.pdf Output/Sequenza/Bowtie2-violin-ADC.Smoking.pdf
step13-gene: Output/Sequenza/BWA-gene-SQC.tar Output/Sequenza/BWA-gene-ADC.tar Output/Sequenza/Bowtie2-gene-SQC.tar Output/Sequenza/Bowtie2-gene-ADC.tar
step313-gene-clinical: Output/Sequenza/BWA-gene-SQC.Recurrence.tar Output/Sequenza/BWA-gene-ADC.Recurrence.tar Output/Sequenza/Bowtie2-gene-SQC.Recurrence.tar Output/Sequenza/Bowtie2-gene-ADC.Recurrence.tar Output/Sequenza/BWA-gene-SQC.Smoking.tar Output/Sequenza/BWA-gene-ADC.Smoking.tar Output/Sequenza/Bowtie2-gene-SQC.Smoking.tar Output/Sequenza/Bowtie2-gene-ADC.Smoking.tar
step13: step13-Sequenza step13-default step13-genome step13-simple step13-genome-Clinical step13-simple-Clinical step13-violin step13-violin-Clinical step13-gene step313-gene-clinical
.PHONY: step13-Sequenza step13-default step13-genome step13-simple step13-genome-Clinical step13-simple-Clinical step13-violin step13-violin-Clinical step13-gene step313-gene-clinical step13

# Step 14 (Mutect2)
Output/Mutect2/BWA Output/Mutect2/Bowtie2 Output/MutationSharedProportion:
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
	$(DOCKER) $(BASH) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2/%.Bowtie2.PASS.vcf: Bash/select_pass.bash Output/Mutect2/Bowtie2/%.Bowtie2.filter.vcf | build.log
	$(DOCKER) $(BASH) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA/%.BWA.PASS.maf: $(REFERENCE_FASTA) Output/Mutect2/BWA/%.BWA.PASS.vcf | build.log
	$(DOCKER) $(VCF2MAF) --vep-path "/Tools/ensembl-vep-release-103" --vep-data "/Data" --vep-forks $(CPUS) --ncbi-build "GRCh38" --input-vcf $(addprefix /,$(filter %.vcf,$^)) --output $(addprefix /,$@) --tumor-id "$*.BWA.sam" --normal-id "$(PAIRED_NORMAL).BWA.sam" --ref-fasta $(addprefix /,$(filter %.fasta,$^)) --tmp-dir $(TMPFS) 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2/%.Bowtie2.PASS.maf: $(REFERENCE_FASTA) Output/Mutect2/Bowtie2/%.Bowtie2.PASS.vcf | build.log
	$(DOCKER) $(VCF2MAF) --vep-path "/Tools/ensembl-vep-release-103" --vep-data "/Data" --vep-forks $(CPUS) --ncbi-build "GRCh38" --input-vcf $(addprefix /,$(filter %.vcf,$^)) --output $(addprefix /,$@) --tumor-id "$*.Bowtie2.sam" --normal-id "$(PAIRED_NORMAL).Bowtie2.sam" --ref-fasta $(addprefix /,$(filter %.fasta,$^)) --tmp-dir $(TMPFS) 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-SQC.histogram.pdf: Python/aggregate_mutect_histogram.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-ADC.histogram.pdf: Python/aggregate_mutect_histogram.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-SQC.histogram.pdf: Python/aggregate_mutect_histogram.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-ADC.histogram.pdf: Python/aggregate_mutect_histogram.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-SQC.heatmap.pdf: Python/aggregate_mutect_heatmap.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-ADC.heatmap.pdf: Python/aggregate_mutect_heatmap.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-SQC.heatmap.pdf: Python/aggregate_mutect_heatmap.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-ADC.heatmap.pdf: Python/aggregate_mutect_heatmap.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-SQC.clustermap.pdf: Python/aggregate_mutect_cluster.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-ADC.clustermap.pdf: Python/aggregate_mutect_cluster.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-SQC.clustermap.pdf: Python/aggregate_mutect_cluster.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-ADC.clustermap.pdf: Python/aggregate_mutect_cluster.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-SQC.venn.tar: Python/aggregate_mutect_venn.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-ADC.venn.tar: Python/aggregate_mutect_venn.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-SQC.venn.tar: Python/aggregate_mutect_venn.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-ADC.venn.tar: Python/aggregate_mutect_venn.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/MutationSharedProportion/BWA.tsv: Python/calculate_MutationSharedProportion.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/MutationSharedProportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationSharedProportion/Bowtie2.tsv: Python/calculate_MutationSharedProportion.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/MutationSharedProportion
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step14-BWA: $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.vcf,$(WES_TUMOR_ID))) $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.filter.vcf,$(WES_TUMOR_ID))) $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.vcf,$(WES_TUMOR_ID))) $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID)))
step14-Bowtie2: $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.vcf,$(WES_TUMOR_ID))) $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.filter.vcf,$(WES_TUMOR_ID))) $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.vcf,$(WES_TUMOR_ID))) $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID)))
step14-histogram: Output/Mutect2/BWA-SQC.histogram.pdf Output/Mutect2/BWA-ADC.histogram.pdf Output/Mutect2/Bowtie2-SQC.histogram.pdf Output/Mutect2/Bowtie2-ADC.histogram.pdf
step14-heatmap: Output/Mutect2/BWA-SQC.heatmap.pdf Output/Mutect2/BWA-ADC.heatmap.pdf Output/Mutect2/Bowtie2-SQC.heatmap.pdf Output/Mutect2/Bowtie2-ADC.heatmap.pdf
step14-clustermap: Output/Mutect2/BWA-SQC.clustermap.pdf Output/Mutect2/BWA-ADC.clustermap.pdf Output/Mutect2/Bowtie2-SQC.clustermap.pdf Output/Mutect2/Bowtie2-ADC.clustermap.pdf
step14-venn: Output/Mutect2/BWA-SQC.venn.tar Output/Mutect2/BWA-ADC.venn.tar Output/Mutect2/Bowtie2-SQC.venn.tar Output/Mutect2/Bowtie2-ADC.venn.tar
step14-Sharing: Output/MutationSharedProportion/BWA.tsv Output/MutationSharedProportion/Bowtie2.tsv
step14: step14-BWA step14-Bowtie2 step14-histogram step14-heatmap step14-clustermap step14-venn step14-Sharing
.PHONY: step14-BWA step14-Bowtie2 step14-histogram step14-heatmap step14-clustermap step14-venn step14-Sharing step14

# Step 15 (HaplotypeCaller)
Output/HaplotypeCaller/BWA Output/HaplotypeCaller/Bowtie2:
	mkdir -p $@

Output/HaplotypeCaller/BWA/%.g.vcf.gz: $(REFERENCE_FASTA) Output/WES/%.BWA.Sort.MarkDuplicate.BQSR.bam $(REFERENCE_DBSNP_VCF) | Output/HaplotypeCaller/BWA build.log
	$(DOCKER) $(GATK) HaplotypeCaller --reference $(addprefix /,$(filter %.fasta,$^)) --input $(addprefix /,$(filter %.bam,$^)) --output $(addprefix /,$@) --emit-ref-confidence GVCF --dbsnp $(addprefix /,$(filter %.vcf,$^)) --native-pair-hmm-threads $(CPUS) --native-pair-hmm-use-double-precision true 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/Bowtie2/%.g.vcf.gz: $(REFERENCE_FASTA) Output/WES/%.Bowtie2.Sort.MarkDuplicate.BQSR.bam $(REFERENCE_DBSNP_VCF) | Output/HaplotypeCaller/Bowtie2 build.log
	$(DOCKER) $(GATK) HaplotypeCaller --reference $(addprefix /,$(filter %.fasta,$^)) --input $(addprefix /,$(filter %.bam,$^)) --output $(addprefix /,$@) --emit-ref-confidence GVCF --dbsnp $(addprefix /,$(filter %.vcf,$^)) --native-pair-hmm-threads $(CPUS) --native-pair-hmm-use-double-precision true 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/BWA/%.sample_map.tsv: Python/make_sample_map_haplotypecaller.py $(addprefix Output/HaplotypeCaller/BWA/,$(addsuffix .g.vcf.gz,$(WES_ID))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --sample "$*" 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/Bowtie2/%.sample_map.tsv: Python/make_sample_map_haplotypecaller.py $(addprefix Output/HaplotypeCaller/Bowtie2/,$(addsuffix .g.vcf.gz,$(WES_ID))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --sample "$*" 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/BWA/%.genomicsdb/callset.json: $(REFERENCE_FASTA) Output/HaplotypeCaller/BWA/%.sample_map.tsv | build.log
	$(DOCKER) $(GATK) GenomicsDBImport --reference $(addprefix /,$(filter %.fasta,$^)) --genomicsdb-workspace-path $(addprefix /,$(@D)) --sample-name-map $(addprefix /,$(filter %.tsv,$^)) --reader-threads $(CPUS) --max-num-intervals-to-import-in-parallel $(CPUS) --overwrite-existing-genomicsdb-workspace true --merge-input-intervals true $(foreach chr,$(CHROMOSOMES),--intervals $(chr)) 1> $(@D).stdout 2> $(@D).stderr
	$(DOCKER) chmod --recursive 777 $(addprefix /,$(@D))

Output/HaplotypeCaller/Bowtie2/%.genomicsdb/callset.json: $(REFERENCE_FASTA) Output/HaplotypeCaller/Bowtie2/%.sample_map.tsv | build.log
	$(DOCKER) $(GATK) GenomicsDBImport --reference $(addprefix /,$(filter %.fasta,$^)) --genomicsdb-workspace-path $(addprefix /,$(@D)) --sample-name-map $(addprefix /,$(filter %.tsv,$^)) --reader-threads $(CPUS) --max-num-intervals-to-import-in-parallel $(CPUS) --overwrite-existing-genomicsdb-workspace true --merge-input-intervals true $(foreach chr,$(CHROMOSOMES),--intervals $(chr)) 1> $(@D).stdout 2> $(@D).stderr
	$(DOCKER) chmod --recursive 777 $(addprefix /,$(@D))

Output/HaplotypeCaller/BWA/%.vcf.gz: $(REFERENCE_FASTA) Output/HaplotypeCaller/BWA/%.genomicsdb/callset.json $(REFERENCE_DBSNP_VCF) | build.log
	$(DOCKER) $(GATK) GenotypeGVCFs --reference $(addprefix /,$(filter %.fasta,$^)) --variant "gendb://$(addprefix /,$(dir $(filter %.json,$^)))" --output $(addprefix /,$@) --dbsnp $(addprefix /,$(filter %.vcf,$^)) 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/Bowtie2/%.vcf.gz: $(REFERENCE_FASTA) Output/HaplotypeCaller/Bowtie2/%.genomicsdb/callset.json $(REFERENCE_DBSNP_VCF) | build.log
	$(DOCKER) $(GATK) GenotypeGVCFs --reference $(addprefix /,$(filter %.fasta,$^)) --variant "gendb://$(addprefix /,$(dir $(filter %.json,$^)))" --output $(addprefix /,$@) --dbsnp $(addprefix /,$(filter %.vcf,$^)) 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/BWA/%.VQSR.recal Output/HaplotypeCaller/BWA/%.VQSR.tranches: $(REFERENCE_FASTA) $(REFERENCE_HAPMAP_VCF) $(REFERENCE_OMNI_VCF) $(REFERENCE_1000G_VCF) $(REFERENCE_DBSNP_VCF) Output/HaplotypeCaller/BWA/%.vcf.gz | build.log
	$(DOCKER) $(GATK) VariantRecalibrator --reference $(addprefix /,$(filter %.fasta,$^)) --output $(addprefix /,$(basename $@).recal) --tranches-file $(addprefix /,$(basename $@).tranches) --variant $(addprefix /,$(lastword $^)) --use-annotation MQ --use-annotation MQRankSum --use-annotation ReadPosRankSum --use-annotation FS --use-annotation SOR --mode BOTH --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $(addprefix /,$(REFERENCE_HAPMAP_VCF)) --resource:omni,known=false,training=true,truth=false,prior=12.0 $(addprefix /,$(REFERENCE_OMNI_VCF)) --resource:1000G,known=false,training=true,truth=false,prior=10.0 $(addprefix /,$(REFERENCE_1000G_VCF)) --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $(addprefix /,$(REFERENCE_DBSNP_VCF)) 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/Bowtie2/%.VQSR.recal Output/HaplotypeCaller/Bowtie2/%.VQSR.tranches: $(REFERENCE_FASTA) $(REFERENCE_HAPMAP_VCF) $(REFERENCE_OMNI_VCF) $(REFERENCE_1000G_VCF) $(REFERENCE_DBSNP_VCF) Output/HaplotypeCaller/Bowtie2/%.vcf.gz | build.log
	$(DOCKER) $(GATK) VariantRecalibrator --reference $(addprefix /,$(filter %.fasta,$^)) --output $(addprefix /,$(basename $@).recal) --tranches-file $(addprefix /,$(basename $@).tranches) --variant $(addprefix /,$(lastword $^)) --use-annotation MQ --use-annotation MQRankSum --use-annotation ReadPosRankSum --use-annotation FS --use-annotation SOR --mode BOTH --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $(addprefix /,$(REFERENCE_HAPMAP_VCF)) --resource:omni,known=false,training=true,truth=false,prior=12.0 $(addprefix /,$(REFERENCE_OMNI_VCF)) --resource:1000G,known=false,training=true,truth=false,prior=10.0 $(addprefix /,$(REFERENCE_1000G_VCF)) --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $(addprefix /,$(REFERENCE_DBSNP_VCF)) 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/BWA/%.VQSR.vcf: $(REFERENCE_FASTA) Output/HaplotypeCaller/BWA/%.vcf.gz Output/HaplotypeCaller/BWA/%.VQSR.recal Output/HaplotypeCaller/BWA/%.VQSR.tranches | build.log
	$(DOCKER) $(GATK) ApplyVQSR --reference $(addprefix /,$(filter %.fasta,$^)) --variant $(addprefix /,$(filter %.vcf.gz,$^)) --recal-file $(addprefix /,$(filter %.recal,$^)) --tranches-file $(addprefix /,$(filter %.tranches,$^)) --output $(addprefix /,$@) --truth-sensitivity-filter-level 99.0 --mode BOTH 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/Bowtie2/%.VQSR.vcf: $(REFERENCE_FASTA) Output/HaplotypeCaller/Bowtie2/%.vcf.gz Output/HaplotypeCaller/Bowtie2/%.VQSR.recal Output/HaplotypeCaller/Bowtie2/%.VQSR.tranches | build.log
	$(DOCKER) $(GATK) ApplyVQSR --reference $(addprefix /,$(filter %.fasta,$^)) --variant $(addprefix /,$(filter %.vcf.gz,$^)) --recal-file $(addprefix /,$(filter %.recal,$^)) --tranches-file $(addprefix /,$(filter %.tranches,$^)) --output $(addprefix /,$@) --truth-sensitivity-filter-level 99.0 --mode BOTH 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/BWA/%.VQSR.PASS.vcf: Bash/select_pass.bash Output/HaplotypeCaller/BWA/%.VQSR.vcf | build.log
	$(DOCKER) $(BASH) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/Bowtie2/%.VQSR.PASS.vcf: Bash/select_pass.bash Output/HaplotypeCaller/Bowtie2/%.VQSR.vcf | build.log
	$(DOCKER) $(BASH) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/BWA/%.VQSR.PASS.maf: $(REFERENCE_FASTA) Output/HaplotypeCaller/BWA/%.VQSR.PASS.vcf | build.log
	$(DOCKER) $(VCF2MAF) --vep-path "/Tools/ensembl-vep-release-103" --vep-data "/Data" --vep-forks $(CPUS) --ncbi-build "GRCh38" --input-vcf $(addprefix /,$(filter %.vcf,$^)) --output $(addprefix /,$@) --tumor-id "$*" --normal-id "$(PAIRED_NORMAL)" --ref-fasta $(addprefix /,$(filter %.fasta,$^)) --tmp-dir $(TMPFS) 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/Bowtie2/%.VQSR.PASS.maf: $(REFERENCE_FASTA) Output/HaplotypeCaller/Bowtie2/%.VQSR.PASS.vcf | build.log
	$(DOCKER) $(VCF2MAF) --vep-path "/Tools/ensembl-vep-release-103" --vep-data "/Data" --vep-forks $(CPUS) --ncbi-build "GRCh38" --input-vcf $(addprefix /,$(filter %.vcf,$^)) --output $(addprefix /,$@) --tumor-id "$*" --normal-id "$(PAIRED_NORMAL)" --ref-fasta $(addprefix /,$(filter %.fasta,$^)) --tmp-dir $(TMPFS) 1> $@.stdout 2> $@.stderr

Output/HaplotypeCaller/BWA-SQC.pdf: Python/aggregate_haplotypecaller.py $(addprefix Output/HaplotypeCaller/BWA/,$(addsuffix .VQSR.PASS.maf,$(WES_TUMOR_ID))) $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --include --patient

step15-BWA: $(addprefix Output/HaplotypeCaller/BWA/,$(addsuffix .g.vcf.gz,$(WES_ID))) $(addprefix Output/HaplotypeCaller/BWA/,$(addsuffix .sample_map.tsv,$(WES_TUMOR_ID))) $(addprefix Output/HaplotypeCaller/BWA/,$(addsuffix .genomicsdb/callset.json,$(WES_TUMOR_ID))) $(addprefix Output/HaplotypeCaller/BWA/,$(addsuffix .vcf.gz,$(WES_TUMOR_ID))) $(addprefix Output/HaplotypeCaller/BWA/,$(addsuffix .VQSR.recal,$(WES_TUMOR_ID))) $(addprefix Output/HaplotypeCaller/BWA/,$(addsuffix .VQSR.vcf,$(WES_TUMOR_ID))) $(addprefix Output/HaplotypeCaller/BWA/,$(addsuffix .VQSR.PASS.vcf,$(WES_TUMOR_ID))) $(addprefix Output/HaplotypeCaller/BWA/,$(addsuffix .VQSR.PASS.maf,$(WES_TUMOR_ID)))
step15-Bowtie2: $(addprefix Output/HaplotypeCaller/Bowtie2/,$(addsuffix .g.vcf.gz,$(WES_ID))) $(addprefix Output/HaplotypeCaller/Bowtie2/,$(addsuffix .sample_map.tsv,$(WES_TUMOR_ID))) $(addprefix Output/HaplotypeCaller/Bowtie2/,$(addsuffix .genomicsdb/callset.json,$(WES_TUMOR_ID))) $(addprefix Output/HaplotypeCaller/Bowtie2/,$(addsuffix .vcf.gz,$(WES_TUMOR_ID))) $(addprefix Output/HaplotypeCaller/Bowtie2/,$(addsuffix .VQSR.recal,$(WES_TUMOR_ID))) $(addprefix Output/HaplotypeCaller/Bowtie2/,$(addsuffix .VQSR.vcf,$(WES_TUMOR_ID))) $(addprefix Output/HaplotypeCaller/Bowtie2/,$(addsuffix .VQSR.PASS.vcf,$(WES_TUMOR_ID))) $(addprefix Output/HaplotypeCaller/Bowtie2/,$(addsuffix .VQSR.PASS.maf,$(WES_TUMOR_ID)))
step15: step15-BWA step15-Bowtie2
.PHONY: step15-BWA step15-Bowtie2 step15

# Step 16 (MutEnricher)
Output/MutEnricher Output/MutEnricher/Reference Output/MutEnricher/BWA Output/MutEnricher/Bowtie2:
	mkdir -p $@

Output/MutEnricher/Reference/annovar_downdb.log Output/MutEnricher/Reference/hg38_refGene.txt Output/MutEnricher/Reference/hg38_refGeneMrna.fa Output/MutEnricher/Reference/hg38_refGeneVersion.txt: | Output/MutEnricher/Reference build.log
	$(DOCKER) annotate_variation.pl --verbose --buildver hg38 --downdb --webfrom annovar refGene $(addprefix /,$(@D)) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/BWA/%.hg38_multianno.vcf: Output/Mutect2/BWA/%.BWA.PASS.vcf Output/MutEnricher/Reference/hg38_refGeneMrna.fa | build.log Output/MutEnricher/BWA
	$(DOCKER) table_annovar.pl $(addprefix /,$(filter %.vcf,$^)) $(addprefix /,$(dir $(filter %.fa,$^))) --verbose --buildver "hg38" --remove --vcfinput --protocol "refGene" --operation g --nastring "." --outfile $(addprefix /,$(@D)/$*) --thread $(CPUS) --maxgenethread $(CPUS) --tempdir $(TMPFS) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/BWA/%.vcf.gz: Bash/prepare_VCF_annotation.bash Output/MutEnricher/BWA/%.hg38_multianno.vcf | build.log
	$(DOCKER) $(BASH) $(addprefix /,$^ $@)

Output/MutEnricher/BWA/%.vcf.gz.tbi: Output/MutEnricher/BWA/%.vcf.gz | build.log
	$(DOCKER) tabix --preset "vcf" --force $(addprefix /,$^) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/BWA/SQC.tsv: Python/make_list_MutEnricher.py $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz,$(WES_TUMOR_ID))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz.tbi,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^) $(filter %.vcf.gz,$^) $(filter %.csv,$^) $@) --SQC 1> $@.stdout 2> $@.stderr

Output/MutEnricher/BWA/ADC.tsv: Python/make_list_MutEnricher.py $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz,$(WES_TUMOR_ID))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz.tbi,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^) $(filter %.vcf.gz,$^) $(filter %.csv,$^) $@) --ADC 1> $@.stdout 2> $@.stderr

Output/MutEnricher/Bowtie2/%.hg38_multianno.vcf: Output/Mutect2/Bowtie2/%.Bowtie2.PASS.vcf Output/MutEnricher/Reference/hg38_refGeneMrna.fa | build.log Output/MutEnricher/Bowtie2
	$(DOCKER) table_annovar.pl $(addprefix /,$(filter %.vcf,$^)) $(addprefix /,$(dir $(filter %.fa,$^))) --verbose --buildver "hg38" --remove --vcfinput --protocol "refGene" --operation g --nastring "." --outfile $(addprefix /,$(@D)/$*) --thread $(CPUS) --maxgenethread $(CPUS) --tempdir $(TMPFS) 1> $@.stdout 2> $@.stderr

Output/MutEnricher/Bowtie2/%.vcf.gz: Bash/prepare_VCF_annotation.bash Output/MutEnricher/Bowtie2/%.hg38_multianno.vcf | build.log Output/MutEnricher/Bowtie2
	$(DOCKER) $(BASH) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

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

Output/Mutect2/BWA-SQC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --include 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-SQC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --include 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-ADC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --include 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-ADC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --include 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-SQC.Both.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_BOTH_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --include 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-SQC.Both.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_BOTH_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --include 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-ADC.Both.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_BOTH_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --include 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-ADC.Both.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_BOTH_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --include 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-SQC.non-CGC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_BOTH_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --exclude 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-SQC.non-CGC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_BOTH_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --exclude 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-ADC.non-CGC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_BOTH_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --exclude 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-ADC.non-CGC.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_BOTH_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --exclude 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-SQC.stage.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --include --stage 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-SQC.stage.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --include --stage 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-ADC.stage.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --include --stage 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-ADC.stage.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --include --stage 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-SQC.gene.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --include --gene 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-SQC.gene.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --include --gene 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-ADC.gene.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --include --gene 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-ADC.gene.pdf: Python/aggregate_mutect.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --include --gene 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-SQC.Recur.pdf: Python/aggregate_mutect_2.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-ADC.Recur.pdf: Python/aggregate_mutect_2.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-SQC.Recur.pdf: Python/aggregate_mutect_2.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-ADC.Recur.pdf: Python/aggregate_mutect_2.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-SQC.Smoking.pdf: Python/aggregate_mutect_2.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Mutect2/BWA-ADC.Smoking.pdf: Python/aggregate_mutect_2.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-SQC.Smoking.pdf: Python/aggregate_mutect_2.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Mutect2/Bowtie2-ADC.Smoking.pdf: Python/aggregate_mutect_2.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

step16-preprocess: Output/MutEnricher/gene_covariates.txt $(addprefix Output/MutEnricher/BWA/,$(addsuffix .hg38_multianno.vcf,$(WES_TUMOR_ID))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz,$(WES_TUMOR_ID))) $(addprefix Output/MutEnricher/BWA/,$(addsuffix .vcf.gz.tbi,$(WES_TUMOR_ID))) Output/MutEnricher/BWA/SQC.tsv Output/MutEnricher/BWA/ADC.tsv $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .hg38_multianno.vcf,$(WES_TUMOR_ID))) $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .vcf.gz,$(WES_TUMOR_ID))) $(addprefix Output/MutEnricher/Bowtie2/,$(addsuffix .vcf.gz.tbi,$(WES_TUMOR_ID))) Output/MutEnricher/Bowtie2/SQC.tsv Output/MutEnricher/Bowtie2/ADC.tsv
step16-MutEnricher: Output/MutEnricher/BWA/SQC/MutEnricher.log Output/MutEnricher/BWA/ADC/MutEnricher.log Output/MutEnricher/Bowtie2/ADC/MutEnricher.log Output/MutEnricher/Bowtie2/SQC/MutEnricher.log
step16-tier1: Output/Mutect2/BWA-SQC.pdf Output/Mutect2/BWA-ADC.pdf Output/Mutect2/Bowtie2-SQC.pdf Output/Mutect2/Bowtie2-ADC.pdf
step16-both: Output/Mutect2/BWA-SQC.Both.pdf Output/Mutect2/BWA-ADC.Both.pdf Output/Mutect2/Bowtie2-SQC.Both.pdf Output/Mutect2/Bowtie2-ADC.Both.pdf
step16-non: Output/Mutect2/BWA-SQC.non-CGC.pdf Output/Mutect2/BWA-ADC.non-CGC.pdf Output/Mutect2/Bowtie2-SQC.non-CGC.pdf Output/Mutect2/Bowtie2-ADC.non-CGC.pdf
step16-stage: Output/Mutect2/BWA-SQC.stage.pdf Output/Mutect2/Bowtie2-SQC.stage.pdf Output/Mutect2/BWA-ADC.stage.pdf Output/Mutect2/Bowtie2-ADC.stage.pdf
step16-gene: Output/Mutect2/BWA-SQC.gene.pdf Output/Mutect2/Bowtie2-SQC.gene.pdf Output/Mutect2/BWA-ADC.stage.pdf Output/Mutect2/Bowtie2-ADC.gene.pdf
step16-recur: Output/Mutect2/BWA-SQC.Recur.pdf Output/Mutect2/BWA-ADC.Recur.pdf Output/Mutect2/Bowtie2-SQC.Recur.pdf Output/Mutect2/Bowtie2-ADC.Recur.pdf
step16-smoking: Output/Mutect2/BWA-SQC.Smoking.pdf Output/Mutect2/BWA-ADC.Smoking.pdf Output/Mutect2/Bowtie2-SQC.Smoking.pdf Output/Mutect2/Bowtie2-ADC.Smoking.pdf
step16: step16-preprocess step16-MutEnricher step16-tier1 step16-both step16-non step16-stage step16-gene step16-recur step16-smoking
.PHONY: step16-preprocess step16-MutEnricher step16-tier1 step16-both step16-non step16-stage step16-gene step16-recur step16-smoking step16

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
Output/RSEM/Bowtie2 Output/RSEM/STAR Output/RSEM/DEG Output/RSEM/DEG-List:
	mkdir -p $@

Output/RSEM/Bowtie2/%.genome.bam Output/RSEM/Bowtie2/%.genes.results: Data/WTS/%_*.fastq.gz Output/RSEM/Reference/hg38.grp | Output/RSEM/Bowtie2 build.log
	$(DOCKER) rsem-calculate-expression --num-threads $(CPUS) --estimate-rspd --ci-memory $(shell echo "$(MEMS) * 1024" | bc) --bowtie2 --output-genome-bam --paired-end $(addprefix /,$(sort $(filter %.fastq.gz,$^)) $(basename $(filter %.grp,$^)) $(@D)/$*) 1> $@.stdout 2> $@.stderr

Output/RSEM/STAR/%.genome.bam Output/RSEM/STAR/%.genes.results: Data/WTS/%_*.fastq.gz Output/RSEM/Reference/hg38.grp | Output/RSEM/STAR build.log
	$(DOCKER) rsem-calculate-expression --num-threads $(CPUS) --estimate-rspd --ci-memory $(shell echo "$(MEMS) * 1024" | bc) --star --star-gzipped-read-file --output-genome-bam --paired-end $(addprefix /,$(sort $(filter %.fastq.gz,$^)) $(basename $(filter %.grp,$^)) $(@D)/$*) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/Bowtie2.tsv: Python/aggregate_RSEM.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(WTS_ID))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/Bowtie2.TPM.tsv: Python/aggregate_RSEM_4.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(WTS_ID))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/Bowtie2.FPKM.tsv: Python/aggregate_RSEM_5.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(WTS_ID))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/STAR.tsv: Python/aggregate_RSEM.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(WTS_ID))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/STAR.TPM.tsv: Python/aggregate_RSEM_4.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(WTS_ID))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/STAR.FPKM.tsv: Python/aggregate_RSEM_5.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(WTS_ID))) $(GENCODE_GTF) $(TREMBL_TSV) | Output/RSEM/DEG build.log
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

Output/RSEM/DEG/%.SQC.Smoking.Normal-CIS.tsv Output/RSEM/DEG/%.SQC.Smoking.Normal-CIS.coldata: Python/aggregate_RSEM_3.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --exact "Smoking" "Ever" --compare "stage" "Normal" "CIS+AIS" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Smoking.Normal-Primary.tsv Output/RSEM/DEG/%.SQC.Smoking.Normal-Primary.coldata: Python/aggregate_RSEM_3.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --exact "Smoking" "Ever" --compare "stage" "Normal" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Smoking.CIS-Primary.tsv Output/RSEM/DEG/%.SQC.Smoking.CIS-Primary.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --exact "Smoking" "Ever" --compare "stage" "CIS+AIS" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Nonsmoking.Normal-CIS.tsv Output/RSEM/DEG/%.SQC.Nonsmoking.Normal-CIS.coldata: Python/aggregate_RSEM_3.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --exact "Smoking" "Never" --compare "stage" "Normal" "CIS+AIS" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Nonsmoking.Normal-Primary.tsv Output/RSEM/DEG/%.SQC.Nonsmoking.Normal-Primary.coldata: Python/aggregate_RSEM_3.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --exact "Smoking" "Never" --compare "stage" "Normal" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Nonsmoking.CIS-Primary.tsv Output/RSEM/DEG/%.SQC.NonsmokingCIS-Primary.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --exact "Smoking" "Never" --compare "stage" "CIS+AIS" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.ADC.Nonsmoking.Normal-AIS.tsv Output/RSEM/DEG/%.ADC.Nonsmoking.Normal-AIS.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --ADC --exact "Smoking" "Never" --compare "stage" "Normal" "CIS+AIS" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.ADC.Nonsmoking.Normal-Primary.tsv Output/RSEM/DEG/%.ADC.Nonsmoking.Normal-Primary.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --ADC --exact "Smoking" "Never" --compare "stage" "Normal" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.ADC.Nonsmoking.AIS-Primary.tsv Output/RSEM/DEG/%.ADC.Nonsmoking.AIS-Primary.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --ADC --exact "Smoking" "Never" --compare "stage" "CIS+AIS" "Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.CIS.Nonsmoking-Smoking.tsv Output/RSEM/DEG/%.SQC.CIS.Nonsmoking-Smoking.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --exact "stage" "CIS+AIS" --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.SQC.Primary.Nonsmoking-Smoking.tsv Output/RSEM/DEG/%.SQC.Primary.Nonsmoking-Smoking.coldata: Python/aggregate_RSEM_2.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --SQC --exact "stage" "Primary" --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG/%.DEG.tsv: R/run_DESeq.R Output/RSEM/DEG/%.tsv Output/RSEM/DEG/%.coldata | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --readcount $(addprefix /,$(filter %.tsv,$^)) --coldata $(addprefix /,$(filter %.coldata,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG-List/%.Up.tex Output/RSEM/DEG-List/%.Up.tsv: Python/select_DEGs.py Output/RSEM/DEG/%.tsv | Output/RSEM/DEG-List build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --up 1> $@.stdout 2> $@.stderr

Output/RSEM/DEG-List/%.Down.tex Output/RSEM/DEG-List/%.Down.tsv: Python/select_DEGs.py Output/RSEM/DEG/%.tsv | Output/RSEM/DEG-List build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --down 1> $@.stdout 2> $@.stderr

step18-preparation: $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(WTS_ID))) Output/RSEM/DEG/Bowtie2.tsv Output/RSEM/DEG/Bowtie2.TPM.tsv Output/RSEM/DEG/Bowtie2.FPKM.tsv $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(WTS_ID))) Output/RSEM/DEG/STAR.tsv Output/RSEM/DEG/STAR.TPM.tsv Output/RSEM/DEG/STAR.FPKM.tsv
step18-SQC: Output/RSEM/DEG/Bowtie2.SQC.Normal-Dysplasia.tsv Output/RSEM/DEG/Bowtie2.SQC.Normal-CIS.tsv Output/RSEM/DEG/Bowtie2.SQC.Normal-Primary.tsv Output/RSEM/DEG/Bowtie2.SQC.Dysplasia-CIS.tsv Output/RSEM/DEG/Bowtie2.SQC.Dysplasia-Primary.tsv Output/RSEM/DEG/Bowtie2.SQC.CIS-Primary.tsv Output/RSEM/DEG/STAR.SQC.Normal-Dysplasia.tsv Output/RSEM/DEG/STAR.SQC.Normal-CIS.tsv Output/RSEM/DEG/STAR.SQC.Normal-Primary.tsv Output/RSEM/DEG/STAR.SQC.Dysplasia-CIS.tsv Output/RSEM/DEG/STAR.SQC.Dysplasia-Primary.tsv Output/RSEM/DEG/STAR.SQC.CIS-Primary.tsv
step18-SQC-DEG: Output/RSEM/DEG/Bowtie2.SQC.Normal-Dysplasia.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Normal-CIS.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Normal-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Dysplasia-CIS.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Dysplasia-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.CIS-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.Normal-Dysplasia.DEG.tsv Output/RSEM/DEG/STAR.SQC.Normal-CIS.DEG.tsv Output/RSEM/DEG/STAR.SQC.Normal-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.Dysplasia-CIS.DEG.tsv Output/RSEM/DEG/STAR.SQC.Dysplasia-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.CIS-Primary.DEG.tsv Output/RSEM/DEG-List/Bowtie2.SQC.Normal-Dysplasia.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.SQC.Normal-CIS.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.SQC.Normal-Primary.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.SQC.Dysplasia-CIS.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.SQC.Dysplasia-Primary.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.SQC.CIS-Primary.DEG.Up.tex Output/RSEM/DEG-List/STAR.SQC.Normal-Dysplasia.DEG.Up.tex Output/RSEM/DEG-List/STAR.SQC.Normal-CIS.DEG.Up.tex Output/RSEM/DEG-List/STAR.SQC.Normal-Primary.DEG.Up.tex Output/RSEM/DEG-List/STAR.SQC.Dysplasia-CIS.DEG.Up.tex Output/RSEM/DEG-List/STAR.SQC.Dysplasia-Primary.DEG.Up.tex Output/RSEM/DEG-List/STAR.SQC.CIS-Primary.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.SQC.Normal-Dysplasia.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.SQC.Normal-CIS.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.SQC.Normal-Primary.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.SQC.Dysplasia-CIS.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.SQC.Dysplasia-Primary.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.SQC.CIS-Primary.DEG.Down.tex Output/RSEM/DEG-List/STAR.SQC.Normal-Dysplasia.DEG.Down.tex Output/RSEM/DEG-List/STAR.SQC.Normal-CIS.DEG.Down.tex Output/RSEM/DEG-List/STAR.SQC.Normal-Primary.DEG.Down.tex Output/RSEM/DEG-List/STAR.SQC.Dysplasia-CIS.DEG.Down.tex Output/RSEM/DEG-List/STAR.SQC.Dysplasia-Primary.DEG.Down.tex Output/RSEM/DEG-List/STAR.SQC.CIS-Primary.DEG.Down.tex
step18-SQC-Recur: Output/RSEM/DEG/Bowtie2.SQC.Recur.Normal-CIS.tsv Output/RSEM/DEG/Bowtie2.SQC.Recur.Normal-Primary.tsv Output/RSEM/DEG/Bowtie2.SQC.Recur.CIS-Primary.tsv Output/RSEM/DEG/STAR.SQC.Recur.Normal-CIS.tsv Output/RSEM/DEG/STAR.SQC.Recur.Normal-Primary.tsv Output/RSEM/DEG/STAR.SQC.Recur.CIS-Primary.tsv Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.Normal-CIS.tsv Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.Normal-Primary.tsv Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.CIS-Primary.tsv Output/RSEM/DEG/STAR.SQC.Nonrecur.Normal-CIS.tsv Output/RSEM/DEG/STAR.SQC.Nonrecur.Normal-Primary.tsv Output/RSEM/DEG/STAR.SQC.Nonrecur.CIS-Primary.tsv
step18-SQC-Recur-DEG: Output/RSEM/DEG/Bowtie2.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Recur.CIS-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/STAR.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.Recur.CIS-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Nonrecur.CIS-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/STAR.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.Nonrecur.CIS-Primary.DEG.tsv Output/RSEM/DEG-List/Bowtie2.SQC.Recur.Normal-CIS.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.SQC.Recur.Normal-Primary.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.SQC.Recur.CIS-Primary.DEG.Up.tex Output/RSEM/DEG-List/STAR.SQC.Recur.Normal-CIS.DEG.Up.tex Output/RSEM/DEG-List/STAR.SQC.Recur.Normal-Primary.DEG.Up.tex Output/RSEM/DEG-List/STAR.SQC.Recur.CIS-Primary.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.SQC.Nonrecur.Normal-CIS.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.SQC.Nonrecur.Normal-Primary.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.SQC.Nonrecur.CIS-Primary.DEG.Up.tex Output/RSEM/DEG-List/STAR.SQC.Nonrecur.Normal-CIS.DEG.Up.tex Output/RSEM/DEG-List/STAR.SQC.Nonrecur.Normal-Primary.DEG.Up.tex Output/RSEM/DEG-List/STAR.SQC.Nonrecur.CIS-Primary.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.SQC.Recur.Normal-CIS.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.SQC.Recur.Normal-Primary.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.SQC.Recur.CIS-Primary.DEG.Down.tex Output/RSEM/DEG-List/STAR.SQC.Recur.Normal-CIS.DEG.Down.tex Output/RSEM/DEG-List/STAR.SQC.Recur.Normal-Primary.DEG.Down.tex Output/RSEM/DEG-List/STAR.SQC.Recur.CIS-Primary.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.SQC.Nonrecur.Normal-CIS.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.SQC.Nonrecur.Normal-Primary.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.SQC.Nonrecur.CIS-Primary.DEG.Down.tex Output/RSEM/DEG-List/STAR.SQC.Nonrecur.Normal-CIS.DEG.Down.tex Output/RSEM/DEG-List/STAR.SQC.Nonrecur.Normal-Primary.DEG.Down.tex Output/RSEM/DEG-List/STAR.SQC.Nonrecur.CIS-Primary.DEG.Down.tex
step18-SQC-Smoking: Output/RSEM/DEG/Bowtie2.SQC.Smoking.Normal-CIS.tsv Output/RSEM/DEG/Bowtie2.SQC.Smoking.Normal-Primary.tsv Output/RSEM/DEG/Bowtie2.SQC.Smoking.CIS-Primary.tsv Output/RSEM/DEG/Bowtie2.SQC.Nonsmoking.Normal-CIS.tsv Output/RSEM/DEG/Bowtie2.SQC.Nonsmoking.Normal-Primary.tsv Output/RSEM/DEG/Bowtie2.SQC.CIS.Nonsmoking-Smoking.tsv Output/RSEM/DEG/Bowtie2.SQC.Primary.Nonsmoking-Smoking.tsv Output/RSEM/DEG/STAR.SQC.Smoking.Normal-CIS.tsv Output/RSEM/DEG/STAR.SQC.Smoking.Normal-Primary.tsv Output/RSEM/DEG/STAR.SQC.Smoking.CIS-Primary.tsv Output/RSEM/DEG/STAR.SQC.Nonsmoking.Normal-CIS.tsv Output/RSEM/DEG/STAR.SQC.Nonsmoking.Normal-Primary.tsv Output/RSEM/DEG/STAR.SQC.CIS.Nonsmoking-Smoking.tsv Output/RSEM/DEG/STAR.SQC.Primary.Nonsmoking-Smoking.tsv
step18-SQC-Smoking-DEG: Output/RSEM/DEG/Bowtie2.SQC.Smoking.Normal-CIS.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Smoking.Normal-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Smoking.CIS-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Nonsmoking.Normal-CIS.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Nonsmoking.Normal-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.CIS.Nonsmoking-Smoking.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Primary.Nonsmoking-Smoking.DEG.tsv Output/RSEM/DEG/STAR.SQC.Smoking.Normal-CIS.DEG.tsv Output/RSEM/DEG/STAR.SQC.Smoking.Normal-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.Smoking.CIS-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.Nonsmoking.Normal-CIS.DEG.tsv Output/RSEM/DEG/STAR.SQC.Nonsmoking.Normal-Primary.DEG.tsv Output/RSEM/DEG/STAR.SQC.CIS.Nonsmoking-Smoking.DEG.tsv Output/RSEM/DEG/STAR.SQC.Primary.Nonsmoking-Smoking.DEG.tsv
step18-SQC-Stage: Output/RSEM/DEG/Bowtie2.SQC.CIS.Nonrecur-Recur.tsv Output/RSEM/DEG/Bowtie2.SQC.Primary.Nonrecur-Recur.tsv Output/RSEM/DEG/STAR.SQC.CIS.Nonrecur-Recur.tsv Output/RSEM/DEG/STAR.SQC.Primary.Nonrecur-Recur.tsv
step18-SQC-Stage-DEG: Output/RSEM/DEG/Bowtie2.SQC.CIS.Nonrecur-Recur.DEG.tsv Output/RSEM/DEG/Bowtie2.SQC.Primary.Nonrecur-Recur.DEG.tsv Output/RSEM/DEG/STAR.SQC.CIS.Nonrecur-Recur.DEG.tsv Output/RSEM/DEG/STAR.SQC.Primary.Nonrecur-Recur.DEG.tsv Output/RSEM/DEG-List/Bowtie2.SQC.CIS.Nonrecur-Recur.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.SQC.Primary.Nonrecur-Recur.DEG.Up.tex Output/RSEM/DEG-List/STAR.SQC.CIS.Nonrecur-Recur.DEG.Up.tex Output/RSEM/DEG-List/STAR.SQC.Primary.Nonrecur-Recur.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.SQC.CIS.Nonrecur-Recur.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.SQC.Primary.Nonrecur-Recur.DEG.Down.tex Output/RSEM/DEG-List/STAR.SQC.CIS.Nonrecur-Recur.DEG.Down.tex Output/RSEM/DEG-List/STAR.SQC.Primary.Nonrecur-Recur.DEG.Down.tex
step18-ADC: Output/RSEM/DEG/Bowtie2.ADC.Normal-AAH.tsv Output/RSEM/DEG/Bowtie2.ADC.Normal-AIS.tsv Output/RSEM/DEG/Bowtie2.ADC.Normal-Primary.tsv Output/RSEM/DEG/Bowtie2.ADC.AAH-AIS.tsv Output/RSEM/DEG/Bowtie2.ADC.AAH-Primary.tsv Output/RSEM/DEG/Bowtie2.ADC.AIS-Primary.tsv Output/RSEM/DEG/STAR.ADC.Normal-AAH.tsv Output/RSEM/DEG/STAR.ADC.Normal-AIS.tsv Output/RSEM/DEG/STAR.ADC.Normal-Primary.tsv Output/RSEM/DEG/STAR.ADC.AAH-AIS.tsv Output/RSEM/DEG/STAR.ADC.AAH-Primary.tsv Output/RSEM/DEG/STAR.ADC.AIS-Primary.tsv
step18-ADC-DEG: Output/RSEM/DEG/Bowtie2.ADC.Normal-AAH.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.Normal-AIS.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.Normal-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.AAH-AIS.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.AAH-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.AIS-Primary.DEG.tsv Output/RSEM/DEG/STAR.ADC.Normal-AAH.DEG.tsv Output/RSEM/DEG/STAR.ADC.Normal-AIS.DEG.tsv Output/RSEM/DEG/STAR.ADC.Normal-Primary.DEG.tsv Output/RSEM/DEG/STAR.ADC.AAH-AIS.DEG.tsv Output/RSEM/DEG/STAR.ADC.AAH-Primary.DEG.tsv Output/RSEM/DEG/STAR.ADC.AIS-Primary.DEG.tsv Output/RSEM/DEG-List/Bowtie2.ADC.Normal-AAH.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.ADC.Normal-AIS.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.ADC.Normal-Primary.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.ADC.AAH-AIS.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.ADC.AAH-Primary.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.ADC.AIS-Primary.DEG.Up.tex Output/RSEM/DEG-List/STAR.ADC.Normal-AAH.DEG.Up.tex Output/RSEM/DEG-List/STAR.ADC.Normal-AIS.DEG.Up.tex Output/RSEM/DEG-List/STAR.ADC.Normal-Primary.DEG.Up.tex Output/RSEM/DEG-List/STAR.ADC.AAH-AIS.DEG.Up.tex Output/RSEM/DEG-List/STAR.ADC.AAH-Primary.DEG.Up.tex Output/RSEM/DEG-List/STAR.ADC.AIS-Primary.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.ADC.Normal-AAH.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.ADC.Normal-AIS.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.ADC.Normal-Primary.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.ADC.AAH-AIS.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.ADC.AAH-Primary.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.ADC.AIS-Primary.DEG.Down.tex Output/RSEM/DEG-List/STAR.ADC.Normal-AAH.DEG.Down.tex Output/RSEM/DEG-List/STAR.ADC.Normal-AIS.DEG.Down.tex Output/RSEM/DEG-List/STAR.ADC.Normal-Primary.DEG.Down.tex Output/RSEM/DEG-List/STAR.ADC.AAH-AIS.DEG.Down.tex Output/RSEM/DEG-List/STAR.ADC.AAH-Primary.DEG.Down.tex Output/RSEM/DEG-List/STAR.ADC.AIS-Primary.DEG.Down.tex
step18-ADC-Recur: Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.Normal-AIS.tsv Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.Normal-Primary.tsv Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.AIS-Primary.tsv Output/RSEM/DEG/STAR.ADC.Nonrecur.Normal-AIS.tsv Output/RSEM/DEG/STAR.ADC.Nonrecur.Normal-Primary.tsv Output/RSEM/DEG/STAR.ADC.Nonrecur.AIS-Primary.tsv
step18-ADC-Recur-DEG: Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.Normal-AIS.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.Nonrecur.AIS-Primary.DEG.tsv Output/RSEM/DEG/STAR.ADC.Nonrecur.Normal-AIS.DEG.tsv Output/RSEM/DEG/STAR.ADC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/STAR.ADC.Nonrecur.AIS-Primary.DEG.tsv Output/RSEM/DEG-List/Bowtie2.ADC.Nonrecur.Normal-AIS.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.ADC.Nonrecur.Normal-Primary.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.ADC.Nonrecur.AIS-Primary.DEG.Up.tex Output/RSEM/DEG-List/STAR.ADC.Nonrecur.Normal-AIS.DEG.Up.tex Output/RSEM/DEG-List/STAR.ADC.Nonrecur.Normal-Primary.DEG.Up.tex Output/RSEM/DEG-List/STAR.ADC.Nonrecur.AIS-Primary.DEG.Up.tex Output/RSEM/DEG-List/Bowtie2.ADC.Nonrecur.Normal-AIS.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.ADC.Nonrecur.Normal-Primary.DEG.Down.tex Output/RSEM/DEG-List/Bowtie2.ADC.Nonrecur.AIS-Primary.DEG.Down.tex Output/RSEM/DEG-List/STAR.ADC.Nonrecur.Normal-AIS.DEG.Down.tex Output/RSEM/DEG-List/STAR.ADC.Nonrecur.Normal-Primary.DEG.Down.tex Output/RSEM/DEG-List/STAR.ADC.Nonrecur.AIS-Primary.DEG.Down.tex
step18-ADC-Smoking: Output/RSEM/DEG/Bowtie2.ADC.Nonsmoking.Normal-AIS.tsv Output/RSEM/DEG/Bowtie2.ADC.Nonsmoking.Normal-Primary.tsv Output/RSEM/DEG/Bowtie2.ADC.Nonsmoking.AIS-Primary.tsv Output/RSEM/DEG/STAR.ADC.Nonsmoking.Normal-AIS.tsv Output/RSEM/DEG/STAR.ADC.Nonsmoking.Normal-Primary.tsv Output/RSEM/DEG/STAR.ADC.Nonsmoking.AIS-Primary.tsv
step18-ADC-Smoking-DEG: Output/RSEM/DEG/Bowtie2.ADC.Nonsmoking.Normal-AIS.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.Nonsmoking.Normal-Primary.DEG.tsv Output/RSEM/DEG/Bowtie2.ADC.Nonsmoking.AIS-Primary.DEG.tsv Output/RSEM/DEG/STAR.ADC.Nonsmoking.Normal-AIS.DEG.tsv Output/RSEM/DEG/STAR.ADC.Nonsmoking.Normal-Primary.DEG.tsv Output/RSEM/DEG/STAR.ADC.Nonsmoking.AIS-Primary.DEG.tsv
step18: step18-preparation step18-SQC step18-SQC-DEG step18-SQC-Recur step18-SQC-Recur-DEG step18-SQC-Stage step18-SQC-Stage-DEG step18-SQC-Smoking step18-SQC-Smoking-DEG step18-ADC step18-ADC-DEG step18-ADC-Recur step18-ADC-Recur-DEG step18-ADC-Smoking step18-ADC-Smoking-DEG
.PHONY: step18-preparation step18-SQC step18-SQC-DEG step18-SQC-Recur step18-SQC-Recur-DEG step18-SQC-Stage step18-SQC-Stage-DEG step18-SQC-Smoking step18-SQC-Smoking-DEG step18-ADC step18-ADC-DEG step18-ADC-Recur step18-ADC-Recur-DEG step18-ADC-Smoking step18-ADC-Smoking-DEG step18

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

Output/RSEM/Plots/Venn/%.ADC.venn.pdf: Python/draw_venn_plot_3.py Output/RSEM/DEG/%.ADC.Normal-AIS.DEG.tsv Output/RSEM/DEG/%.ADC.Normal-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-AIS" "Normal-Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.ADC.Nonrecur.Up.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.ADC.Nonrecur.Normal-AIS.DEG.tsv Output/RSEM/DEG/%.ADC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.ADC.Nonrecur.AIS-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-AIS" "Normal-Primary" "AIS-Primary" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.ADC.Nonrecur.Down.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.ADC.Nonrecur.Normal-AIS.DEG.tsv Output/RSEM/DEG/%.ADC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.ADC.Nonrecur.AIS-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-AIS" "Normal-Primary" "AIS-Primary" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.ADC.Nonrecur.venn.pdf: Python/draw_venn_plot_3.py Output/RSEM/DEG/%.ADC.Nonrecur.Normal-AIS.DEG.tsv Output/RSEM/DEG/%.ADC.Nonrecur.Normal-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-AIS" "Normal-Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.Up.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Normal-Dysplasia.DEG.tsv Output/RSEM/DEG/%.SQC.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Normal-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-Dysplasia" "Normal-CIS" "Normal-Primary" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.Down.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Normal-Dysplasia.DEG.tsv Output/RSEM/DEG/%.SQC.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Normal-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-Dysplasia" "Normal-CIS" "Normal-Primary" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.venn.pdf: Python/draw_venn_plot_3.py Output/RSEM/DEG/%.SQC.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Normal-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-CIS" "Normal-Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.Recur.Up.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.CIS-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-CIS" "Normal-Primary" "CIS-Primary" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.Recur.Down.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.CIS-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-CIS" "Normal-Primary" "CIS-Primary" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.Recur.venn.pdf: Python/draw_venn_plot_3.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-CIS" "Normal-Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.Nonrecur.Up.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.CIS-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-CIS" "Normal-Primary" "CIS-Primary" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.Nonrecur.Down.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.CIS-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-CIS" "Normal-Primary" "CIS-Primary" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.Nonrecur.venn.pdf: Python/draw_venn_plot_3.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Normal-CIS" "Normal-Primary" 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC-CIS.Recur-Nonrecur.Up.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Recur" "Non-recur" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC-CIS.Recur-Nonrecur.Down.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Recur" "Non-recur" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC-CIS.Recur-Nonrecur.venn.pdf: Python/draw_venn_plot_3.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Recur" "Non-recur" 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC-Primary.Recur-Nonrecur.Up.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Recur" "Non-recur" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC-Primary.Recur-Nonrecur.Down.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Recur" "Non-recur" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC-Primary.Recur-Nonrecur.venn.pdf: Python/draw_venn_plot_3.py Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "Recur" "Non-recur" 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.CIS-Primary.Recur-Nonrecur.Up.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.CIS.Nonrecur-Recur.DEG.tsv Output/RSEM/DEG/%.SQC.Primary.Nonrecur-Recur.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "CIS" "Primary" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Venn/%.SQC.CIS-Primary.Recur-Nonrecur.Down.venn.pdf: Python/draw_venn_plot.py Output/RSEM/DEG/%.SQC.CIS.Nonrecur-Recur.DEG.tsv Output/RSEM/DEG/%.SQC.Primary.Nonrecur-Recur.DEG.tsv | Output/RSEM/Plots/Venn build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "CIS" "Primary" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Volcano/%.volcano.pdf: Python/draw_volcano_plot.py Output/RSEM/DEG/%.DEG.tsv | Output/RSEM/Plots/Volcano build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Violin/%.ADC.violin.tar: Python/draw_violin_plots.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | Output/RSEM/Plots/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Violin/%.SQC.violin.tar: Python/draw_violin_plots.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | Output/RSEM/Plots/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Violin/%.ADC.Recurrence.violin.tar: Python/draw_violin_plots_clinical.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | Output/RSEM/Plots/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Violin/%.SQC.Recurrence.violin.tar: Python/draw_violin_plots_clinical.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | Output/RSEM/Plots/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Violin/%.ADC.Smoking.violin.tar: Python/draw_violin_plots_clinical.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | Output/RSEM/Plots/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/Violin/%.SQC.Smoking.violin.tar: Python/draw_violin_plots_clinical.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | Output/RSEM/Plots/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.Up.KEGG.tsv Output/RSEM/Enrichment/%.Up.KEGG.tex Output/RSEM/Enrichment/%.Up.KEGG.pdf: Python/get_enrichment_enrichr.py Output/RSEM/DEG/%.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.Down.KEGG.tsv Output/RSEM/Enrichment/%.Down.KEGG.tex Output/RSEM/Enrichment/%.Down.KEGG.pdf: Python/get_enrichment_enrichr.py Output/RSEM/DEG/%.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.Up.MSigDB.tsv Output/RSEM/Enrichment/%.Up.MSigDB.tex Output/RSEM/Enrichment/%.Up.MSigDB.pdf: Python/get_enrichment_enrichr.py Output/RSEM/DEG/%.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --up 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.Down.MSigDB.tsv Output/RSEM/Enrichment/%.Down.MSigDB.tex Output/RSEM/Enrichment/%.Down.MSigDB.pdf: Python/get_enrichment_enrichr.py Output/RSEM/DEG/%.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --down 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Recur.Up.KEGG.tsv Output/RSEM/Enrichment/%.SQC-CIS.Recur.Up.KEGG.tex Output/RSEM/Enrichment/%.SQC-CIS.Recur.Up.KEGG.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --up --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Recur.Down.KEGG.tsv Output/RSEM/Enrichment/%.SQC-CIS.Recur.Down.KEGG.tex Output/RSEM/Enrichment/%.SQC-CIS.Recur.Down.KEGG.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --down --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Recur.Up.MSigDB.tsv Output/RSEM/Enrichment/%.SQC-CIS.Recur.Up.MSigDB.tex Output/RSEM/Enrichment/%.SQC-CIS.Recur.Up.MSigDB.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --up --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Recur.Down.MSigDB.tsv Output/RSEM/Enrichment/%.SQC-CIS.Recur.Down.MSigDB.tex Output/RSEM/Enrichment/%.SQC-CIS.Recur.Down.MSigDB.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --down --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Nonrecur.Up.KEGG.tsv Output/RSEM/Enrichment/%.SQC-CIS.Nonrecur.Up.KEGG.tex Output/RSEM/Enrichment/%.SQC-CIS.Nonrecur.Up.KEGG.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --up --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Nonrecur.Down.KEGG.tsv Output/RSEM/Enrichment/%.SQC-CIS.Nonrecur.Down.KEGG.tex Output/RSEM/Enrichment/%.SQC-CIS.Nonrecur.Down.KEGG.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --down --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Nonrecur.Up.MSigDB.tsv Output/RSEM/Enrichment/%.SQC-CIS.Nonrecur.Up.MSigDB.tex Output/RSEM/Enrichment/%.SQC-CIS.Nonrecur.Up.MSigDB.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --up --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Nonrecur.Down.MSigDB.tsv Output/RSEM/Enrichment/%.SQC-CIS.Nonrecur.Down.MSigDB.tex Output/RSEM/Enrichment/%.SQC-CIS.Nonrecur.Down.MSigDB.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --down --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Both.Up.KEGG.tsv Output/RSEM/Enrichment/%.SQC-CIS.Both.Up.KEGG.tex Output/RSEM/Enrichment/%.SQC-CIS.Both.Up.KEGG.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --up --intersect 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Both.Down.KEGG.tsv Output/RSEM/Enrichment/%.SQC-CIS.Both.Down.KEGG.tex Output/RSEM/Enrichment/%.SQC-CIS.Both.Down.KEGG.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --down --intersect 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Both.Up.MSigDB.tsv Output/RSEM/Enrichment/%.SQC-CIS.Both.Up.MSigDB.tex Output/RSEM/Enrichment/%.SQC-CIS.Both.Up.MSigDB.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --up --intersect 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-CIS.Both.Down.MSigDB.tsv Output/RSEM/Enrichment/%.SQC-CIS.Both.Down.MSigDB.tex Output/RSEM/Enrichment/%.SQC-CIS.Both.Down.MSigDB.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-CIS.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-CIS.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --down --intersect 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Recur.Up.KEGG.tsv Output/RSEM/Enrichment/%.SQC-Primary.Recur.Up.KEGG.tex Output/RSEM/Enrichment/%.SQC-Primary.Recur.Up.KEGG.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --up --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Recur.Down.KEGG.tsv Output/RSEM/Enrichment/%.SQC-Primary.Recur.Down.KEGG.tex Output/RSEM/Enrichment/%.SQC-Primary.Recur.Down.KEGG.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --down --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Recur.Up.MSigDB.tsv Output/RSEM/Enrichment/%.SQC-Primary.Recur.Up.MSigDB.tex Output/RSEM/Enrichment/%.SQC-Primary.Recur.Up.MSigDB.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --up --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Recur.Down.MSigDB.tsv Output/RSEM/Enrichment/%.SQC-Primary.Recur.Down.MSigDB.tex Output/RSEM/Enrichment/%.SQC-Primary.Recur.Down.MSigDB.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --down --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Nonrecur.Up.KEGG.tsv Output/RSEM/Enrichment/%.SQC-Primary.Nonrecur.Up.KEGG.tex Output/RSEM/Enrichment/%.SQC-Primary.Nonrecur.Up.KEGG.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --up --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Nonrecur.Down.KEGG.tsv Output/RSEM/Enrichment/%.SQC-Primary.Nonrecur.Down.KEGG.tex Output/RSEM/Enrichment/%.SQC-Primary.Nonrecur.Down.KEGG.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --down --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Nonrecur.Up.MSigDB.tsv Output/RSEM/Enrichment/%.SQC-Primary.Nonrecur.Up.MSigDB.tex Output/RSEM/Enrichment/%.SQC-Primary.Nonrecur.Up.MSigDB.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --up --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Nonrecur.Down.MSigDB.tsv Output/RSEM/Enrichment/%.SQC-Primary.Nonrecur.Down.MSigDB.tex Output/RSEM/Enrichment/%.SQC-Primary.Nonrecur.Down.MSigDB.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --down --subtract 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Both.Up.KEGG.tsv Output/RSEM/Enrichment/%.SQC-Primary.Both.Up.KEGG.tex Output/RSEM/Enrichment/%.SQC-Primary.Both.Up.KEGG.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --up --intersect 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Both.Down.KEGG.tsv Output/RSEM/Enrichment/%.SQC-Primary.Both.Down.KEGG.tex Output/RSEM/Enrichment/%.SQC-Primary.Both.Down.KEGG.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "KEGG_2021_Human" --down --intersect 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Both.Up.MSigDB.tsv Output/RSEM/Enrichment/%.SQC-Primary.Both.Up.MSigDB.tex Output/RSEM/Enrichment/%.SQC-Primary.Both.Up.MSigDB.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --up --intersect 1> $@.stdout 2> $@.stderr

Output/RSEM/Enrichment/%.SQC-Primary.Both.Down.MSigDB.tsv Output/RSEM/Enrichment/%.SQC-Primary.Both.Down.MSigDB.tex Output/RSEM/Enrichment/%.SQC-Primary.Both.Down.MSigDB.pdf: Python/get_enrichment_enrichr_2.py Output/RSEM/DEG/%.SQC.Nonrecur.Normal-Primary.DEG.tsv Output/RSEM/DEG/%.SQC.Recur.Normal-Primary.DEG.tsv | Output/RSEM/Enrichment build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --DB "MSigDB_Oncogenic_Signatures" --down --intersect 1> $@.stdout 2> $@.stderr

Output/RSEM/Plots/SQC.%.pair.pdf: Python/draw_pair_plot.py Output/RSEM/DEG/SQC.%.tsv Output/RSEM/DEG/SQC.%.DEG.tsv | Output/RSEM/Plots build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC

step19-cluster: Output/RSEM/Plots/Cluster/Bowtie2.SQC.pdf Output/RSEM/Plots/Cluster/STAR.SQC.pdf Output/RSEM/Plots/Cluster/Bowtie2.ADC.pdf Output/RSEM/Plots/Cluster/STAR.ADC.pdf Output/RSEM/Plots/Cluster/Bowtie2.ADC-SQC.pdf Output/RSEM/Plots/Cluster/STAR.ADC-SQC.pdf
step19-venn-ADC: Output/RSEM/Plots/Venn/Bowtie2.ADC.Up.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.ADC.Down.venn.pdf Output/RSEM/Plots/Venn/STAR.ADC.Up.venn.pdf Output/RSEM/Plots/Venn/STAR.ADC.Down.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.ADC.venn.pdf Output/RSEM/Plots/Venn/STAR.ADC.venn.pdf
step19-venn-ADC-Recur: Output/RSEM/Plots/Venn/Bowtie2.ADC.Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.ADC.Nonrecur.Down.venn.pdf Output/RSEM/Plots/Venn/STAR.ADC.Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/STAR.ADC.Nonrecur.Down.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.ADC.Nonrecur.venn.pdf Output/RSEM/Plots/Venn/STAR.ADC.Nonrecur.venn.pdf
step19-venn-SQC: Output/RSEM/Plots/Venn/Bowtie2.SQC.Up.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC.Down.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.Up.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.Down.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.venn.pdf
step19-venn-SQC-Recur: Output/RSEM/Plots/Venn/Bowtie2.SQC.Recur.Up.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC.Recur.Down.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.Recur.Up.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.Recur.Down.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC.Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC.Nonrecur.Down.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.Nonrecur.Down.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC-CIS.Recur-Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC-CIS.Recur-Nonrecur.Down.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC-Primary.Recur-Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC-Primary.Recur-Nonrecur.Down.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC-CIS.Recur-Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC-CIS.Recur-Nonrecur.Down.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC-Primary.Recur-Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC-Primary.Recur-Nonrecur.Down.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC.CIS-Primary.Recur-Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC.CIS-Primary.Recur-Nonrecur.Down.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.CIS-Primary.Recur-Nonrecur.Up.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.CIS-Primary.Recur-Nonrecur.Down.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC.Recur.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.Recur.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC.Nonrecur.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC.Nonrecur.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC-CIS.Recur-Nonrecur.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC-CIS.Recur-Nonrecur.venn.pdf Output/RSEM/Plots/Venn/Bowtie2.SQC-Primary.Recur-Nonrecur.venn.pdf Output/RSEM/Plots/Venn/STAR.SQC-Primary.Recur-Nonrecur.venn.pdf
step19-volcano-ADC: Output/RSEM/Plots/Volcano/Bowtie2.ADC.Normal-AAH.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.ADC.Normal-AIS.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.ADC.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.ADC.AAH-AIS.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.ADC.AAH-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.ADC.AIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.Normal-AAH.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.Normal-AIS.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.AAH-AIS.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.AAH-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.AIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.ADC.Nonrecur.Normal-AIS.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.ADC.Nonrecur.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.ADC.Nonrecur.AIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.Nonrecur.Normal-AIS.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.Nonrecur.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.ADC.Nonrecur.AIS-Primary.volcano.pdf
step19-volcano-SQC: Output/RSEM/Plots/Volcano/Bowtie2.SQC.Normal-Dysplasia.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Normal-CIS.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Dysplasia-CIS.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Dysplasia-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.CIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Normal-Dysplasia.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Normal-CIS.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Dysplasia-CIS.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Dysplasia-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.CIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Recur.Normal-CIS.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Recur.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Recur.Normal-CIS.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Recur.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Recur.Normal-CIS.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Recur.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Recur.CIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Recur.Normal-CIS.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Recur.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Recur.CIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Nonrecur.Normal-CIS.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Nonrecur.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Nonrecur.CIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Nonrecur.Normal-CIS.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Nonrecur.Normal-Primary.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Nonrecur.CIS-Primary.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.CIS.Nonrecur-Recur.volcano.pdf Output/RSEM/Plots/Volcano/Bowtie2.SQC.Primary.Nonrecur-Recur.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.CIS.Nonrecur-Recur.volcano.pdf Output/RSEM/Plots/Volcano/STAR.SQC.Primary.Nonrecur-Recur.volcano.pdf
step19-violin-ADC: Output/RSEM/Plots/Violin/Bowtie2.TPM.ADC.violin.tar Output/RSEM/Plots/Violin/STAR.TPM.ADC.violin.tar Output/RSEM/Plots/Violin/Bowtie2.FPKM.ADC.violin.tar Output/RSEM/Plots/Violin/STAR.FPKM.ADC.violin.tar Output/RSEM/Plots/Violin/Bowtie2.TPM.ADC.Recurrence.violin.tar Output/RSEM/Plots/Violin/STAR.TPM.ADC.Recurrence.violin.tar Output/RSEM/Plots/Violin/Bowtie2.FPKM.ADC.Recurrence.violin.tar Output/RSEM/Plots/Violin/STAR.FPKM.ADC.Recurrence.violin.tar Output/RSEM/Plots/Violin/Bowtie2.TPM.ADC.Smoking.violin.tar Output/RSEM/Plots/Violin/STAR.TPM.ADC.Smoking.violin.tar Output/RSEM/Plots/Violin/Bowtie2.FPKM.ADC.Smoking.violin.tar Output/RSEM/Plots/Violin/STAR.FPKM.ADC.Smoking.violin.tar
step19-violin-SQC: Output/RSEM/Plots/Violin/Bowtie2.TPM.SQC.violin.tar Output/RSEM/Plots/Violin/STAR.TPM.SQC.violin.tar Output/RSEM/Plots/Violin/Bowtie2.FPKM.SQC.violin.tar Output/RSEM/Plots/Violin/STAR.FPKM.SQC.violin.tar Output/RSEM/Plots/Violin/Bowtie2.TPM.SQC.Recurrence.violin.tar Output/RSEM/Plots/Violin/STAR.TPM.SQC.Recurrence.violin.tar Output/RSEM/Plots/Violin/Bowtie2.FPKM.SQC.Recurrence.violin.tar Output/RSEM/Plots/Violin/STAR.FPKM.SQC.Recurrence.violin.tar Output/RSEM/Plots/Violin/Bowtie2.TPM.SQC.Smoking.violin.tar Output/RSEM/Plots/Violin/STAR.TPM.SQC.Smoking.violin.tar Output/RSEM/Plots/Violin/Bowtie2.FPKM.SQC.Smoking.violin.tar Output/RSEM/Plots/Violin/STAR.FPKM.SQC.Smoking.violin.tar
step19-KEGG-ADC: Output/RSEM/Enrichment/Bowtie2.ADC.Normal-AAH.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-AIS.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AAH-AIS.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AAH-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AIS-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-AAH.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-AIS.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AAH-AIS.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AAH-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AIS-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-AAH.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-AIS.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.AAH-AIS.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.AAH-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.AIS-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-AAH.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-AIS.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.AAH-AIS.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.AAH-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.AIS-Primary.Down.KEGG.tsv
step19-MSigDB-ADC: Output/RSEM/Enrichment/Bowtie2.ADC.Normal-AAH.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-AIS.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AAH-AIS.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AAH-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AIS-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-AAH.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-AIS.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AAH-AIS.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AAH-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.AIS-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-AAH.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-AIS.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.AAH-AIS.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.AAH-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.AIS-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-AAH.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-AIS.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.AAH-AIS.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.AAH-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.AIS-Primary.Down.MSigDB.tsv
step19-KEGG-ADC-Recur: Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.Normal-AIS.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.AIS-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.Normal-AIS.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.AIS-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.Normal-AIS.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.AIS-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.Normal-AIS.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.AIS-Primary.Down.KEGG.tsv
step19-MSigDB-ADC-Recur: Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.Normal-AIS.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.AIS-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.Normal-AIS.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.ADC.Nonrecur.AIS-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.Normal-AIS.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.AIS-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.Normal-AIS.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.ADC.Nonrecur.AIS-Primary.Down.MSigDB.tsv
step19-KEGG-SQC: Output/RSEM/Enrichment/Bowtie2.SQC.Normal-Dysplasia.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-CIS.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Dysplasia-CIS.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Dysplasia-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.CIS-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-Dysplasia.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-CIS.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Dysplasia-CIS.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Dysplasia-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.CIS-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-Dysplasia.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-CIS.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Dysplasia-CIS.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Dysplasia-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.CIS-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-Dysplasia.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-CIS.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Dysplasia-CIS.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Dysplasia-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.CIS-Primary.Down.KEGG.tsv
step19-MSigDB-SQC: Output/RSEM/Enrichment/Bowtie2.SQC.Normal-Dysplasia.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-CIS.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Dysplasia-CIS.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Dysplasia-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.CIS-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-Dysplasia.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-CIS.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Dysplasia-CIS.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Dysplasia-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.CIS-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-Dysplasia.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-CIS.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Dysplasia-CIS.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Dysplasia-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.CIS-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-Dysplasia.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-CIS.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Dysplasia-CIS.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Dysplasia-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.CIS-Primary.Down.MSigDB.tsv
step19-KEGG-SQC-Recur: Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Recur.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Recur.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Nonrecur.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Nonrecur.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Recur.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Recur.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Nonrecur.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Nonrecur.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Recur.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Recur.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Nonrecur.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Nonrecur.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Recur.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Recur.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Nonrecur.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Nonrecur.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Recur.Normal-CIS.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Recur.Normal-CIS.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Recur.Normal-CIS.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Recur.Normal-CIS.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Nonrecur.Normal-CIS.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Nonrecur.Normal-CIS.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Nonrecur.Normal-CIS.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Nonrecur.Normal-CIS.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Recur.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Recur.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Recur.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Recur.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Nonrecur.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Nonrecur.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Nonrecur.Normal-Primary.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC.Nonrecur.Normal-Primary.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Both.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Both.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Both.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Both.Down.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Both.Up.KEGG.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Both.Down.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Both.Up.KEGG.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Both.Down.KEGG.tsv
step19-MSigDB-SQC-Recur: Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Recur.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Recur.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Nonrecur.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Nonrecur.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Recur.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Recur.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Nonrecur.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Nonrecur.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Recur.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Recur.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Nonrecur.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Nonrecur.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Recur.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Recur.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Nonrecur.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Nonrecur.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Recur.Normal-CIS.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Recur.Normal-CIS.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Recur.Normal-CIS.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Recur.Normal-CIS.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Nonrecur.Normal-CIS.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Nonrecur.Normal-CIS.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Nonrecur.Normal-CIS.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Nonrecur.Normal-CIS.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Recur.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Recur.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Recur.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Recur.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Nonrecur.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC.Nonrecur.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Nonrecur.Normal-Primary.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC.Nonrecur.Normal-Primary.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Both.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-CIS.Both.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Both.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-CIS.Both.Down.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Both.Up.MSigDB.tsv Output/RSEM/Enrichment/Bowtie2.SQC-Primary.Both.Down.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Both.Up.MSigDB.tsv Output/RSEM/Enrichment/STAR.SQC-Primary.Both.Down.MSigDB.tsv
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
	$(DOCKER) $(BASH) $(addprefix /,$< Output/PyClone/BWA/$*.tsv Output/PyClone/BWA/$(PAIRED_PRIMARY).tsv Output/PyClone/BWA/$*.Cellularity.txt Output/PyClone/BWA/$(PAIRED_PRIMARY).Cellularity.txt $(@D)) 1> $(@D).stdout 2> $(@D).stderr || true

# I know that this recipe contains unnecessary dependencies; but I cannot see the solution.
Output/PyClone/Bowtie2/%/loci.tsv: Bash/run_pyclone.bash $(addprefix Output/PyClone/Bowtie2/,$(addsuffix .tsv,$(WES_ID))) $(addprefix Output/PyClone/Bowtie2/,$(addsuffix .Cellularity.txt,$(WES_ID))) | build.log
	$(DOCKER) $(BASH) $(addprefix /,$< Output/PyClone/Bowtie2/$*.tsv Output/PyClone/Bowtie2/$(PAIRED_PRIMARY).tsv Output/PyClone/Bowtie2/$*.Cellularity.txt Output/PyClone/Bowtie2/$(PAIRED_PRIMARY).Cellularity.txt $(@D)) 1> $(@D).stdout 2> $(@D).stderr || true

Output/PyClone/BWA-Plots/%.CCF.pdf: Python/draw_output_pyclone.py Output/PyClone/BWA/%/loci.tsv $(GENE_CENCUS_ONE_CSV) $(KNOWN_GENE_GTF) | build.log Output/PyClone/BWA-Plots
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --CCF 1> $@.stdout 2> $@.stderr

Output/PyClone/BWA-Plots/%.VAF.pdf: Python/draw_output_pyclone.py Output/PyClone/BWA/%/loci.tsv $(GENE_CENCUS_ONE_CSV) $(KNOWN_GENE_GTF) | build.log Output/PyClone/BWA-Plots
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --VAF 1> $@.stdout 2> $@.stderr

Output/PyClone/Bowtie2-Plots/%.CCF.pdf: Python/draw_output_pyclone.py Output/PyClone/Bowtie2/%/loci.tsv $(GENE_CENCUS_ONE_CSV) $(KNOWN_GENE_GTF) | build.log Output/PyClone/Bowtie2-Plots
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --CCF 1> $@.stdout 2> $@.stderr

Output/PyClone/Bowtie2-Plots/%.VAF.pdf: Python/draw_output_pyclone.py Output/PyClone/Bowtie2/%/loci.tsv $(GENE_CENCUS_ONE_CSV) $(KNOWN_GENE_GTF) | build.log Output/PyClone/Bowtie2-Plots
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

Output/CIBERSORTx/GSE131907_signature.tsv: Python/make_signature_cibersort.py $(LUNG_GSE131907_ANNOTATION_TXT) $(LUNG_GSE131907_NORMALIZED_MATRIX_TXT) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/CIBERSORTx/%/CIBERSORTxGEP_NA_Weights.txt: Output/CIBERSORTx/%.tsv Output/CIBERSORTx/GSE131907_signature.tsv | build.log
	$(DOCKER) rm -rfv $(@D) && mkdir -p $(@D)
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

Output/VAF/BWA.SQC.tar: Python/draw_vaf_plots.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_ONE_CSV) | Output/VAF build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/VAF/BWA.ADC.tar: Python/draw_vaf_plots.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_ONE_CSV) | Output/VAF build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/VAF/Bowtie2.SQC.tar: Python/draw_vaf_plots.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_ONE_CSV) | Output/VAF build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/VAF/Bowtie2.ADC.tar: Python/draw_vaf_plots.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_ONE_CSV) | Output/VAF build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/VAF-tSNE/BWA.SQC.pdf: Python/draw_vaf_tsne.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/VAF-tSNE build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/VAF-tSNE/BWA.ADC.pdf: Python/draw_vaf_tsne.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_ONE_CSV) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt | Output/VAF-tSNE build.log
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

Output/Revolver/BWA/SQC.tsv: Python/make_input_revolver.py $(addprefix Output/PyClone/BWA/,$(addsuffix /loci.tsv,$(filter-out %P,$(filter-out %P1,$(WES_TUMOR_ID))))) $(CLINICAL_DATA_CSV) Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) | build.log Output/Revolver/BWA
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Revolver/BWA/ADC.tsv: Python/make_input_revolver.py $(addprefix Output/PyClone/BWA/,$(addsuffix /loci.tsv,$(filter-out %P,$(filter-out %P1,$(WES_TUMOR_ID))))) $(CLINICAL_DATA_CSV) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) | build.log Output/Revolver/BWA
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Revolver/Bowtie2/SQC.tsv: Python/make_input_revolver.py $(addprefix Output/PyClone/Bowtie2/,$(addsuffix /loci.tsv,$(filter-out %P,$(filter-out %P1,$(WES_TUMOR_ID))))) $(CLINICAL_DATA_CSV) Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) | build.log Output/Revolver/Bowtie2
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Revolver/Bowtie2/ADC.tsv: Python/make_input_revolver.py $(addprefix Output/PyClone/Bowtie2/,$(addsuffix /loci.tsv,$(filter-out %P,$(filter-out %P1,$(WES_TUMOR_ID))))) $(CLINICAL_DATA_CSV) Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) | build.log Output/Revolver/Bowtie2
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
	$(WORKDIR_DOCKER) $(BASH) /Tools/arriba_v2.1.0/download_references.sh "hg38+RefSeq" 1> $@.stdout 2> $@.stderr

Output/Arriba/Raw/%/fusions.tsv: $(GENCODE_GTF) $(REFERENCE_FASTA) Data/WTS/%_*.fastq.gz | Output/Arriba/Raw build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D)) && mkdir -p $(@D)
	$(WORKDIR_DOCKER) run_arriba.sh $(addprefix /,$(ARRIBA_STAR_DIR) $(filter %.gtf,$^) $(filter %.fasta,$^) $(ARRIBA_BLACKLIST_TSV) $(ARRIBA_KNOWN_FUSIONS_TSV) $(ARRIBA_PROTEIN_DOMAIN_GFF3)) $(CPUS) $(addprefix /,$(sort $(filter %.fastq.gz,$^))) 1> $(@D).stdout 2> $(@D).stderr

# I know that this recipe contains unnecessary dependencies; but I cannot see the solution.
# It is working even if there is no matched normal sample.
Output/Arriba/Filtered/%.tsv: Python/filter_arriba.py Output/Arriba/Raw/%/fusions.tsv $(addprefix Output/Arriba/Raw/,$(addsuffix /fusions.tsv,$(WTS_NORMAL_ID))) | build.log Output/Arriba/Filtered
	$(DOCKER) $(PYTHON) $(addprefix /,$< Output/Arriba/Raw/$*/fusions.tsv Output/Arriba/Raw/$(PAIRED_NORMAL)/fusions.tsv $@) 1> $@.stdout 2> $@.stderr

Output/Arriba/ADC.pdf: Python/aggregate_arriba.py $(addprefix Output/Arriba/Filtered/,$(addsuffix .tsv,$(WTS_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_ONE_CSV) | Output/Arriba build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Arriba/SQC.pdf: Python/aggregate_arriba.py $(addprefix Output/Arriba/Filtered/,$(addsuffix .tsv,$(WTS_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_ONE_CSV) | Output/Arriba build.log
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
	$(DOCKER) $(BASH) $(addprefix /,$^ $(@D)) 1> $(@D).stdout 2> $(@D).stderr

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

Output/FusionCatcher/Plot/ADC.STAR.heatmap.pdf: Python/aggregate_fusioncatcher.py $(addprefix Output/FusionCatcher/STAR/,$(addsuffix .tsv,$(WTS_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/FusionCatcher/Plot build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/FusionCatcher/Plot/SQC.STAR.heatmap.pdf: Python/aggregate_fusioncatcher.py $(addprefix Output/FusionCatcher/STAR/,$(addsuffix .tsv,$(WTS_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/FusionCatcher/Plot build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/FusionCatcher/Plot/ADC.Bowtie2.heatmap.pdf: Python/aggregate_fusioncatcher.py $(addprefix Output/FusionCatcher/Bowtie2/,$(addsuffix .tsv,$(WTS_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/FusionCatcher/Plot build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/FusionCatcher/Plot/SQC.Bowtie2.heatmap.pdf: Python/aggregate_fusioncatcher.py $(addprefix Output/FusionCatcher/Bowtie2/,$(addsuffix .tsv,$(WTS_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/FusionCatcher/Plot build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/FusionCatcher/Plot/ADC.STAR.circos.pdf: Python/aggregate_fusioncatcher_circos.py $(addprefix Output/FusionCatcher/STAR/,$(addsuffix .tsv,$(WTS_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | Output/FusionCatcher/Plot build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/FusionCatcher/Plot/SQC.STAR.circos.pdf: Python/aggregate_fusioncatcher_circos.py $(addprefix Output/FusionCatcher/STAR/,$(addsuffix .tsv,$(WTS_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV)| Output/FusionCatcher/Plot build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/FusionCatcher/Plot/ADC.Bowtie2.circos.pdf: Python/aggregate_fusioncatcher_circos.py $(addprefix Output/FusionCatcher/Bowtie2/,$(addsuffix .tsv,$(WTS_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | Output/FusionCatcher/Plot build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/FusionCatcher/Plot/SQC.Bowtie2.circos.pdf: Python/aggregate_fusioncatcher_circos.py $(addprefix Output/FusionCatcher/Bowtie2/,$(addsuffix .tsv,$(WTS_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | Output/FusionCatcher/Plot build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

step26-prepar: Output/FusionCatcher/Reference/human_v102/genes.fa
step26-STAR: $(addprefix Output/FusionCatcher/STAR/,$(addsuffix /normal.tsv,$(WTS_TUMOR_ID))) $(addprefix Output/FusionCatcher/STAR/,$(addsuffix .tsv,$(WTS_TUMOR_ID)))
step26-Bowtie: $(addprefix Output/FusionCatcher/Bowtie2/,$(addsuffix /normal.tsv,$(WTS_TUMOR_ID))) $(addprefix Output/FusionCatcher/Bowtie2/,$(addsuffix .tsv,$(WTS_TUMOR_ID)))
step26-heatmap: Output/FusionCatcher/Plot/ADC.STAR.heatmap.pdf Output/FusionCatcher/Plot/SQC.STAR.heatmap.pdf Output/FusionCatcher/Plot/ADC.Bowtie2.heatmap.pdf Output/FusionCatcher/Plot/SQC.Bowtie2.heatmap.pdf
step26-circos: Output/FusionCatcher/Plot/ADC.STAR.circos.pdf Output/FusionCatcher/Plot/SQC.STAR.circos.pdf Output/FusionCatcher/Plot/ADC.Bowtie2.circos.pdf Output/FusionCatcher/Plot/SQC.Bowtie2.circos.pdf
step26: step26-prepare step26-STAR step26-Bowtie step26-heatmap step26-circos
.PHONY: step26-prepare step26-STAR step26-Bowtie step26-heatmap step26-circos step26

# Step 27 (STAR-Fusion)
Output/STAR-Fusion:
	mkdir -p $@

Output/STAR-Fusion/%/star-fusion.fusion_predictions.tsv Output/STAR-Fusion/%/star-fusion.fusion_predictions.abridged.tsv: Data/WTS/%_*.fastq.gz | Output/STAR-Fusion build.log
	$(DOCKER) rm -rf $(addsuffix /,$(@D)) && mkdir -p $(@D)
	$(DOCKER) STAR-Fusion --genome_lib_dir $(addprefix /,$(CTAT_RESOURCE_LIB_DIR)) --left_fq $(addprefix /,$(word 1,$(sort $(filter %.fastq.gz,$^)))) --right_fq $(addprefix /,$(word 2,$(sort $(filter %.fastq.gz,$^)))) --CPU $(CPUS) --output_dir $(addprefix /,$(@D)) 1> $(@D).stdout 2> $(@D).stderr

step27: $(addprefix Output/STAR-Fusion/,$(addsuffix /star-fusion.fusion_predictions.tsv,$(WTS_TUMOR_ID)))
.PHONY: step27

# Step 28 (MuSiC)
Output/Deconvolution Output/MuSiC:
	mkdir -p $@

Output/Deconvolution/Bowtie2.SQC.RSEM.tsv: Python/aggregate_RSEM_deconvolution.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(WTS_ID))) $(CLINICAL_DATA_CSV) $(GENCODE_GTF) $(TREMBL_TSV) | Output/Deconvolution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution/Bowtie2.ADC.RSEM.tsv: Python/aggregate_RSEM_deconvolution.py $(addprefix Output/RSEM/Bowtie2/,$(addsuffix .genes.results,$(WTS_ID))) $(CLINICAL_DATA_CSV) $(GENCODE_GTF) $(TREMBL_TSV) | Output/Deconvolution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Deconvolution/STAR.SQC.RSEM.tsv: Python/aggregate_RSEM_deconvolution.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(WTS_ID))) $(CLINICAL_DATA_CSV) $(GENCODE_GTF) $(TREMBL_TSV) | Output/Deconvolution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution/STAR.ADC.RSEM.tsv: Python/aggregate_RSEM_deconvolution.py $(addprefix Output/RSEM/STAR/,$(addsuffix .genes.results,$(WTS_ID))) $(CLINICAL_DATA_CSV) $(GENCODE_GTF) $(TREMBL_TSV) | Output/Deconvolution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/MuSiC/%.GSE131907.tsv: R/run_MuSiC.R Output/Deconvolution/%.RSEM.tsv $(LUNG_GSE131907_ANNOTATION_TXT) $(LUNG_GSE131907_RAW_MATRIX_RDS) | build.log Output/MuSiC
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.tsv,$^)) --annotation $(addprefix /,$(filter %.txt,$^)) --reference $(addprefix /,$(filter %.rds,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/MuSiC/%.GSE162498.tsv: R/run_MuSiC_2.R Output/Deconvolution/%.RSEM.tsv $(LUNG_GSE162498_11_TUMORS_RDS) | build.log Output/MuSiC
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.tsv,$^)) --reference $(addprefix /,$(filter %.Rds,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/MuSiC/%.heatmap.pdf: Python/draw_heatmap_MuSiC.py Output/MuSiC/%.tsv | build
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/MuSiC/%.violin.tar: Python/draw_violin_plot_MuSiC.py Output/MuSiC/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MuSiC/%.cluster.pdf: Python/draw_clustermap_MuSiC.py Output/MuSiC/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/MuSiC/%.Tcell_density.pdf: Python/draw_Tcell_density_violin_MuSiC.py Output/MuSiC/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step28-prepare: Output/Deconvolution/Bowtie2.SQC.RSEM.tsv Output/Deconvolution/Bowtie2.ADC.RSEM.tsv Output/Deconvolution/STAR.SQC.RSEM.tsv Output/Deconvolution/STAR.ADC.RSEM.tsv
step28-MuSiC: Output/MuSiC/Bowtie2.SQC.GSE131907.tsv Output/MuSiC/Bowtie2.ADC.GSE131907.tsv Output/MuSiC/STAR.SQC.GSE131907.tsv Output/MuSiC/STAR.ADC.GSE131907.tsv Output/MuSiC/Bowtie2.SQC.GSE162498.tsv Output/MuSiC/Bowtie2.ADC.GSE162498.tsv Output/MuSiC/STAR.SQC.GSE162498.tsv Output/MuSiC/STAR.ADC.GSE162498.tsv
step28-heatmap: Output/MuSiC/Bowtie2.SQC.GSE131907.heatmap.pdf Output/MuSiC/Bowtie2.ADC.GSE131907.heatmap.pdf Output/MuSiC/STAR.SQC.GSE131907.heatmap.pdf Output/MuSiC/STAR.ADC.GSE131907.heatmap.pdf Output/MuSiC/Bowtie2.SQC.GSE162498.heatmap.pdf Output/MuSiC/Bowtie2.ADC.GSE162498.heatmap.pdf Output/MuSiC/STAR.SQC.GSE162498.heatmap.pdf Output/MuSiC/STAR.ADC.GSE162498.heatmap.pdf
step28-violin: Output/MuSiC/Bowtie2.SQC.GSE131907.violin.tar Output/MuSiC/Bowtie2.ADC.GSE131907.violin.tar Output/MuSiC/STAR.SQC.GSE131907.violin.tar Output/MuSiC/STAR.ADC.GSE131907.violin.tar Output/MuSiC/Bowtie2.SQC.GSE162498.violin.tar Output/MuSiC/Bowtie2.ADC.GSE162498.violin.tar Output/MuSiC/STAR.SQC.GSE162498.violin.tar Output/MuSiC/STAR.ADC.GSE162498.violin.tar
step28-Tcell: Output/MuSiC/STAR.SQC.GSE131907.Tcell_density.pdf Output/MuSiC/STAR.ADC.GSE131907.Tcell_density.pdf Output/MuSiC/Bowtie2.SQC.GSE131907.Tcell_density.pdf Output/MuSiC/Bowtie2.ADC.GSE131907.Tcell_density.pdf
step28: step28-prepare step28-MuSiC step28-heatmap step28-violin step28-cluster step28-Tcell
.PHONY: step28-prepare step28-MuSiC step28-heatmap step28-violin step28-cluster step28-Tcell step28

# Step 29 (BisqueRNA)
Output/BisqueRNA:
	mkdir -p $@

Output/BisqueRNA/%.GSE131907.tsv: R/run_BisqueRNA.R Output/Deconvolution/%.RSEM.tsv $(LUNG_GSE131907_ANNOTATION_TXT) $(LUNG_GSE131907_RAW_MATRIX_RDS) | Output/BisqueRNA build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.tsv,$^)) --annotation $(addprefix /,$(filter %.txt,$^)) --reference $(addprefix /,$(filter %.rds,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/BisqueRNA/%.GSE162498.tsv: R/run_BisqueRNA_2.R Output/Deconvolution/%.RSEM.tsv $(LUNG_GSE162498_11_TUMORS_RDS) | build.log Output/BisqueRNA
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.tsv,$^)) --reference $(addprefix /,$(filter %.Rds,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/BisqueRNA/%.heatmap.pdf: Python/draw_heatmap_BisqueRNA.py Output/BisqueRNA/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/BisqueRNA/%.violin.tar: Python/draw_violin_plot_BisqueRNA.py Output/BisqueRNA/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/BisqueRNA/%.cluster.pdf: Python/draw_clustermap_BisqueRNA.py Output/BisqueRNA/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/BisqueRNA/%.Tcell_density.pdf: Python/draw_Tcell_density_violin_BisqueRNA.py Output/BisqueRNA/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step29-BisqueRNA: Output/BisqueRNA/Bowtie2.SQC.GSE131907.tsv Output/BisqueRNA/Bowtie2.ADC.GSE131907.tsv Output/BisqueRNA/STAR.SQC.GSE131907.tsv Output/BisqueRNA/STAR.ADC.GSE131907.tsv Output/BisqueRNA/Bowtie2.SQC.GSE162498.tsv Output/BisqueRNA/Bowtie2.ADC.GSE162498.tsv Output/BisqueRNA/STAR.SQC.GSE162498.tsv Output/BisqueRNA/STAR.ADC.GSE162498.tsv
step29-heatmap: Output/BisqueRNA/Bowtie2.SQC.GSE131907.heatmap.pdf Output/BisqueRNA/Bowtie2.ADC.GSE131907.heatmap.pdf Output/BisqueRNA/STAR.SQC.GSE131907.heatmap.pdf Output/BisqueRNA/STAR.ADC.GSE131907.heatmap.pdf Output/BisqueRNA/Bowtie2.SQC.GSE162498.heatmap.pdf Output/BisqueRNA/Bowtie2.ADC.GSE162498.heatmap.pdf Output/BisqueRNA/STAR.SQC.GSE162498.heatmap.pdf Output/BisqueRNA/STAR.ADC.GSE162498.heatmap.pdf
step29-violin: Output/BisqueRNA/Bowtie2.SQC.GSE131907.violin.tar Output/BisqueRNA/Bowtie2.ADC.GSE131907.violin.tar Output/BisqueRNA/STAR.SQC.GSE131907.violin.tar Output/BisqueRNA/STAR.ADC.GSE131907.violin.tar Output/BisqueRNA/Bowtie2.SQC.GSE162498.violin.tar Output/BisqueRNA/Bowtie2.ADC.GSE162498.violin.tar Output/BisqueRNA/STAR.SQC.GSE162498.violin.tar Output/BisqueRNA/STAR.ADC.GSE162498.violin.tar
step29-cluster: Output/BisqueRNA/Bowtie2.SQC.GSE131907.cluster.pdf Output/BisqueRNA/Bowtie2.ADC.GSE131907.cluster.pdf Output/BisqueRNA/STAR.SQC.GSE131907.cluster.pdf Output/BisqueRNA/STAR.ADC.GSE131907.cluster.pdf Output/BisqueRNA/Bowtie2.SQC.GSE162498.cluster.pdf Output/BisqueRNA/Bowtie2.ADC.GSE162498.cluster.pdf Output/BisqueRNA/STAR.SQC.GSE162498.cluster.pdf Output/BisqueRNA/STAR.ADC.GSE162498.cluster.pdf
step29-Tcell: Output/BisqueRNA/STAR.SQC.GSE131907.Tcell_density.pdf Output/BisqueRNA/STAR.ADC.GSE131907.Tcell_density.pdf Output/BisqueRNA/Bowtie2.SQC.GSE131907.Tcell_density.pdf Output/BisqueRNA/Bowtie2.ADC.GSE131907.Tcell_density.pdf
step29: step29-BisqueRNA step29-heatmap step29-violin step29-cluster step29-Tcell
.PHONY: step29-BisqueRNA step29-heatmap step29-violin step29-cluster step29-Tcell step29

# Step 30 (SCDC)
Output/SCDC:
	mkdir -p $@

Output/SCDC/%.GSE131907.tsv: R/run_SCDC.R Output/Deconvolution/%.RSEM.tsv $(LUNG_GSE131907_ANNOTATION_TXT) $(LUNG_GSE131907_RAW_MATRIX_RDS) | Output/SCDC build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.tsv,$^)) --annotation $(addprefix /,$(filter %.txt,$^)) --reference $(addprefix /,$(filter %.rds,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/SCDC/%.GSE162498.tsv: R/run_SCDC_2.R Output/Deconvolution/%.RSEM.tsv $(LUNG_GSE162498_11_TUMORS_RDS) | build.log Output/SCDC
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.tsv,$^)) --reference $(addprefix /,$(filter %.Rds,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/SCDC/%.heatmap.pdf: Python/draw_heatmap_MuSiC.py Output/SCDC/%.tsv | build
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/SCDC/%.violin.tar: Python/draw_violin_plot_MuSiC.py Output/SCDC/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/SCDC/%.cluster.pdf: Python/draw_clustermap_MuSiC.py Output/SCDC/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/SCDC/%.Tcell_density.pdf: Python/draw_Tcell_density_violin_MuSiC.py Output/SCDC/%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step30-SCDC: Output/SCDC/Bowtie2.SQC.GSE131907.tsv Output/SCDC/Bowtie2.ADC.GSE131907.tsv Output/SCDC/STAR.SQC.GSE131907.tsv Output/SCDC/STAR.ADC.GSE131907.tsv Output/SCDC/Bowtie2.SQC.GSE162498.tsv Output/SCDC/Bowtie2.ADC.GSE162498.tsv Output/SCDC/STAR.SQC.GSE162498.tsv Output/SCDC/STAR.ADC.GSE162498.tsv
step30-heatmap: Output/SCDC/Bowtie2.SQC.GSE131907.heatmap.pdf Output/SCDC/Bowtie2.ADC.GSE131907.heatmap.pdf Output/SCDC/STAR.SQC.GSE131907.heatmap.pdf Output/SCDC/STAR.ADC.GSE131907.heatmap.pdf Output/SCDC/Bowtie2.SQC.GSE162498.heatmap.pdf Output/SCDC/Bowtie2.ADC.GSE162498.heatmap.pdf Output/SCDC/STAR.SQC.GSE162498.heatmap.pdf Output/SCDC/STAR.ADC.GSE162498.heatmap.pdf
step30-violin: Output/SCDC/Bowtie2.SQC.GSE131907.violin.tar Output/SCDC/Bowtie2.ADC.GSE131907.violin.tar Output/SCDC/STAR.SQC.GSE131907.violin.tar Output/SCDC/STAR.ADC.GSE131907.violin.tar Output/SCDC/Bowtie2.SQC.GSE162498.violin.tar Output/SCDC/Bowtie2.ADC.GSE162498.violin.tar Output/SCDC/STAR.SQC.GSE162498.violin.tar Output/SCDC/STAR.ADC.GSE162498.violin.tar
step30-cluster: Output/SCDC/Bowtie2.SQC.GSE131907.cluster.pdf Output/SCDC/Bowtie2.ADC.GSE131907.cluster.pdf Output/SCDC/STAR.SQC.GSE131907.cluster.pdf Output/SCDC/STAR.ADC.GSE131907.cluster.pdf Output/SCDC/Bowtie2.SQC.GSE162498.cluster.pdf Output/SCDC/Bowtie2.ADC.GSE162498.cluster.pdf Output/SCDC/STAR.SQC.GSE162498.cluster.pdf Output/SCDC/STAR.ADC.GSE162498.cluster.pdf
step30-Tcell: Output/SCDC/STAR.SQC.GSE131907.Tcell_density.pdf Output/SCDC/STAR.ADC.GSE131907.Tcell_density.pdf Output/SCDC/Bowtie2.SQC.GSE131907.Tcell_density.pdf Output/SCDC/Bowtie2.ADC.GSE131907.Tcell_density.pdf
step30: step30-SCDC step30-heatmap step30-violin step30-cluster step30-Tcell
.PHONY: step30-SCDC step30-heatmap step30-violin step30-cluster step30-Tcell step30

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
	$(DOCKER) $(BASH) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step31-prepare: Output/TIDE/Bowtie2.SQC.tsv Output/TIDE/Bowtie2.ADC.tsv Output/TIDE/STAR.SQC.tsv Output/TIDE/STAR.ADC.tsv
step31-tide: Output/TIDE/Bowtie2.SQC.TIDE.tsv Output/TIDE/Bowtie2.ADC.TIDE.tsv Output/TIDE/STAR.SQC.TIDE.tsv Output/TIDE/STAR.ADC.TIDE.tsv
step31: step31-prepare step31-tide
.PHONY: step31-prepare step31-tide step31

# Step 32 (Data Composition)
Output/Data_Composition:
	mkdir -p $@

Output/Data_Composition/WES_samples.tex: Python/get_data_composition.py $(CLINICAL_DATA_CSV) | build.log Output/Data_Composition
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^)) $(WES_ID) $(addprefix /,$(filter %.csv,$^) $@) 1> $@.stdout 2> $@.stderr

Output/Data_Composition/WTS_samples.tex: Python/get_data_composition.py $(CLINICAL_DATA_CSV) | build.log Output/Data_Composition
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^)) $(WTS_ID) $(addprefix /,$(filter %.csv,$^) $@) 1> $@.stdout 2> $@.stderr

Output/Data_Composition/WES_SQC_recurrence.tex: Python/get_data_composition_recurrence.py $(CLINICAL_DATA_CSV) | build.log Output/Data_Composition
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^)) $(WES_ID) $(addprefix /,$(filter %.csv,$^) $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Data_Composition/WES_ADC_recurrence.tex: Python/get_data_composition_recurrence.py $(CLINICAL_DATA_CSV) | build.log Output/Data_Composition
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^)) $(WES_ID) $(addprefix /,$(filter %.csv,$^) $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Data_Composition/WTS_SQC_recurrence.tex: Python/get_data_composition_recurrence.py $(CLINICAL_DATA_CSV) | build.log Output/Data_Composition
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^)) $(WTS_ID) $(addprefix /,$(filter %.csv,$^) $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Data_Composition/WTS_ADC_recurrence.tex: Python/get_data_composition_recurrence.py $(CLINICAL_DATA_CSV) | build.log Output/Data_Composition
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^)) $(WTS_ID) $(addprefix /,$(filter %.csv,$^) $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Data_Composition/WES_SQC_smoking.tex: Python/get_data_composition_smoking.py $(CLINICAL_DATA_CSV) | build.log Output/Data_Composition
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^)) $(WES_ID) $(addprefix /,$(filter %.csv,$^) $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Data_Composition/WES_ADC_smoking.tex: Python/get_data_composition_smoking.py $(CLINICAL_DATA_CSV) | build.log Output/Data_Composition
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^)) $(WES_ID) $(addprefix /,$(filter %.csv,$^) $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Data_Composition/WTS_SQC_smoking.tex: Python/get_data_composition_smoking.py $(CLINICAL_DATA_CSV) | build.log Output/Data_Composition
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^)) $(WTS_ID) $(addprefix /,$(filter %.csv,$^) $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Data_Composition/WTS_ADC_smoking.tex: Python/get_data_composition_smoking.py $(CLINICAL_DATA_CSV) | build.log Output/Data_Composition
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^)) $(WTS_ID) $(addprefix /,$(filter %.csv,$^) $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Data_Composition/WES_SQC_%_BWA.tex: Python/get_data_composition_gene.py $(CLINICAL_DATA_CSV) $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) | build.log Output/Data_Composition
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^)) $(WES_ID) $(addprefix /,$(filter %.csv,$^) $@) --mutect $(addprefix /,$(filter %.maf,$^)) --gene "$*" --SQC 1> $@.stdout 2> $@.stderr

Output/Data_Composition/WES_SQC_%_Bowtie2.tex: Python/get_data_composition_gene.py $(CLINICAL_DATA_CSV) $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) | build.log Output/Data_Composition
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^)) $(WES_ID) $(addprefix /,$(filter %.csv,$^) $@) --mutect $(addprefix /,$(filter %.maf,$^)) --gene "$*" --SQC 1> $@.stdout 2> $@.stderr

Output/Data_Composition/WES_ADC_%_BWA.tex: Python/get_data_composition_gene.py $(CLINICAL_DATA_CSV) $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) | build.log Output/Data_Composition
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^)) $(WES_ID) $(addprefix /,$(filter %.csv,$^) $@) --mutect $(addprefix /,$(filter %.maf,$^)) --gene "$*" --ADC 1> $@.stdout 2> $@.stderr

Output/Data_Composition/WES_ADC_%_Bowtie2.tex: Python/get_data_composition_gene.py $(CLINICAL_DATA_CSV) $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) | build.log Output/Data_Composition
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter %.py,$^)) $(WES_ID) $(addprefix /,$(filter %.csv,$^) $@) --mutect $(addprefix /,$(filter %.maf,$^)) --gene "$*" --ADC 1> $@.stdout 2> $@.stderr

step32-all: Output/Data_Composition/WES_samples.tex Output/Data_Composition/WTS_samples.tex
step32-Recur: Output/Data_Composition/WES_SQC_recurrence.tex Output/Data_Composition/WES_ADC_recurrence.tex Output/Data_Composition/WTS_SQC_recurrence.tex Output/Data_Composition/WTS_ADC_recurrence.tex
step32-Smoking: Output/Data_Composition/WES_SQC_smoking.tex Output/Data_Composition/WES_ADC_smoking.tex Output/Data_Composition/WTS_SQC_smoking.tex Output/Data_Composition/WTS_ADC_smoking.tex
step32-gene: Output/Data_Composition/WES_SQC_EGFR_BWA.tex Output/Data_Composition/WES_SQC_EGFR_Bowtie2.tex Output/Data_Composition/WES_ADC_EGFR_BWA.tex Output/Data_Composition/WES_ADC_EGFR_Bowtie2.tex Output/Data_Composition/WES_SQC_KRAS_BWA.tex Output/Data_Composition/WES_SQC_KRAS_Bowtie2.tex Output/Data_Composition/WES_ADC_KRAS_BWA.tex Output/Data_Composition/WES_ADC_KRAS_Bowtie2.tex Output/Data_Composition/WES_SQC_TP53_BWA.tex Output/Data_Composition/WES_SQC_TP53_Bowtie2.tex Output/Data_Composition/WES_ADC_TP53_BWA.tex Output/Data_Composition/WES_ADC_TP53_Bowtie2.tex
step32: step32-all step32-Recur step32-Smoking step32-gene
.PHONY: step32-all step32-Recur step32-Smoking step32-gene step32

# Step 33 (SigProfilerExtractor)
Output/SigProfiler/BWA Output/SigProfiler/Bowtie2 Output/SigProfiler/BWA-Violin Output/SigProfiler/Bowtie2-Violin Output/SigProfiler/BWA-Bar Output/SigProfiler/Bowtie2-Bar Output/SigProfiler/BWA-Distribution Output/SigProfiler/Bowtie2-Distribution:
	mkdir -p $@

Output/SigProfiler/BWA/COSMIC_DBS78_Activities_refit.txt Output/SigProfiler/BWA/DBS78_De-Novo_Activities_refit.txt Output/SigProfiler/BWA/COSMIC_ID83_Activities_refit.txt Output/SigProfiler/BWA/ID83_De-Novo_Activities_refit.txt Output/SigProfiler/BWA/COSMIC_SBS96_Activities_refit.txt Output/SigProfiler/BWA/SBS96_De-Novo_Activities_refit.txt: Python/run_SigProfiler.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.vcf,$(WES_TUMOR_ID))) | Output/SigProfiler/BWA build.log
	$(DOCKER) rm -rfv $(addprefix /,$(@D))
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(@D)) --cpus $(CPUS) 1> $(@D).stdout 2> $(@D).stderr
	cp -v $(@D)/**/*Activities_refit.txt $(@D)

Output/SigProfiler/Bowtie2/COSMIC_DBS78_Activities_refit.txt Output/SigProfiler/Bowtie2/DBS78_De-Novo_Activities_refit.txt Output/SigProfiler/Bowtie2/COSMIC_ID83_Activities_refit.txt Output/SigProfiler/Bowtie2/ID83_De-Novo_Activities_refit.txt Output/SigProfiler/Bowtie2/COSMIC_SBS96_Activities_refit.txt Output/SigProfiler/Bowtie2/SBS96_De-Novo_Activities_refit.txt: Python/run_SigProfiler.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.vcf,$(WES_TUMOR_ID))) | Output/SigProfiler/Bowtie2 build.log
	$(DOCKER) rm -rfv $(addprefix /,$(@D))
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(@D)) --cpus $(CPUS) 1> $(@D).stdout 2> $(@D).stderr
	cp -v $(@D)/**/*Activities_refit.txt $(@D)

Output/SigProfiler/%-Violin/COSMIC-SingleBase.SQC.Subtype.Absolute.tar: Python/plot_signature.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-SingleBase.ADC.Subtype.Absolute.tar: Python/plot_signature.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-DoubleBase.SQC.Subtype.Absolute.tar: Python/plot_signature.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-DoubleBase.ADC.Subtype.Absolute.tar: Python/plot_signature.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-Indel.SQC.Subtype.Absolute.tar: Python/plot_signature.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-Indel.ADC.Subtype.Absolute.tar: Python/plot_signature.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-SingleBase.SQC.Subtype.Relative.tar: Python/plot_signature_violin_relative.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-SingleBase.ADC.Subtype.Relative.tar: Python/plot_signature_violin_relative.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-DoubleBase.SQC.Subtype.Relative.tar: Python/plot_signature_violin_relative.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-DoubleBase.ADC.Subtype.Relative.tar: Python/plot_signature_violin_relative.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-Indel.SQC.Subtype.Relative.tar: Python/plot_signature_violin_relative.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-Indel.ADC.Subtype.Relative.tar: Python/plot_signature_violin_relative.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-SingleBase.SQC.Recurrence.Absolute.tar: Python/plot_signature_2.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-SingleBase.ADC.Recurrence.Absolute.tar: Python/plot_signature_2.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-DoubleBase.SQC.Recurrence.Absolute.tar: Python/plot_signature_2.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-DoubleBase.ADC.Recurrence.Absolute.tar: Python/plot_signature_2.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-Indel.SQC.Recurrence.Absolute.tar: Python/plot_signature_2.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-Indel.ADC.Recurrence.Absolute.tar: Python/plot_signature_2.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-SingleBase.SQC.Smoking.Absolute.tar: Python/plot_signature_2.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-SingleBase.ADC.Smoking.Absolute.tar: Python/plot_signature_2.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-DoubleBase.SQC.Smoking.Absolute.tar: Python/plot_signature_2.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-DoubleBase.ADC.Smoking.Absolute.tar: Python/plot_signature_2.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-Indel.SQC.Smoking.Absolute.tar: Python/plot_signature_2.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-Indel.ADC.Smoking.Absolute.tar: Python/plot_signature_2.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-SingleBase.SQC.Smoking.Relative.tar: Python/plot_signature_2_relative.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-SingleBase.ADC.Smoking.Relative.tar: Python/plot_signature_2_relative.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-DoubleBase.SQC.Smoking.Relative.tar: Python/plot_signature_2_relative.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-DoubleBase.ADC.Smoking.Relative.tar: Python/plot_signature_2_relative.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-Indel.SQC.Smoking.Relative.tar: Python/plot_signature_2_relative.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-Indel.ADC.Smoking.Relative.tar: Python/plot_signature_2_relative.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-SingleBase.SQC.Recurrence.Relative.tar: Python/plot_signature_2_relative.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-SingleBase.ADC.Recurrence.Relative.tar: Python/plot_signature_2_relative.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-DoubleBase.SQC.Recurrence.Relative.tar: Python/plot_signature_2_relative.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-DoubleBase.ADC.Recurrence.Relative.tar: Python/plot_signature_2_relative.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-Indel.SQC.Recurrence.Relative.tar: Python/plot_signature_2_relative.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Violin/COSMIC-Indel.ADC.Recurrence.Relative.tar: Python/plot_signature_2_relative.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Violin
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-SingleBase.SQC.Absolute.pdf: Python/plot_signature_bar.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-SingleBase.SQC.Relative.pdf: Python/plot_signature_bar.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-DoubleBase.SQC.Absolute.pdf: Python/plot_signature_bar.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-DoubleBase.SQC.Relative.pdf: Python/plot_signature_bar.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-Indel.SQC.Absolute.pdf: Python/plot_signature_bar.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-Indel.SQC.Relative.pdf: Python/plot_signature_bar.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-SingleBase.ADC.Absolute.pdf: Python/plot_signature_bar.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-SingleBase.ADC.Relative.pdf: Python/plot_signature_bar.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-DoubleBase.ADC.Absolute.pdf: Python/plot_signature_bar.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-DoubleBase.ADC.Relative.pdf: Python/plot_signature_bar.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-Indel.ADC.Absolute.pdf: Python/plot_signature_bar.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-Indel.ADC.Relative.pdf: Python/plot_signature_bar.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-SingleBase.SQC.Subtype.Absolute.pdf: Python/plot_signature_bar_subtype.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-SingleBase.SQC.Subtype.Relative.pdf: Python/plot_signature_bar_subtype.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-DoubleBase.SQC.Subtype.Absolute.pdf: Python/plot_signature_bar_subtype.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-DoubleBase.SQC.Subtype.Relative.pdf: Python/plot_signature_bar_subtype.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-Indel.SQC.Subtype.Absolute.pdf: Python/plot_signature_bar_subtype.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-Indel.SQC.Subtype.Relative.pdf: Python/plot_signature_bar_subtype.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-SingleBase.ADC.Subtype.Absolute.pdf: Python/plot_signature_bar_subtype.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-SingleBase.ADC.Subtype.Relative.pdf: Python/plot_signature_bar_subtype.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-DoubleBase.ADC.Subtype.Absolute.pdf: Python/plot_signature_bar_subtype.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-DoubleBase.ADC.Subtype.Relative.pdf: Python/plot_signature_bar_subtype.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-Indel.ADC.Subtype.Absolute.pdf: Python/plot_signature_bar_subtype.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-Indel.ADC.Subtype.Relative.pdf: Python/plot_signature_bar_subtype.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-SingleBase.SQC.Smoking.Absolute.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-SingleBase.SQC.Smoking.Relative.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-DoubleBase.SQC.Smoking.Absolute.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-DoubleBase.SQC.Smoking.Relative.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-Indel.SQC.Smoking.Absolute.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-Indel.SQC.Smoking.Relative.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-SingleBase.ADC.Smoking.Absolute.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-SingleBase.ADC.Smoking.Relative.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-DoubleBase.ADC.Smoking.Absolute.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-DoubleBase.ADC.Smoking.Relative.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-Indel.ADC.Smoking.Absolute.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-Indel.ADC.Smoking.Relative.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-SingleBase.SQC.Recurrence.Absolute.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-SingleBase.SQC.Recurrence.Relative.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-DoubleBase.SQC.Recurrence.Absolute.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-DoubleBase.SQC.Recurrence.Relative.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-Indel.SQC.Recurrence.Absolute.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-Indel.SQC.Recurrence.Relative.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-SingleBase.ADC.Recurrence.Absolute.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-SingleBase.ADC.Recurrence.Relative.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-DoubleBase.ADC.Recurrence.Absolute.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-DoubleBase.ADC.Recurrence.Relative.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-Indel.ADC.Recurrence.Absolute.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Bar/COSMIC-Indel.ADC.Recurrence.Relative.pdf: Python/plot_signature_bar_subtype_clinical.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/SigProfiler/%-Bar
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Distribution/COSMIC-SingleBase.SQC.Absolute.pdf: Python/plot_signature_distribution.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | Output/SigProfiler/%-Distribution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Distribution/COSMIC-SingleBase.SQC.Relative.pdf: Python/plot_signature_distribution_relative.py Output/SigProfiler/%/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | Output/SigProfiler/%-Distribution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Distribution/COSMIC-DoubleBase.SQC.Absolute.pdf: Python/plot_signature_distribution.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | Output/SigProfiler/%-Distribution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Distribution/COSMIC-DoubleBase.SQC.Relative.pdf: Python/plot_signature_distribution_relative.py Output/SigProfiler/%/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | Output/SigProfiler/%-Distribution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Distribution/COSMIC-Indel.SQC.Absolute.pdf: Python/plot_signature_distribution.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | Output/SigProfiler/%-Distribution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/SigProfiler/%-Distribution/COSMIC-Indel.SQC.Relative.pdf: Python/plot_signature_distribution_relative.py Output/SigProfiler/%/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | Output/SigProfiler/%-Distribution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

step33-SigProfiler: Output/SigProfiler/BWA/COSMIC_DBS78_Activities_refit.txt Output/SigProfiler/Bowtie2/COSMIC_DBS78_Activities_refit.txt
step33-Violin-Subtype: Output/SigProfiler/BWA-Violin/COSMIC-SingleBase.SQC.Subtype.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-SingleBase.ADC.Subtype.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-SingleBase.SQC.Subtype.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-SingleBase.ADC.Subtype.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-DoubleBase.SQC.Subtype.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-DoubleBase.ADC.Subtype.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-DoubleBase.SQC.Subtype.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-DoubleBase.ADC.Subtype.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-Indel.SQC.Subtype.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-Indel.ADC.Subtype.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-Indel.SQC.Subtype.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-Indel.ADC.Subtype.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-SingleBase.SQC.Subtype.Relative.tar Output/SigProfiler/BWA-Violin/COSMIC-SingleBase.ADC.Subtype.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-SingleBase.SQC.Subtype.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-SingleBase.ADC.Subtype.Relative.tar Output/SigProfiler/BWA-Violin/COSMIC-DoubleBase.SQC.Subtype.Relative.tar Output/SigProfiler/BWA-Violin/COSMIC-DoubleBase.ADC.Subtype.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-DoubleBase.SQC.Subtype.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-DoubleBase.ADC.Subtype.Relative.tar Output/SigProfiler/BWA-Violin/COSMIC-Indel.SQC.Subtype.Relative.tar Output/SigProfiler/BWA-Violin/COSMIC-Indel.ADC.Subtype.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-Indel.SQC.Subtype.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-Indel.ADC.Subtype.Relative.tar
step33-Violin-Recurrence: Output/SigProfiler/BWA-Violin/COSMIC-SingleBase.SQC.Recurrence.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-SingleBase.ADC.Recurrence.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-DoubleBase.SQC.Recurrence.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-DoubleBase.ADC.Recurrence.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-Indel.SQC.Recurrence.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-Indel.ADC.Recurrence.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-SingleBase.SQC.Recurrence.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-SingleBase.ADC.Recurrence.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-DoubleBase.SQC.Recurrence.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-DoubleBase.ADC.Recurrence.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-Indel.SQC.Recurrence.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-Indel.ADC.Recurrence.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-SingleBase.SQC.Recurrence.Relative.tar Output/SigProfiler/BWA-Violin/COSMIC-SingleBase.ADC.Recurrence.Relative.tar Output/SigProfiler/BWA-Violin/COSMIC-DoubleBase.SQC.Recurrence.Relative.tar Output/SigProfiler/BWA-Violin/COSMIC-DoubleBase.ADC.Recurrence.Relative.tar Output/SigProfiler/BWA-Violin/COSMIC-Indel.SQC.Recurrence.Relative.tar Output/SigProfiler/BWA-Violin/COSMIC-Indel.ADC.Recurrence.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-SingleBase.SQC.Recurrence.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-SingleBase.ADC.Recurrence.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-DoubleBase.SQC.Recurrence.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-DoubleBase.ADC.Recurrence.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-Indel.SQC.Recurrence.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-Indel.ADC.Recurrence.Relative.tar
step33-Violin-Smoking: Output/SigProfiler/BWA-Violin/COSMIC-SingleBase.SQC.Smoking.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-SingleBase.ADC.Smoking.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-DoubleBase.SQC.Smoking.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-DoubleBase.ADC.Smoking.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-Indel.SQC.Smoking.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-Indel.ADC.Smoking.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-SingleBase.SQC.Smoking.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-SingleBase.ADC.Smoking.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-DoubleBase.SQC.Smoking.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-DoubleBase.ADC.Smoking.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-Indel.SQC.Smoking.Absolute.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-Indel.ADC.Smoking.Absolute.tar Output/SigProfiler/BWA-Violin/COSMIC-SingleBase.SQC.Smoking.Relative.tar Output/SigProfiler/BWA-Violin/COSMIC-SingleBase.ADC.Smoking.Relative.tar Output/SigProfiler/BWA-Violin/COSMIC-DoubleBase.SQC.Smoking.Relative.tar Output/SigProfiler/BWA-Violin/COSMIC-DoubleBase.ADC.Smoking.Relative.tar Output/SigProfiler/BWA-Violin/COSMIC-Indel.SQC.Smoking.Relative.tar Output/SigProfiler/BWA-Violin/COSMIC-Indel.ADC.Smoking.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-SingleBase.SQC.Smoking.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-SingleBase.ADC.Smoking.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-DoubleBase.SQC.Smoking.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-DoubleBase.ADC.Smoking.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-Indel.SQC.Smoking.Relative.tar Output/SigProfiler/Bowtie2-Violin/COSMIC-Indel.ADC.Smoking.Relative.tar
step33-Bar-SQC: Output/SigProfiler/BWA-Bar/COSMIC-SingleBase.SQC.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-SingleBase.SQC.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-DoubleBase.SQC.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-DoubleBase.SQC.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-Indel.SQC.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-Indel.SQC.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-SingleBase.SQC.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-SingleBase.SQC.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-DoubleBase.SQC.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-DoubleBase.SQC.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-Indel.SQC.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-Indel.SQC.Relative.pdf
step33-Bar-ADC: Output/SigProfiler/BWA-Bar/COSMIC-SingleBase.ADC.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-SingleBase.ADC.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-DoubleBase.ADC.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-DoubleBase.ADC.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-Indel.ADC.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-Indel.ADC.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-SingleBase.ADC.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-SingleBase.ADC.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-DoubleBase.ADC.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-DoubleBase.ADC.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-Indel.ADC.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-Indel.ADC.Relative.pdf
step33-Bar-Subtype: Output/SigProfiler/BWA-Bar/COSMIC-SingleBase.SQC.Subtype.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-SingleBase.SQC.Subtype.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-DoubleBase.SQC.Subtype.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-DoubleBase.SQC.Subtype.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-Indel.SQC.Subtype.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-Indel.SQC.Subtype.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-SingleBase.ADC.Subtype.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-SingleBase.ADC.Subtype.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-DoubleBase.ADC.Subtype.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-DoubleBase.ADC.Subtype.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-Indel.ADC.Subtype.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-Indel.ADC.Subtype.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-SingleBase.SQC.Subtype.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-SingleBase.SQC.Subtype.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-DoubleBase.SQC.Subtype.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-DoubleBase.SQC.Subtype.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-Indel.SQC.Subtype.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-Indel.SQC.Subtype.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-SingleBase.ADC.Subtype.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-SingleBase.ADC.Subtype.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-DoubleBase.ADC.Subtype.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-DoubleBase.ADC.Subtype.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-Indel.ADC.Subtype.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-Indel.ADC.Subtype.Relative.pdf
step33-Bar-Smoking: Output/SigProfiler/BWA-Bar/COSMIC-SingleBase.SQC.Smoking.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-SingleBase.SQC.Smoking.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-DoubleBase.SQC.Smoking.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-DoubleBase.SQC.Smoking.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-Indel.SQC.Smoking.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-Indel.SQC.Smoking.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-SingleBase.ADC.Smoking.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-SingleBase.ADC.Smoking.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-DoubleBase.ADC.Smoking.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-DoubleBase.ADC.Smoking.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-Indel.ADC.Smoking.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-Indel.ADC.Smoking.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-SingleBase.SQC.Smoking.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-SingleBase.SQC.Smoking.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-DoubleBase.SQC.Smoking.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-DoubleBase.SQC.Smoking.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-Indel.SQC.Smoking.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-Indel.SQC.Smoking.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-SingleBase.ADC.Smoking.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-SingleBase.ADC.Smoking.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-DoubleBase.ADC.Smoking.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-DoubleBase.ADC.Smoking.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-Indel.ADC.Smoking.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-Indel.ADC.Smoking.Relative.pdf
step33-Bar-Recurrence: Output/SigProfiler/BWA-Bar/COSMIC-SingleBase.SQC.Recurrence.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-SingleBase.SQC.Recurrence.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-DoubleBase.SQC.Recurrence.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-DoubleBase.SQC.Recurrence.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-Indel.SQC.Recurrence.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-Indel.SQC.Recurrence.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-SingleBase.ADC.Recurrence.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-SingleBase.ADC.Recurrence.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-DoubleBase.ADC.Recurrence.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-DoubleBase.ADC.Recurrence.Relative.pdf Output/SigProfiler/BWA-Bar/COSMIC-Indel.ADC.Recurrence.Absolute.pdf Output/SigProfiler/BWA-Bar/COSMIC-Indel.ADC.Recurrence.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-SingleBase.SQC.Recurrence.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-SingleBase.SQC.Recurrence.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-DoubleBase.SQC.Recurrence.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-DoubleBase.SQC.Recurrence.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-Indel.SQC.Recurrence.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-Indel.SQC.Recurrence.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-SingleBase.ADC.Recurrence.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-SingleBase.ADC.Recurrence.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-DoubleBase.ADC.Recurrence.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-DoubleBase.ADC.Recurrence.Relative.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-Indel.ADC.Recurrence.Absolute.pdf Output/SigProfiler/Bowtie2-Bar/COSMIC-Indel.ADC.Recurrence.Relative.pdf
step33-Distribution: Output/SigProfiler/BWA-Distribution/COSMIC-SingleBase.SQC.Absolute.pdf Output/SigProfiler/BWA-Distribution/COSMIC-SingleBase.SQC.Relative.pdf Output/SigProfiler/BWA-Distribution/COSMIC-DoubleBase.SQC.Absolute.pdf Output/SigProfiler/BWA-Distribution/COSMIC-DoubleBase.SQC.Relative.pdf Output/SigProfiler/BWA-Distribution/COSMIC-Indel.SQC.Absolute.pdf Output/SigProfiler/BWA-Distribution/COSMIC-Indel.SQC.Relative.pdf Output/SigProfiler/Bowtie2-Distribution/COSMIC-SingleBase.SQC.Absolute.pdf Output/SigProfiler/Bowtie2-Distribution/COSMIC-SingleBase.SQC.Relative.pdf Output/SigProfiler/Bowtie2-Distribution/COSMIC-DoubleBase.SQC.Absolute.pdf Output/SigProfiler/Bowtie2-Distribution/COSMIC-DoubleBase.SQC.Relative.pdf Output/SigProfiler/Bowtie2-Distribution/COSMIC-Indel.SQC.Absolute.pdf Output/SigProfiler/Bowtie2-Distribution/COSMIC-Indel.SQC.Relative.pdf
step33: step33-SigProfiler step33-Violin-Subtype step33-Violin-Recurrence step33-Violin-Smoking step33-Bar-SQC step33-Bar-ADC step33-Bar-Subtype step33-Bar-Smoking step33-Bar-Recurrence step33-Distribution
.PHONY: step33-SigProfiler step33-Violin-Subtype step33-Violin-Recurrence step33-Violin-Smoking step33-Bar-SQC step33-Bar-ADC step33-Bar-Subtype step33-Bar-Smoking step33-Bar-Recurrence step33-Distribution step33

# Step 34 (Short variants analysis upon clinical data with WES data)
Output/Gene_Clinical:
	mkdir -p $@

Output/Gene_Clinical/BWA-SQC.Recur.tsv Output/Gene_Clinical/BWA-SQC.Recur.pdf: Python/plot_gene_clinical.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/BWA-ADC.Recur.tsv Output/Gene_Clinical/BWA-ADC.Recur.pdf: Python/plot_gene_clinical.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/Bowtie2-SQC.Recur.tsv Output/Gene_Clinical/Bowtie2-SQC.Recur.pdf: Python/plot_gene_clinical.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/Bowtie2-ADC.Recur.tsv Output/Gene_Clinical/Bowtie2-ADC.Recur.pdf: Python/plot_gene_clinical.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/BWA-SQC.CIS.Recur.tsv Output/Gene_Clinical/BWA-SQC.CIS.Recur.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "CIS+AIS" --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/BWA-SQC.Precancer.Recur.tsv Output/Gene_Clinical/BWA-SQC.Precancer.Recur.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "Precancer" --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/BWA-SQC.Primary.Recur.tsv Output/Gene_Clinical/BWA-SQC.Primary.Recur.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "Primary" --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/BWA-ADC.AAH.Recur.tsv Output/Gene_Clinical/BWA-ADC.AAH.Recur.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "AAH" --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/BWA-ADC.Precancer.Recur.tsv Output/Gene_Clinical/BWA-ADC.Precancer.Recur.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "Precancer" --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/BWA-ADC.Primary.Recur.tsv Output/Gene_Clinical/BWA-ADC.Primary.Recur.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "Primary" --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/Bowtie2-SQC.CIS.Recur.tsv Output/Gene_Clinical/Bowtie2-SQC.CIS.Recur.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "CIS+AIS" --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/Bowtie2-SQC.Precancer.Recur.tsv Output/Gene_Clinical/Bowtie2-SQC.Precancer.Recur.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "Precancer" --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/Bowtie2-SQC.Primary.Recur.tsv Output/Gene_Clinical/Bowtie2-SQC.Primary.Recur.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "Primary" --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/Bowtie2-ADC.AAH.Recur.tsv Output/Gene_Clinical/Bowtie2-ADC.AAH.Recur.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "AAH" --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/Bowtie2-ADC.Precancer.Recur.tsv Output/Gene_Clinical/Bowtie2-ADC.Precancer.Recur.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "Precancer" --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/Bowtie2-ADC.Primary.Recur.tsv Output/Gene_Clinical/Bowtie2-ADC.Primary.Recur.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "Primary" --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/BWA-SQC.Smoking.tsv Output/Gene_Clinical/BWA-SQC.Smoking.pdf: Python/plot_gene_clinical.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/BWA-ADC.Smoking.tsv Output/Gene_Clinical/BWA-ADC.Smoking.pdf: Python/plot_gene_clinical.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --ADC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/Bowtie2-SQC.Smoking.tsv Output/Gene_Clinical/Bowtie2-SQC.Smoking.pdf: Python/plot_gene_clinical.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/Bowtie2-ADC.Smoking.tsv Output/Gene_Clinical/Bowtie2-ADC.Smoking.pdf: Python/plot_gene_clinical.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --ADC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/BWA-SQC.CIS.Smoking.tsv Output/Gene_Clinical/BWA-SQC.CIS.Smoking.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "CIS+AIS" --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/BWA-SQC.Precancer.Smoking.tsv Output/Gene_Clinical/BWA-SQC.Precancer.Smoking.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "Precancer" --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/BWA-SQC.Primary.Smoking.tsv Output/Gene_Clinical/BWA-SQC.Primary.Smoking.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "Primary" --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/BWA-ADC.Precancer.Smoking.tsv Output/Gene_Clinical/BWA-ADC.Precancer.Smoking.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "Precancer" --cpus $(CPUS) --ADC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/BWA-ADC.Primary.Smoking.tsv Output/Gene_Clinical/BWA-ADC.Primary.Smoking.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "Primary" --cpus $(CPUS) --ADC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/Bowtie2-SQC.CIS.Smoking.tsv Output/Gene_Clinical/Bowtie2-SQC.CIS.Smoking.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "CIS+AIS" --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/Bowtie2-SQC.Precancer.Smoking.tsv Output/Gene_Clinical/Bowtie2-SQC.Precancer.Smoking.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "Precancer" --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/Bowtie2-SQC.Primary.Smoking.tsv Output/Gene_Clinical/Bowtie2-SQC.Primary.Smoking.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "Primary" --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical/Bowtie2-ADC.Precancer.Smoking.tsv Output/Gene_Clinical/Bowtie2-ADC.Precancer.Smoking.pdf: Python/plot_gene_clinical_2.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --stage "Precancer" --cpus $(CPUS) --ADC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

step34-BWA-Recur: Output/Gene_Clinical/BWA-SQC.Recur.pdf Output/Gene_Clinical/BWA-ADC.Recur.pdf
step34-Bowtie2-Recur: Output/Gene_Clinical/Bowtie2-SQC.Recur.pdf Output/Gene_Clinical/Bowtie2-ADC.Recur.pdf
step34-BWA-Stage-Recur: Output/Gene_Clinical/BWA-SQC.CIS.Recur.pdf Output/Gene_Clinical/BWA-SQC.Precancer.Recur.pdf Output/Gene_Clinical/BWA-SQC.Primary.Recur.pdf Output/Gene_Clinical/BWA-ADC.AAH.Recur.pdf Output/Gene_Clinical/BWA-ADC.Precancer.Recur.pdf Output/Gene_Clinical/BWA-ADC.Primary.Recur.pdf
step34-Bowtie2-Stage-Recur: Output/Gene_Clinical/Bowtie2-SQC.CIS.Recur.pdf Output/Gene_Clinical/Bowtie2-SQC.Precancer.Recur.pdf Output/Gene_Clinical/Bowtie2-SQC.Primary.Recur.pdf Output/Gene_Clinical/Bowtie2-ADC.AAH.Recur.pdf Output/Gene_Clinical/Bowtie2-ADC.Precancer.Recur.pdf Output/Gene_Clinical/Bowtie2-ADC.Primary.Recur.pdf
step34-BWA-Smoking: Output/Gene_Clinical/BWA-SQC.Smoking.pdf Output/Gene_Clinical/BWA-ADC.Smoking.pdf
step34-Bowtie2-Smoking: Output/Gene_Clinical/Bowtie2-SQC.Smoking.pdf Output/Gene_Clinical/Bowtie2-ADC.Smoking.pdf
step34-BWA-Stage-Smoking: Output/Gene_Clinical/BWA-SQC.CIS.Smoking.pdf Output/Gene_Clinical/BWA-SQC.Precancer.Smoking.pdf Output/Gene_Clinical/BWA-SQC.Primary.Smoking.pdf Output/Gene_Clinical/BWA-ADC.Precancer.Smoking.pdf Output/Gene_Clinical/BWA-ADC.Primary.Smoking.pdf
step34-Bowtie2-Stage-Smoking: Output/Gene_Clinical/Bowtie2-SQC.CIS.Smoking.pdf Output/Gene_Clinical/Bowtie2-SQC.Precancer.Smoking.pdf Output/Gene_Clinical/Bowtie2-SQC.Primary.Smoking.pdf Output/Gene_Clinical/Bowtie2-ADC.Precancer.Smoking.pdf Output/Gene_Clinical/Bowtie2-ADC.Primary.Smoking.pdf
step34: step34-BWA-Recur step34-Bowtie2-Recur step34-BWA-Stage-Recur step34-Bowtie2-Stage-Recur step34-BWA-Smoking step34-Bowtie2-Smoking step34-BWA-Stage-Smoking step34-Bowtie2-Stage-Smoking
.PHONY: step34-BWA-Recur step34-Bowtie2-Recur step34-BWA-Stage-Recur step34-Bowtie2-Stage-Recur step34-BWA-Smoking step34-Bowtie2-Smoking step34-BWA-Stage-Smoking step34-Bowtie2-Stage-Smoking step34

# Step 35 (Gene expression changing upon clinical data)
Output/GeneExpression:
	mkdir -p $@

Output/GeneExpression/%.SQC.Subtype.tar: Python/plot_expression_clinical.py Output/RSEM/DEG/%.TPM.tsv $(CLINICAL_DATA_CSV) | build.log Output/GeneExpression
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/GeneExpression/%.ADC.Subtype.tar: Python/plot_expression_clinical.py Output/RSEM/DEG/%.TPM.tsv $(CLINICAL_DATA_CSV) | build.log Output/GeneExpression
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/GeneExpression/%.SQC.Recurrence.tar: Python/plot_expression_clinical_2.py Output/RSEM/DEG/%.TPM.tsv $(CLINICAL_DATA_CSV) | build.log Output/GeneExpression
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/GeneExpression/%.ADC.Recurrence.tar: Python/plot_expression_clinical_2.py Output/RSEM/DEG/%.TPM.tsv $(CLINICAL_DATA_CSV) | build.log Output/GeneExpression
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

step35-Subtype: Output/GeneExpression/STAR.SQC.Subtype.tar Output/GeneExpression/STAR.ADC.Subtype.tar Output/GeneExpression/Bowtie2.SQC.Subtype.tar Output/GeneExpression/Bowtie2.ADC.Subtype.tar
step35-Recur: Output/GeneExpression/STAR.SQC.Recurrence.tar Output/GeneExpression/STAR.ADC.Recurrence.tar Output/GeneExpression/Bowtie2.SQC.Recurrence.tar Output/GeneExpression/Bowtie2.ADC.Recurrence.tar
step35: step35-Subtype step35-Recur
.PHONY: step35-Subtype step35-Recur step35

# Step 36 (Picard - CollectMultipleMetrics)
Output/CollectMultipleMetrics/BWA Output/CollectMultipleMetrics/Bowtie2:
	mkdir -p $@

Output/CollectMultipleMetrics/BWA/%.alignment_summary_metrics.txt Output/CollectMultipleMetrics/BWA/%.base_distribution_by_cycle_metrics.txt Output/CollectMultipleMetrics/BWA/%.base_distribution_by_cycle.pdf Output/CollectMultipleMetrics/BWA/%.insert_size_histogram.pdf Output/CollectMultipleMetrics/BWA/%.insert_size_metrics.txt Output/CollectMultipleMetrics/BWA/%.quality_by_cycle_metrics.txt Output/CollectMultipleMetrics/BWA/%.quality_by_cycle.pdf Output/CollectMultipleMetrics/BWA/%.quality_distribution_metrics.txt Output/CollectMultipleMetrics/BWA/%.quality_distribution.pdf Output/CollectMultipleMetrics/BWA/%.read_length_histogram.pdf: $(REFERENCE_FASTA) $(REFERENCE_DBSNP_VCF) Output/WES/%.BWA.Sort.MarkDuplicate.BQSR.bam | build.log Output/CollectMultipleMetrics/BWA
	$(DOCKER) $(PICARD) CollectMultipleMetrics --INPUT $(addprefix /,$(filter %.bam,$^)) --OUTPUT $(addprefix /,$(@D)/$*) --DB_SNP $(addprefix /,$(filter %.vcf,$^)) --FILE_EXTENSION ".txt" --COMPRESSION_LEVEL 9 --REFERENCE_SEQUENCE $(addprefix /,$(filter %.fasta,$^)) --VALIDATION_STRINGENCY SILENT 1> $(@D)/$*.stdout 2> $(@D)/$*.stderr

Output/CollectMultipleMetrics/Bowtie2/%.alignment_summary_metrics.txt Output/CollectMultipleMetrics/Bowtie2/%.base_distribution_by_cycle_metrics.txt Output/CollectMultipleMetrics/Bowtie2/%.base_distribution_by_cycle.pdf Output/CollectMultipleMetrics/Bowtie2/%.insert_size_histogram.pdf Output/CollectMultipleMetrics/Bowtie2/%.insert_size_metrics.txt Output/CollectMultipleMetrics/Bowtie2/%.quality_by_cycle_metrics.txt Output/CollectMultipleMetrics/Bowtie2/%.quality_by_cycle.pdf Output/CollectMultipleMetrics/Bowtie2/%.quality_distribution_metrics.txt Output/CollectMultipleMetrics/Bowtie2/%.quality_distribution.pdf Output/CollectMultipleMetrics/Bowtie2/%.read_length_histogram.pdf: $(REFERENCE_FASTA) $(REFERENCE_DBSNP_VCF) Output/WES/%.Bowtie2.Sort.MarkDuplicate.BQSR.bam | build.log Output/CollectMultipleMetrics/Bowtie2
	$(DOCKER) $(PICARD) CollectMultipleMetrics --INPUT $(addprefix /,$(filter %.bam,$^)) --OUTPUT $(addprefix /,$(@D)/$*) --DB_SNP $(addprefix /,$(filter %.vcf,$^)) --FILE_EXTENSION ".txt" --COMPRESSION_LEVEL 9 --REFERENCE_SEQUENCE $(addprefix /,$(filter %.fasta,$^)) --VALIDATION_STRINGENCY SILENT 1> $(@D)/$*.stdout 2> $(@D)/$*.stderr

Output/CollectMultipleMetrics/BWA.pdf: Python/analysis_quality_distribution.py $(addprefix Output/CollectMultipleMetrics/BWA/,$(addsuffix .quality_distribution_metrics.txt,$(WES_ID))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/CollectMultipleMetrics/Bowtie2.pdf: Python/analysis_quality_distribution.py $(addprefix Output/CollectMultipleMetrics/Bowtie2/,$(addsuffix .quality_distribution_metrics.txt,$(WES_ID))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step36-BWA: $(addprefix Output/CollectMultipleMetrics/BWA/,$(addsuffix .quality_distribution_metrics.txt,$(WES_ID))) Output/CollectMultipleMetrics/BWA.pdf
step36-Bowtie2: $(addprefix Output/CollectMultipleMetrics/Bowtie2/,$(addsuffix .quality_distribution_metrics.txt,$(WES_ID))) Output/CollectMultipleMetrics/Bowtie2.pdf
step36: step36-BWA step36-Bowtie2
.PHONY: step36-BWA step36-Bowtie2 step36

# Step 37 (PureCN)
# https://www.bioconductor.org/packages/devel/bioc/vignettes/PureCN/inst/doc/Quick.html
Output/PureCN/Reference Output/PureCN/BWA Output/PureCN/Bowtie2:
	mkdir -p $@

Output/PureCN/Reference/hg38.txt Output/PureCN/Reference/hg38.bed: $(REFERENCE_FASTA) $(LIBRARY_BED) $(MAPPABILITY_BW) | Output/PureCN/Reference build.log
	$(DOCKER) $(R) $(PURECN_EXTDATA_DIR)/IntervalFile.R --fasta=$(addprefix /,$(filter %.fasta,$^)) --infile=$(addprefix /,$(filter %.bed,$^)) --offtarget --mappability=$(addprefix /,$(filter %.bw,$^)) --genome="hg38" --outfile=$(addprefix /,$(@D)/hg38.txt) --export=$(addprefix /,$(@D)/hg38.bed) --force 1> $@.stdout 2> $@.stderr

Output/PureCN/Reference/Mappability.bed: $(MAPPABILITY_BW) | build.log
	$(DOCKER) $(R) $(PURECN_EXTDATA_DIR)/FilterCallableLoci.R --in-file $(addprefix /,$(filter %.bw,$^)) --genome="hg38" --out-file=$(addprefix /,$@) --force 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA/%.BWA.Sort.MarkDuplicate.BQSR_coverage.txt.gz: Output/WES/%.BWA.Sort.MarkDuplicate.BQSR.bam Output/PureCN/Reference/hg38.txt | Output/PureCN/BWA build.log
	$(DOCKER) $(R) $(PURECN_EXTDATA_DIR)/Coverage.R --bam=$(addprefix /,$(filter %.bam,$^)) --intervals=$(addprefix /,$(filter %.txt,$^)) --out-dir=$(addprefix /,$(@D)) --force --skip-gc-norm --parallel --cores=$(CPUS) --seed="12345" 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2/%.Bowtie2.Sort.MarkDuplicate.BQSR_coverage.txt.gz: Output/WES/%.Bowtie2.Sort.MarkDuplicate.BQSR.bam Output/PureCN/Reference/hg38.txt | Output/PureCN/Bowtie2 build.log
	$(DOCKER) $(R) $(PURECN_EXTDATA_DIR)/Coverage.R --bam=$(addprefix /,$(filter %.bam,$^)) --intervals=$(addprefix /,$(filter %.txt,$^)) --out-dir=$(addprefix /,$(@D)) --force --skip-gc-norm 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA/Normal_coverages.list: $(addprefix Output/PureCN/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR_coverage.txt.gz,$(WES_NORMAL_ID))) | build.log
	$(DOCKER) ls -1 $(addprefix /,$(sort $^)) > $@

Output/PureCN/Bowtie2/Normal_coverages.list: $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR_coverage.txt.gz,$(WES_NORMAL_ID))) | build.log
	$(DOCKER) ls -1 $(addprefix /,$(sort $^)) > $@

Output/PureCN/BWA/interval_weights_hg38.png Output/PureCN/BWA/normalDB_hg38.rds Output/PureCN/BWA/low_coverage_targets_hg38.bed: Output/PureCN/BWA/Normal_coverages.list | build.log
	$(DOCKER) $(R) $(PURECN_EXTDATA_DIR)/NormalDB.R --coverage-files=$(addprefix /,$(filter %.list,$^)) --genome="hg38" --out-dir=$(addprefix /,$(@D)) --force 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2/interval_weights_hg38.png Output/PureCN/Bowtie2/normalDB_hg38.rds Output/PureCN/Bowtie2/low_coverage_targets_hg38.bed: Output/PureCN/Bowtie2/Normal_coverages.list | build.log
	$(DOCKER) $(R) $(PURECN_EXTDATA_DIR)/NormalDB.R --coverage-files=$(addprefix /,$(filter %.list,$^)) --genome="hg38" --out-dir=$(addprefix /,$(@D)) --force 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA/%_segmentation.pdf Output/PureCN/BWA/%.csv Output/PureCN/BWA/%.pdf Output/PureCN/BWA/%_local_optima.pdf Output/PureCN/BWA/%_dnacopy.seg Output/PureCN/BWA/%.log Output/PureCN/BWA/%.rds: Output/PureCN/Reference/hg38.txt Output/PureCN/BWA/%.BWA.Sort.MarkDuplicate.BQSR_coverage.txt.gz Output/PureCN/BWA/normalDB_hg38.rds | build.log
	$(DOCKER) $(R) $(PURECN_EXTDATA_DIR)/PureCN.R --sampleid=$* --tumor=$(addprefix /,$(filter %.txt.gz,$^)) --normaldb=$(addprefix /,$(filter %.rds,$^)) --genome="hg38" --parallel --cores=$(CPUS) --out=$(addprefix /,$(@D)) --force --post-optimize --seed="12345" --intervals=$(addprefix /,$(filter %.txt,$^)) --model="betabin" --out-vcf 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2/%_segmentation.pdf Output/PureCN/Bowtie2/%.csv Output/PureCN/Bowtie2/%.pdf Output/PureCN/Bowtie2/%_local_optima.pdf Output/PureCN/Bowtie2/%_dnacopy.seg Output/PureCN/Bowtie2/%.log Output/PureCN/Bowtie2/%.rds: Output/PureCN/Reference/hg38.txt Output/PureCN/Bowtie2/%.Bowtie2.Sort.MarkDuplicate.BQSR_coverage.txt.gz Output/PureCN/Bowtie2/normalDB_hg38.rds | build.log
	$(DOCKER) $(R) $(PURECN_EXTDATA_DIR)/PureCN.R --sampleid=$* --tumor=$(addprefix /,$(filter %.txt.gz,$^)) --normaldb=$(addprefix /,$(filter %.rds,$^)) --genome="hg38" --parallel --cores=$(CPUS) --out=$(addprefix /,$(@D)) --force --post-optimize --seed="12345" --intervals=$(addprefix /,$(filter %.txt,$^)) --model="betabin" --out-vcf 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA/%.seg.tsv: R/get_pureCN.R Output/PureCN/BWA/%.rds | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.rds,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2/%.seg.tsv: R/get_pureCN.R Output/PureCN/Bowtie2/%.rds | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.rds,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-PureCN-SQC.pdf: Python/aggregate_purecn.py $(addprefix Output/PureCN/BWA/,$(addsuffix .csv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-PureCN-ADC.pdf: Python/aggregate_purecn.py $(addprefix Output/PureCN/BWA/,$(addsuffix .csv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-PureCN-SQC.pdf: Python/aggregate_purecn.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .csv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-PureCN-ADC.pdf: Python/aggregate_purecn.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .csv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-genome-SQC.pdf Output/PureCN/BWA-genome-SQC.png: Python/aggregate_PureCN_genome.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-genome-ADC.pdf Output/PureCN/BWA-genome-ADC.png: Python/aggregate_PureCN_genome.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-genome-SQC.pdf Output/PureCN/Bowtie2-genome-SQC.png: Python/aggregate_PureCN_genome.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-genome-ADC.pdf Output/PureCN/Bowtie2-genome-ADC.png: Python/aggregate_PureCN_genome.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-genome-SQC.patient.pdf Output/PureCN/BWA-genome-SQC.patient.png: Python/aggregate_PureCN_genome.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-genome-ADC.patient.pdf Output/PureCN/BWA-genome-ADC.patient.png: Python/aggregate_PureCN_genome.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-genome-SQC.patient.pdf Output/PureCN/Bowtie2-genome-SQC.patient.png: Python/aggregate_PureCN_genome.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-genome-ADC.patient.pdf Output/PureCN/Bowtie2-genome-ADC.patient.png: Python/aggregate_PureCN_genome.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-simple-SQC.pdf: Python/aggregate_PureCN_simple.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-simple-ADC.pdf: Python/aggregate_PureCN_simple.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-simple-SQC.pdf: Python/aggregate_PureCN_simple.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-simple-ADC.pdf: Python/aggregate_PureCN_simple.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-genome-SQC.Recurrence.pdf Output/PureCN/BWA-genome-SQC.Recurrence.png: Python/aggregate_PureCN_genome_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-genome-ADC.Recurrence.pdf Output/PureCN/BWA-genome-ADC.Recurrence.png: Python/aggregate_PureCN_genome_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-genome-SQC.Recurrence.pdf Output/PureCN/Bowtie2-genome-SQC.Recurrence.png: Python/aggregate_PureCN_genome_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-genome-ADC.Recurrence.pdf Output/PureCN/Bowtie2-genome-ADC.Recurrence.png: Python/aggregate_PureCN_genome_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-genome-SQC.Recurrence.patient.pdf Output/PureCN/BWA-genome-SQC.Recurrence.patient.png: Python/aggregate_PureCN_genome_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-genome-ADC.Recurrence.patient.pdf Output/PureCN/BWA-genome-ADC.Recurrence.patient.png: Python/aggregate_PureCN_genome_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-genome-SQC.Recurrence.patient.pdf Output/PureCN/Bowtie2-genome-SQC.Recurrence.patient.png: Python/aggregate_PureCN_genome_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-genome-ADC.Recurrence.patient.pdf Output/PureCN/Bowtie2-genome-ADC.Recurrence.patient.png: Python/aggregate_PureCN_genome_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-genome-SQC.Smoking.pdf Output/PureCN/BWA-genome-SQC.Smoking.png: Python/aggregate_PureCN_genome_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-genome-ADC.Smoking.pdf Output/PureCN/BWA-genome-ADC.Smoking.png: Python/aggregate_PureCN_genome_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-genome-SQC.Smoking.pdf Output/PureCN/Bowtie2-genome-SQC.Smoking.png: Python/aggregate_PureCN_genome_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-genome-ADC.Smoking.pdf Output/PureCN/Bowtie2-genome-ADC.Smoking.png: Python/aggregate_PureCN_genome_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-genome-SQC.Smoking.patient.pdf Output/PureCN/BWA-genome-SQC.Smoking.patient.png: Python/aggregate_PureCN_genome_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-genome-ADC.Smoking.patient.pdf Output/PureCN/BWA-genome-ADC.Smoking.patient.png: Python/aggregate_PureCN_genome_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-genome-SQC.Smoking.patient.pdf Output/PureCN/Bowtie2-genome-SQC.Smoking.patient.png: Python/aggregate_PureCN_genome_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-genome-ADC.Smoking.patient.pdf Output/PureCN/Bowtie2-genome-ADC.Smoking.patient.png: Python/aggregate_PureCN_genome_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-simple-SQC.Recurrence.pdf: Python/aggregate_PureCN_simple_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-simple-ADC.Recurrence.pdf: Python/aggregate_PureCN_simple_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-simple-SQC.Recurrence.pdf: Python/aggregate_PureCN_simple_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-simple-ADC.Recurrence.pdf: Python/aggregate_PureCN_simple_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-simple-SQC.Smoking.pdf: Python/aggregate_PureCN_simple_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-simple-ADC.Smoking.pdf: Python/aggregate_PureCN_simple_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-simple-SQC.Smoking.pdf: Python/aggregate_PureCN_simple_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-simple-ADC.Smoking.pdf: Python/aggregate_PureCN_simple_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-violin-SQC.pdf: Python/aggregate_PureCN_violin.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-violin-ADC.pdf: Python/aggregate_PureCN_violin.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-violin-SQC.pdf: Python/aggregate_PureCN_violin.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-violin-ADC.pdf: Python/aggregate_PureCN_violin.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-violin-SQC.Recurrence.pdf: Python/aggregate_PureCN_violin_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-violin-ADC.Recurrence.pdf: Python/aggregate_PureCN_violin_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-violin-SQC.Recurrence.pdf: Python/aggregate_PureCN_violin_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-violin-ADC.Recurrence.pdf: Python/aggregate_PureCN_violin_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-violin-SQC.Smoking.pdf: Python/aggregate_PureCN_violin_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-violin-ADC.Smoking.pdf: Python/aggregate_PureCN_violin_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-violin-SQC.Smoking.pdf: Python/aggregate_PureCN_violin_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-violin-ADC.Smoking.pdf: Python/aggregate_PureCN_violin_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-chromosome-SQC.tar: Python/aggregate_PureCN_chromosome.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(GENE_CENCUS_ONE_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-chromosome-ADC.tar: Python/aggregate_PureCN_chromosome.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(GENE_CENCUS_ONE_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-chromosome-SQC.tar: Python/aggregate_PureCN_chromosome.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(GENE_CENCUS_ONE_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-chromosome-ADC.tar: Python/aggregate_PureCN_chromosome.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(GENE_CENCUS_ONE_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-chromosome-SQC.Recurrence.tar: Python/aggregate_PureCN_chromosome_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(GENE_CENCUS_ONE_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-chromosome-ADC.Recurrence.tar: Python/aggregate_PureCN_chromosome_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(GENE_CENCUS_ONE_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-chromosome-SQC.Smoking.tar: Python/aggregate_PureCN_chromosome_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(GENE_CENCUS_ONE_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-chromosome-ADC.Smoking.tar: Python/aggregate_PureCN_chromosome_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(GENE_CENCUS_ONE_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-chromosome-SQC.Recurrence.tar: Python/aggregate_PureCN_chromosome_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(GENE_CENCUS_ONE_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-chromosome-ADC.Recurrence.tar: Python/aggregate_PureCN_chromosome_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(GENE_CENCUS_ONE_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-chromosome-SQC.Smoking.tar: Python/aggregate_PureCN_chromosome_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(GENE_CENCUS_ONE_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-chromosome-ADC.Smoking.tar: Python/aggregate_PureCN_chromosome_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(GENE_CENCUS_ONE_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-gene-SQC.tar: Python/aggregate_PureCN_violin_gene.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-gene-ADC.tar: Python/aggregate_PureCN_violin_gene.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-gene-SQC.tar: Python/aggregate_PureCN_violin_gene.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-gene-ADC.tar: Python/aggregate_PureCN_violin_gene.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-gene-SQC.Recurrence.tar: Python/aggregate_PureCN_violin_gene_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-gene-ADC.Recurrence.tar: Python/aggregate_PureCN_violin_gene_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-gene-SQC.Recurrence.tar: Python/aggregate_PureCN_violin_gene_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-gene-ADC.Recurrence.tar: Python/aggregate_PureCN_violin_gene_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-gene-SQC.Smoking.tar: Python/aggregate_PureCN_violin_gene_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/PureCN/BWA-gene-ADC.Smoking.tar: Python/aggregate_PureCN_violin_gene_clinical.py $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-gene-SQC.Smoking.tar: Python/aggregate_PureCN_violin_gene_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

Output/PureCN/Bowtie2-gene-ADC.Smoking.tar: Python/aggregate_PureCN_violin_gene_clinical.py $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking" "Never" "Ever" 1> $@.stdout 2> $@.stderr

step37-Reference: Output/PureCN/Reference/hg38.txt# Output/PureCN/Referece/Mappability.bed
step37-BWA: $(addprefix Output/PureCN/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR_coverage.txt.gz,$(WES_ID))) Output/PureCN/BWA/Normal_coverages.list Output/PureCN/BWA/normalDB_hg38.rds $(addprefix Output/PureCN/BWA/,$(addsuffix .csv,$(WES_REAL_TUMOR_ID))) $(addprefix Output/PureCN/BWA/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID)))
step37-Bowtie2: $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR_coverage.txt.gz,$(WES_ID))) Output/PureCN/Bowtie2/Normal_coverages.list Output/PureCN/Bowtie2/normalDB_hg38.rds $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .csv,$(WES_REAL_TUMOR_ID))) $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .seg.tsv,$(WES_REAL_TUMOR_ID)))
step37-default: Output/PureCN/BWA-PureCN-SQC.pdf Output/PureCN/BWA-PureCN-ADC.pdf Output/PureCN/Bowtie2-PureCN-SQC.pdf Output/PureCN/Bowtie2-PureCN-ADC.pdf
step37-genome: Output/PureCN/BWA-genome-SQC.pdf Output/PureCN/BWA-genome-ADC.pdf Output/PureCN/Bowtie2-genome-SQC.pdf Output/PureCN/Bowtie2-genome-ADC.pdf Output/PureCN/BWA-genome-SQC.patient.pdf Output/PureCN/BWA-genome-ADC.patient.pdf Output/PureCN/Bowtie2-genome-SQC.patient.pdf Output/PureCN/Bowtie2-genome-ADC.patient.pdf
step37-simple: Output/PureCN/BWA-simple-SQC.pdf Output/PureCN/BWA-simple-ADC.pdf Output/PureCN/Bowtie2-simple-SQC.pdf Output/PureCN/Bowtie2-simple-ADC.pdf
step37-genome-clinical: Output/PureCN/BWA-genome-SQC.Recurrence.pdf Output/PureCN/BWA-genome-ADC.Recurrence.pdf Output/PureCN/Bowtie2-genome-SQC.Recurrence.pdf Output/PureCN/Bowtie2-genome-ADC.Recurrence.pdf Output/PureCN/BWA-genome-SQC.Smoking.pdf Output/PureCN/BWA-genome-ADC.Smoking.pdf Output/PureCN/Bowtie2-genome-SQC.Smoking.pdf Output/PureCN/Bowtie2-genome-ADC.Smoking.pdf Output/PureCN/BWA-genome-SQC.Recurrence.patient.pdf Output/PureCN/BWA-genome-ADC.Recurrence.patient.pdf Output/PureCN/Bowtie2-genome-SQC.Recurrence.patient.pdf Output/PureCN/Bowtie2-genome-ADC.Recurrence.patient.pdf Output/PureCN/BWA-genome-SQC.Smoking.patient.pdf Output/PureCN/BWA-genome-ADC.Smoking.patient.pdf Output/PureCN/Bowtie2-genome-SQC.Smoking.patient.pdf Output/PureCN/Bowtie2-genome-ADC.Smoking.patient.pdf
step37-simple-clinical: Output/PureCN/BWA-simple-SQC.Recurrence.pdf Output/PureCN/BWA-simple-ADC.Recurrence.pdf Output/PureCN/Bowtie2-simple-SQC.Recurrence.pdf Output/PureCN/Bowtie2-simple-ADC.Recurrence.pdf Output/PureCN/BWA-simple-SQC.Smoking.pdf Output/PureCN/BWA-simple-ADC.Smoking.pdf Output/PureCN/Bowtie2-simple-SQC.Smoking.pdf Output/PureCN/Bowtie2-simple-ADC.Smoking.pdf
step37-violin: Output/PureCN/BWA-violin-SQC.pdf Output/PureCN/BWA-violin-ADC.pdf Output/PureCN/Bowtie2-violin-SQC.pdf Output/PureCN/Bowtie2-violin-ADC.pdf
step37-violin-clinical: Output/PureCN/BWA-violin-SQC.Recurrence.pdf Output/PureCN/BWA-violin-ADC.Recurrence.pdf Output/PureCN/Bowtie2-violin-SQC.Recurrence.pdf Output/PureCN/Bowtie2-violin-ADC.pdf Output/PureCN/BWA-violin-SQC.Smoking.pdf Output/PureCN/BWA-violin-ADC.Smoking.pdf Output/PureCN/Bowtie2-violin-SQC.Smoking.pdf Output/PureCN/Bowtie2-violin-ADC.Smoking.pdf
step37-chromosome: Output/PureCN/BWA-chromosome-SQC.tar Output/PureCN/BWA-chromosome-ADC.tar Output/PureCN/Bowtie2-chromosome-SQC.tar Output/PureCN/Bowtie2-chromosome-ADC.tar
step37-chromosome-clinical: Output/PureCN/BWA-chromosome-SQC.Recurrence.tar Output/PureCN/BWA-chromosome-ADC.Recurrence.tar Output/PureCN/BWA-chromosome-SQC.Smoking.tar Output/PureCN/BWA-chromosome-ADC.Smoking.tar Output/PureCN/Bowtie2-chromosome-SQC.Recurrence.tar Output/PureCN/Bowtie2-chromosome-ADC.Recurrence.tar Output/PureCN/Bowtie2-chromosome-SQC.Smoking.tar Output/PureCN/Bowtie2-chromosome-ADC.Smoking.tar
step37-gene: Output/PureCN/BWA-gene-SQC.tar Output/PureCN/BWA-gene-ADC.tar Output/PureCN/Bowtie2-gene-SQC.tar Output/PureCN/Bowtie2-gene-ADC.tar
step37-gene-clinical: Output/PureCN/BWA-gene-SQC.Recurrence.tar Output/PureCN/BWA-gene-ADC.Recurrence.tar Output/PureCN/Bowtie2-gene-SQC.Recurrence.tar Output/PureCN/Bowtie2-gene-ADC.Recurrence.tar Output/PureCN/BWA-gene-SQC.Smoking.tar Output/PureCN/BWA-gene-ADC.Smoking.tar Output/PureCN/Bowtie2-gene-SQC.Smoking.tar Output/PureCN/Bowtie2-gene-ADC.Smoking.tar
step37: step37-Reference step37-BWA step37-Bowtie2 step37-default step37-genome step37-simple step37-genome-clinical step37-simple-clinical step37-violin step37-violin-clinical step37-chromosome step37-chromosome-clinical step37-gene step37-gene-clinical
.PHONY: step37-Reference step37-BWA step37-Bowtie2 step37-default step37-genome step37-simple step37-genome-clinical step37-simple-clinical step37-violin step37-violin-clinical step37-chromosome step37-chromosome-clinical step37-gene step37-gene-clinical step37

# Step 38 (Short variants analysis upon clinical data with WES data as Venn diagram)
Output/Gene_Clinical_Venn:
	mkdir -p $@

Output/Gene_Clinical_Venn/BWA-SQC.pdf Output/Gene_Clinical_Venn/BWA-SQC.tsv: Python/plot_gene_clinical_venn.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical_Venn
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).pdf $(basename $@).tsv) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Venn/BWA-ADC.pdf Output/Gene_Clinical_Venn/BWA-ADC.tsv: Python/plot_gene_clinical_venn.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical_Venn
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).pdf $(basename $@).tsv) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Venn/Bowtie2-SQC.pdf Output/Gene_Clinical_Venn/Bowtie2-SQC.tsv: Python/plot_gene_clinical_venn.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical_Venn
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).pdf $(basename $@).tsv) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Venn/Bowtie2-ADC.pdf Output/Gene_Clinical_Venn/Bowtie2-ADC.tsv: Python/plot_gene_clinical_venn.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Gene_Clinical_Venn
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).pdf $(basename $@).tsv) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

step38: Output/Gene_Clinical_Venn/BWA-SQC.pdf Output/Gene_Clinical_Venn/BWA-ADC.pdf Output/Gene_Clinical_Venn/Bowtie2-SQC.pdf Output/Gene_Clinical_Venn/Bowtie2-ADC.pdf
.PHONY: step38

# Step 39 (CNVkit)
Output/CNVkit/BWA Output/CNVkit/Bowtie2:
	mkdir -p $@

Output/CNVkit/BWA/Reference.tsv: $(REFERENCE_FASTA) $(LIBRARY_BED) $(UCSC_REFFLAT_TXT) $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.bam,$(WES_ID))) | build.log Output/CNVkit/BWA
	$(DOCKER) cnvkit.py batch --processes $(CPUS) --normal $(addprefix /,$(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.bam,$(WES_NORMAL_ID)))) --fasta $(addprefix /,$(filter %.fasta,$^)) --targets $(addprefix /,$(filter %.bed,$^)) --annotate $(addprefix /,$(filter %.txt,$^)) --output-reference $(addprefix /,$@) --scatter --diagram --cluster --output-dir $(addprefix /,$(@D)) $(addprefix /,$(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.bam,$(WES_TUMOR_ID)))) 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2/Reference.tsv: $(REFERENCE_FASTA) $(LIBRARY_BED) $(UCSC_REFFLAT_TXT) $(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.bam,$(WES_ID))) | build.log Output/CNVkit/Bowtie2
	$(DOCKER) cnvkit.py batch --processes $(CPUS) --normal $(addprefix /,$(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.bam,$(WES_NORMAL_ID)))) --fasta $(addprefix /,$(filter %.fasta,$^)) --targets $(addprefix /,$(filter %.bed,$^)) --annotate $(addprefix /,$(filter %.txt,$^)) --output-reference $(addprefix /,$@) --scatter --diagram --cluster --output-dir $(addprefix /,$(@D)) $(addprefix /,$(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.bam,$(WES_TUMOR_ID)))) 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA/%.BWA.Sort.MarkDuplicate.BQSR.cnr Output/CNVkit/BWA/%.BWA.Sort.MarkDuplicate.BQSR.cns: Output/CNVkit/BWA/Reference.tsv
	$(DOCKER) touch $(addprefix /,$@)

Output/CNVkit/Bowtie2/%.Bowtie2.Sort.MarkDuplicate.BQSR.cnr Output/CNVkit/Bowtie2/%.Bowtie2.Sort.MarkDuplicate.BQSR.cns: Output/CNVkit/Bowtie2/Reference.tsv
	$(DOCKER) touch $(addprefix /,$@)

Output/CNVkit/BWA-genome-SQC.pdf Output/CNVkit/BWA-genome-SQC.png: Python/aggregate_CNVkit_genome.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-genome-ADC.pdf Output/CNVkit/BWA-genome-ADC.png: Python/aggregate_CNVkit_genome.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-genome-SQC.pdf Output/CNVkit/Bowtie2-genome-SQC.png: Python/aggregate_CNVkit_genome.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-genome-ADC.pdf Output/CNVkit/Bowtie2-genome-ADC.png: Python/aggregate_CNVkit_genome.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-genome-SQC.patient.pdf Output/CNVkit/BWA-genome-SQC.patient.png: Python/aggregate_CNVkit_genome.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-genome-ADC.patient.pdf Output/CNVkit/BWA-genome-ADC.patient.png: Python/aggregate_CNVkit_genome.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-genome-SQC.patient.pdf Output/CNVkit/Bowtie2-genome-SQC.patient.png: Python/aggregate_CNVkit_genome.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-genome-ADC.patient.pdf Output/CNVkit/Bowtie2-genome-ADC.patient.png: Python/aggregate_CNVkit_genome.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-simple-SQC.pdf Output/CNVkit/BWA-simple-SQC.png: Python/aggregate_CNVkit_simple.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-simple-ADC.pdf Output/CNVkit/BWA-simple-ADC.png: Python/aggregate_CNVkit_simple.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-simple-SQC.pdf Output/CNVkit/Bowtie2-simple-SQC.png: Python/aggregate_CNVkit_simple.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-simple-ADC.pdf Output/CNVkit/Bowtie2-simple-ADC.png: Python/aggregate_CNVkit_simple.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-genome-SQC.Recurrence.pdf Output/CNVkit/BWA-genome-SQC.Recurrence.png: Python/aggregate_CNVkit_genome_clinical.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-genome-ADC.Recurrence.pdf Output/CNVkit/BWA-genome-ADC.Recurrence.png: Python/aggregate_CNVkit_genome_clinical.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-genome-SQC.Recurrence.pdf Output/CNVkit/Bowtie2-genome-SQC.Recurrence.png: Python/aggregate_CNVkit_genome_clinical.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-genome-ADC.Recurrence.pdf Output/CNVkit/Bowtie2-genome-ADC.Recurrence.png: Python/aggregate_CNVkit_genome_clinical.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-genome-SQC.Recurrence.patient.pdf Output/CNVkit/BWA-genome-SQC.Recurrence.patient.png: Python/aggregate_CNVkit_genome_clinical.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-genome-ADC.Recurrence.patient.pdf Output/CNVkit/BWA-genome-ADC.Recurrence.patient.png: Python/aggregate_CNVkit_genome_clinical.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-genome-SQC.Recurrence.patient.pdf Output/CNVkit/Bowtie2-genome-SQC.Recurrence.patient.png: Python/aggregate_CNVkit_genome_clinical.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-genome-ADC.Recurrence.patient.pdf Output/CNVkit/Bowtie2-genome-ADC.Recurrence.patient.png: Python/aggregate_CNVkit_genome_clinical.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-genome-SQC.Smoking.pdf Output/CNVkit/BWA-genome-SQC.Smoking.png: Python/aggregate_CNVkit_genome_clinical.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-genome-ADC.Smoking.pdf Output/CNVkit/BWA-genome-ADC.Smoking.png: Python/aggregate_CNVkit_genome_clinical.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-genome-SQC.Smoking.pdf Output/CNVkit/Bowtie2-genome-SQC.Smoking.png: Python/aggregate_CNVkit_genome_clinical.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC --type 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-genome-ADC.Smoking.pdf Output/CNVkit/Bowtie2-genome-ADC.Smoking.png: Python/aggregate_CNVkit_genome_clinical.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC --type 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-genome-SQC.Smoking.patient.pdf Output/CNVkit/BWA-genome-SQC.Smoking.patient.png: Python/aggregate_CNVkit_genome_clinical.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-genome-ADC.Smoking.patient.pdf Output/CNVkit/BWA-genome-ADC.Smoking.patient.png: Python/aggregate_CNVkit_genome_clinical.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-genome-SQC.Smoking.patient.pdf Output/CNVkit/Bowtie2-genome-SQC.Smoking.patient.png: Python/aggregate_CNVkit_genome_clinical.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC --patient 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-genome-ADC.Smoking.patient.pdf Output/CNVkit/Bowtie2-genome-ADC.Smoking.patient.png: Python/aggregate_CNVkit_genome_clinical.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC --patient 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-violin-SQC.pdf: Python/aggregate_CNVkit_violin.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-violin-ADC.pdf: Python/aggregate_CNVkit_violin.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-violin-SQC.pdf: Python/aggregate_CNVkit_violin.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-violin-ADC.pdf: Python/aggregate_CNVkit_violin.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-violin-SQC.Recurrence.pdf: Python/aggregate_CNVkit_violin_clinical.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-violin-ADC.Recurrence.pdf: Python/aggregate_CNVkit_violin_clinical.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-violin-SQC.Recurrence.pdf: Python/aggregate_CNVkit_violin_clinical.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-violin-ADC.Recurrence.pdf: Python/aggregate_CNVkit_violin_clinical.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-violin-SQC.Smoking.pdf: Python/aggregate_CNVkit_violin_clinical.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-violin-ADC.Smoking.pdf: Python/aggregate_CNVkit_violin_clinical.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-violin-SQC.Smoking.pdf: Python/aggregate_CNVkit_violin_clinical.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-violin-ADC.Smoking.pdf: Python/aggregate_CNVkit_violin_clinical.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --compare "Smoking" "Never" "Ever" --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-chromosome-SQC.tar: Python/aggregate_CNVkit_chromosome.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(GENE_CENCUS_ONE_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/CNVkit/BWA-chromosome-ADC.tar: Python/aggregate_CNVkit_chromosome.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(GENE_CENCUS_ONE_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-chromosome-SQC.tar: Python/aggregate_CNVkit_chromosome.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(GENE_CENCUS_ONE_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/CNVkit/Bowtie2-chromosome-ADC.tar: Python/aggregate_CNVkit_chromosome.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_SIZE_TSV) $(GENE_CENCUS_ONE_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

step39-Reference: Output/CNVkit/BWA/Reference.tsv Output/CNVkit/Bowtie2/Reference.tsv
step39-CNVkit: $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID))) $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cnr,$(WES_REAL_TUMOR_ID)))
step39-genome: Output/CNVkit/BWA-genome-SQC.pdf Output/CNVkit/BWA-genome-ADC.pdf Output/CNVkit/Bowtie2-genome-SQC.pdf Output/CNVkit/Bowtie2-genome-ADC.pdf Output/CNVkit/BWA-genome-SQC.patient.pdf Output/CNVkit/BWA-genome-ADC.patient.pdf Output/CNVkit/Bowtie2-genome-SQC.patient.pdf Output/CNVkit/Bowtie2-genome-ADC.patient.pdf
step39-genome-clinical: Output/CNVkit/BWA-genome-SQC.Recurrence.pdf Output/CNVkit/BWA-genome-ADC.Recurrence.pdf Output/CNVkit/Bowtie2-genome-SQC.Recurrence.pdf Output/CNVkit/Bowtie2-genome-ADC.Recurrence.pdf Output/CNVkit/BWA-genome-SQC.Recurrence.patient.pdf Output/CNVkit/BWA-genome-ADC.Recurrence.patient.pdf Output/CNVkit/Bowtie2-genome-SQC.Recurrence.patient.pdf Output/CNVkit/Bowtie2-genome-ADC.Recurrence.patient.pdf Output/CNVkit/BWA-genome-SQC.Smoking.pdf Output/CNVkit/BWA-genome-ADC.Smoking.pdf Output/CNVkit/Bowtie2-genome-SQC.Smoking.pdf Output/CNVkit/Bowtie2-genome-ADC.Smoking.pdf Output/CNVkit/BWA-genome-SQC.Smoking.patient.pdf Output/CNVkit/BWA-genome-ADC.Smoking.patient.pdf Output/CNVkit/Bowtie2-genome-SQC.Smoking.patient.pdf Output/CNVkit/Bowtie2-genome-ADC.Smoking.patient.pdf
step39-simple: Output/CNVkit/BWA-simple-SQC.pdf Output/CNVkit/BWA-simple-ADC.pdf Output/CNVkit/Bowtie2-simple-SQC.pdf Output/CNVkit/Bowtie2-simple-ADC.pdf
step39-violin: Output/CNVkit/BWA-violin-SQC.pdf Output/CNVkit/BWA-violin-ADC.pdf Output/CNVkit/Bowtie2-violin-SQC.pdf Output/CNVkit/Bowtie2-violin-ADC.pdf
step39-violin-clinical: Output/CNVkit/BWA-violin-SQC.Recurrence.pdf Output/CNVkit/BWA-violin-ADC.Recurrence.pdf Output/CNVkit/Bowtie2-violin-ADC.Recurrence.pdf Output/CNVkit/Bowtie2-violin-ADC.Recurrence.pdf Output/CNVkit/BWA-violin-SQC.Smoking.pdf Output/CNVkit/BWA-violin-ADC.Smoking.pdf Output/CNVkit/Bowtie2-violin-SQC.Smoking.pdf Output/CNVkit/Bowtie2-violin-ADC.Smoking.pdf
step39-chromosome: Output/CNVkit/BWA-chromosome-SQC.tar Output/CNVkit/BWA-chromosome-ADC.tar Output/CNVkit/Bowtie2-chromosome-SQC.tar Output/CNVkit/Bowtie2-chromosome-ADC.tar
step39: step39-Reference step39-CNVkit step39-genome step39-genome-clinical step39-simple step39-violin step39-violin-clinical step39-chromosome
.PHONY: step39-Reference step39-CNVkit step39-genome step39-genome-clinical step39-simple step39-violin step39-violin-clinical step39-chromosome step39

# Step 40 (Gene Location)
Output/Gene_Location:
	mkdir -p $@

Output/Gene_Location/gene.tar: Python/get_gene_location.py $(GENE_CENCUS_ONE_CSV) $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log Output/Gene_Location
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Gene_Location/cytoband.tar: Python/get_gene_location_2.py $(GENE_CENCUS_ONE_CSV) $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log Output/Gene_Location
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step40: Output/Gene_Location/gene.tar Output/Gene_Location/cytoband.tar
.PHONY: step40

# Step 41 (Mobster)
Output/Mobster/BWA Output/Mobster/Bowtie2:
	mkdir -p $@

Output/Mobster/BWA/SQC.txt: Python/make_input_mobster.py $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.bam,$(WES_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Mobster/BWA
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Mobster/BWA/ADC.txt: Python/make_input_mobster.py $(addprefix Output/WES/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.bam,$(WES_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Mobster/BWA
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Mobster/BWA/%_predictions.txt: misc/Mobster.properties Output/Mobster/BWA/%.txt | build.log
	$(WORKDIR_DOCKER) mobster -Xmx$(MEMS)g -properties $(addprefix /,$(filter %.properties,$^)) -in $$(cat $(filter %.txt,$^) | head -n 1) -out $* -sn $$(cat $(filter %.txt,$^) | head -n 2 | tail -n 1) 1> $@.stdout 2> $@.stderr

Output/Mobster/Bowtie2/SQC.txt: Python/make_input_mobster.py $(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.bam,$(WES_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Mobster/Bowtie2
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Mobster/Bowtie2/ADC.txt: Python/make_input_mobster.py $(addprefix Output/WES/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.bam,$(WES_ID))) $(CLINICAL_DATA_CSV) | build.log Output/Mobster/Bowtie2
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Mobster/Bowtie2/%_predictions.txt: misc/Mobster.properties Output/Mobster/Bowtie2/%.txt | build.log
	$(WORKDIR_DOCKER) mobster -Xmx$(MEMS)g -properties $(addprefix /,$(filter %.properties,$^)) -in $$(cat $(filter %.txt,$^) | head -n 1) -out $* -sn $$(cat $(filter %.txt,$^) | head -n 2 | tail -n 1) 1> $@.stdout 2> $@.stderr

step41-BWA: Output/Mobster/BWA/SQC.txt Output/Mobster/BWA/ADC.txt Output/Mobster/BWA/SQC_predictions.txt Output/Mobster/BWA/ADC_predictions.txt
step41-Bowtie2: Output/Mobster/Bowtie2/SQC.txt Output/Mobster/Bowtie2/ADC.txt Output/Mobster/Bowtie2/SQC_predictions.txt Output/Mobster/Bowtie2/ADC_predictions.txt
step41: step41-BWA step41-Bowtie2
.PHONY: step41-BWA step41-Bowtie2 step41

# Step 42 (CollectRnaSeqMetrics)
Output/CollectRnaSeqMetrics/STAR Output/CollectRnaSeqMetrics/Bowtie2:
	mkdir -p $@

Output/CollectRnaSeqMetrics/STAR/%.RNA_metrics: $(UCSC_REFFLAT_TXT) Output/RSEM/Bowtie2/%.genome.bam | Output/CollectRnaSeqMetrics/STAR build.log
	$(DOCKER) $(PICARD) CollectRnaSeqMetrics --INPUT $(addprefix /,$(filter %.bam,$^)) --OUTPUT $(addprefix /,$@) --REF_FLAT $(addprefix /,$(filter %.txt,$^)) --STRAND_SPECIFICITY NONE 1> $@.stdout 2> $@.stderr

Output/CollectRnaSeqMetrics/Bowtie2/%.RNA_metrics: $(UCSC_REFFLAT_TXT) Output/RSEM/Bowtie2/%.genome.bam | Output/CollectRnaSeqMetrics/Bowtie2 build.log
	$(DOCKER) $(PICARD) CollectRnaSeqMetrics --INPUT $(addprefix /,$(filter %.bam,$^)) --OUTPUT $(addprefix /,$@) --REF_FLAT $(addprefix /,$(filter %.txt,$^)) --STRAND_SPECIFICITY NONE 1> $@.stdout 2> $@.stderr

step42-STAR: $(addprefix Output/CollectRnaSeqMetrics/STAR/,$(addsuffix .RNA_metrics,$(WTS_ID)))
step42-Bowtie2: $(addprefix Output/CollectRnaSeqMetrics/Bowtie2/,$(addsuffix .RNA_metrics,$(WTS_ID)))
step42: step42-STAR step42-Bowtie2
.PHONY: step42-STAR step42-Bowtie2 step42

# Step 43 (Gistic2.0)
# http://portals.broadinstitute.org/cgi-bin/cancer/publications/pub_paper.cgi?mode=view&paper_id=216&p=t
# ftp://ftp.broadinstitute.org/pub/GISTIC2.0/INSTALL.txt
Output/Gistic:
	mkdir -p $@

Output/Gistic/hg38.refgene.mat: | Output/Gistic
	wget "ftp://ftp.broadinstitute.org/distribution/GISTIC2.0/refgenes/hg38.UCSC.add_miR.160920.refgene.mat" -O $@
	touch $@

Output/Gistic/BWA-CNVkit.SQC.tsv: Python/merge_CNVkit_Gistic.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC.tsv: Python/merge_CNVkit_Gistic.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC.tsv: Python/merge_CNVkit_Gistic.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC.tsv: Python/merge_CNVkit_Gistic.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Dysplasia.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Dysplasia" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-AAH.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "AAH" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-CIS.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "CIS+AIS" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Precancer.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Primary.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-AAH.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "AAH" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-AIS.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "CIS+AIS" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Precancer.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Primary.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Dysplasia.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Dysplasia" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-AAH.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "AAH" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-CIS.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "CIS+AIS" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Primary.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-AAH.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "AAH" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-AIS.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "CIS+AIS" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Primary.tsv: Python/merge_CNVkit_Gistic_2.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-CIS.Nonrecur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "CIS+AIS" --exact "Recurrence" "NO" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-CIS.Recur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "CIS+AIS" --exact "Recurrence" "YES" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Precancer.Nonrecur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Recurrence" "NO" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Precancer.Recur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Recurrence" "YES" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Primary.Nonrecur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Primary" --exact "Recurrence" "NO" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Primary.Recur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Primary" --exact "Recurrence" "YES" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Nonrecur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "CIS+AIS" --exact "Recurrence" "NO" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Recur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "CIS+AIS" --exact "Recurrence" "YES" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Nonrecur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Recurrence" "NO" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Recur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Recurrence" "YES" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Nonrecur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Primary" --exact "Recurrence" "NO" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Recur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Primary" --exact "Recurrence" "YES" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-AAH.Nonrecur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "AAH" --exact "Recurrence" "NO" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-AAH.Recur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "AAH" --exact "Recurrence" "YES" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Precancer.Nonrecur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Recurrence" "NO" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Precancer.Recur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Recurrence" "YES" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Primary.Nonrecur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Primary" --exact "Recurrence" "NO" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Primary.Recur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Primary" --exact "Recurrence" "YES" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Nonrecur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "AAH" --exact "Recurrence" "NO" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Recur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "AAH" --exact "Recurrence" "YES" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Nonrecur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Recurrence" "NO" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Recur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Recurrence" "YES" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Nonrecur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Primary" --exact "Recurrence" "NO" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Recur.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Primary" --exact "Recurrence" "YES" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Never.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "CIS+AIS" --exact "Smoking-Detail" "Never" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Ex.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "CIS+AIS" --exact "Smoking-Detail" "Ex" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Current.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "CIS+AIS" --exact "Smoking-Detail" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Never.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Smoking-Detail" "Never" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Ex.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Smoking-Detail" "Ex" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Current.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Smoking-Detail" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Never.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Primary" --exact "Smoking-Detail" "Never" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Ex.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Primary" --exact "Smoking-Detail" "Ex" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Current.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Primary" --exact "Smoking-Detail" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking_Never.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "CIS+AIS" --exact "Smoking-Detail" "Never" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking_Ex.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "CIS+AIS" --exact "Smoking-Detail" "Ex" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking_Current.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "CIS+AIS" --exact "Smoking-Detail" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking_Never.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Smoking-Detail" "Never" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking_Ex.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Smoking-Detail" "Ex" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking_Current.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Smoking-Detail" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking_Never.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Primary" --exact "Smoking-Detail" "Never" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking_Ex.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Primary" --exact "Smoking-Detail" "Ex" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking_Current.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --stage "Primary" --exact "Smoking-Detail" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking_Ex.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "AAH" --exact "Smoking-Detail" "Ex" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking_Current.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "AAH" --exact "Smoking-Detail" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking_Ex.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Smoking-Detail" "Ex" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking_Current.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Smoking-Detail" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking_Ex.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Primary" --exact "Smoking-Detail" "Ex" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking_Current.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/BWA/,$(addsuffix .BWA.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Primary" --exact "Smoking-Detail" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking_Ex.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "AAH" --exact "Smoking-Detail" "Ex" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking_Current.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "AAH" --exact "Smoking-Detail" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking_Ex.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Smoking-Detail" "Ex" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking_Current.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Dysplasia" "AAH" "CIS+AIS" "MIA" --exact "Smoking-Detail" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking_Ex.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Primary" --exact "Smoking-Detail" "Ex" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking_Current.tsv: Python/merge_CNVkit_Gistic_3.py $(addprefix Output/CNVkit/Bowtie2/,$(addsuffix .Bowtie2.Sort.MarkDuplicate.BQSR.cns,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) $(REFERENCE_BLACKLIST) | Output/Gistic build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --stage "Primary" --exact "Smoking-Detail" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/%/all_data_by_genes.txt Output/Gistic/%/all_thresholded.by_genes.txt Output/Gistic/%/focal_data_by_genes.txt Output/Gistic/%/amp_genes.conf_99.txt Output/Gistic/%/del_genes.conf_99.txt Output/Gistic/%/all_lesions.conf_99.txt: Output/Gistic/hg38.refgene.mat Output/Gistic/%.tsv | build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D)) && mkdir $(@D)
	$(TOOLS_DOCKER) ./gistic2 -b $(addprefix /,$(@D)) -seg $(addprefix /,$(filter %.tsv,$^)) -refgene $(addprefix /,$(filter %.mat,$^)) -genegistic 1 -smallmem 0 -broad 1 -brlen 0.7 -conf 0.99 -armpeel 1 -savegene 1 -gcm mean -maxseg 2500 -ta 0.2 -td 0.2 -cap 1.5 -js 2 -rx 0 -v 10 1> $(@D).stdout 2> $(@D).stderr

Output/Gistic/BWA-CNVkit.SQC.Amplification.tex Output/Gistic/BWA-CNVkit.SQC.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.SQC-Dysplasia/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Dysplasia" "CIS" "Primary" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC.Deletion.tex Output/Gistic/BWA-CNVkit.SQC.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.SQC-Dysplasia/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Dysplasia" "CIS" "Primary" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-simple.Amplification.tex Output/Gistic/BWA-CNVkit.SQC-simple.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.SQC-Precancer/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Precancer" "Primary" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-simple.Deletion.tex Output/Gistic/BWA-CNVkit.SQC-simple.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.SQC-Precancer/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Precancer" "Primary" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC.Amplification.tex Output/Gistic/BWA-CNVkit.ADC.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.ADC-AAH/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AIS/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "AAH" "AIS" "Primary" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC.Deletion.tex Output/Gistic/BWA-CNVkit.ADC.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.ADC-AAH/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AIS/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "AAH" "AIS" "Primary" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-simple.Amplification.tex Output/Gistic/BWA-CNVkit.ADC-simple.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.ADC-Precancer/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Precancer" "Primary" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC.Deletion-simple.tex Output/Gistic/BWA-CNVkit.ADC-simple.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.ADC-Precancer/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Precancer" "Primary" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC.Amplification.tex Output/Gistic/Bowtie2-CNVkit.SQC.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.SQC-Dysplasia/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Dysplasia" "CIS" "Primary" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC.Deletion.tex Output/Gistic/Bowtie2-CNVkit.SQC.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.SQC-Dysplasia/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Dysplasia" "CIS" "Primary" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-simple.Amplification.tex Output/Gistic/Bowtie2-CNVkit.SQC-simple.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.SQC-Precancer/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Precancer" "Primary" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-simple.Deletion.tex Output/Gistic/Bowtie2-CNVkit.SQC-simple.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.SQC-Precancer/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Precancer" "Primary" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC.Amplification.tex Output/Gistic/Bowtie2-CNVkit.ADC.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.ADC-AAH/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AIS/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "AAH" "AIS" "Primary" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC.Deletion.tex Output/Gistic/Bowtie2-CNVkit.ADC.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.ADC-AAH/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AIS/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "AAH" "AIS" "Primary" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-simple.Amplification.tex Output/Gistic/Bowtie2-CNVkit.ADC-simple.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.ADC-Precancer/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Precancer" "Primary" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-simple.Deletion.tex Output/Gistic/Bowtie2-CNVkit.ADC-simple.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.ADC-Precancer/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Precancer" "Primary" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-CIS.NR-R.Amplification.tex Output/Gistic/BWA-CNVkit.SQC-CIS.NR-R.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.SQC-CIS.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-CIS.NR-R.Deletion.tex Output/Gistic/BWA-CNVkit.SQC-CIS.NR-R.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.SQC-CIS.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Precancer.NR-R.Amplification.tex Output/Gistic/BWA-CNVkit.SQC-Precancer.NR-R.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.SQC-Precancer.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Precancer.NR-R.Deletion.tex Output/Gistic/BWA-CNVkit.SQC-Precancer.NR-R.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.SQC-Precancer.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Primary.NR-R.Amplification.tex Output/Gistic/BWA-CNVkit.SQC-Primary.NR-R.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.SQC-Primary.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Primary.NR-R.Deletion.tex Output/Gistic/BWA-CNVkit.SQC-Primary.NR-R.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.SQC-Primary.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-CIS.NR-R.Amplification.tex Output/Gistic/Bowtie2-CNVkit.SQC-CIS.NR-R.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-CIS.NR-R.Deletion.tex Output/Gistic/Bowtie2-CNVkit.SQC-CIS.NR-R.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.NR-R.Amplification.tex Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.NR-R.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.NR-R.Deletion.tex Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.NR-R.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Primary.NR-R.Amplification.tex Output/Gistic/Bowtie2-CNVkit.SQC-Primary.NR-R.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Primary.NR-R.Deletion.tex Output/Gistic/Bowtie2-CNVkit.SQC-Primary.NR-R.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-AAH.NR-R.Amplification.tex Output/Gistic/BWA-CNVkit.ADC-AAH.NR-R.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.ADC-AAH.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AAH.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-AAH.NR-R.Deletion.tex Output/Gistic/BWA-CNVkit.ADC-AAH.NR-R.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.ADC-AAH.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AAH.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Precancer.NR-R.Amplification.tex Output/Gistic/BWA-CNVkit.ADC-Precancer.NR-R.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.ADC-Precancer.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Precancer.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Precancer.NR-R.Deletion.tex Output/Gistic/BWA-CNVkit.ADC-Precancer.NR-R.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.ADC-Precancer.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Precancer.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Primary.NR-R.Amplification.tex Output/Gistic/BWA-CNVkit.ADC-Primary.NR-R.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.ADC-Primary.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Primary.NR-R.Deletion.tex Output/Gistic/BWA-CNVkit.ADC-Primary.NR-R.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.ADC-Primary.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-AAH.NR-R.Amplification.tex Output/Gistic/Bowtie2-CNVkit.ADC-AAH.NR-R.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-AAH.NR-R.Deletion.tex Output/Gistic/Bowtie2-CNVkit.ADC-AAH.NR-R.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.NR-R.Amplification.tex Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.NR-R.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.NR-R.Deletion.tex Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.NR-R.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Primary.NR-R.Amplification.tex Output/Gistic/Bowtie2-CNVkit.ADC-Primary.NR-R.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Primary.NR-R.Deletion.tex Output/Gistic/Bowtie2-CNVkit.ADC-Primary.NR-R.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Recur/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking.Amplification.tex Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Never/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking.Deletion.tex Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Never/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking.Amplification.tex Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Never/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking.Deletion.tex Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Never/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking.Amplification.tex Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Never/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking.Deletion.tex Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Never/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking.Amplification.tex Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking_Never/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking.Deletion.tex Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking_Never/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking.Amplification.tex Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking_Never/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking.Deletion.tex Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking_Never/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking.Amplification.tex Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking_Never/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking.Deletion.tex Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking_Never/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking.Amplification.tex Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking.Deletion.tex Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking.Amplification.tex Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking.Deletion.tex Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking.Amplification.tex Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking.Deletion.tex Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking.Amplification.tex Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking.Deletion.tex Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking.Amplification.tex Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking.Deletion.tex Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking.Amplification.tex Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking.Amplification.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" --amplification 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking.Deletion.tex Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking.Deletion.pdf: Python/draw_venn_plot_gistic.py Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking_Current/all_lesions.conf_99.txt $(GENE_CENCUS_BOTH_CSV) $(REFERENCE_BAND) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" --deletion 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC.Amplification.genes.tex Output/Gistic/BWA-CNVkit.SQC.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.SQC-Dysplasia/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Dysplasia" "CIS" "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC.Deletion.genes.tex Output/Gistic/BWA-CNVkit.SQC.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.SQC-Dysplasia/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Dysplasia" "CIS" "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-simple.Amplification.genes.tex Output/Gistic/BWA-CNVkit.SQC-simple.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.SQC-Precancer/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Precancer" "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-simple.Deletion.genes.tex Output/Gistic/BWA-CNVkit.SQC-simple.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.SQC-Precancer/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Precancer" "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC.Amplification.genes.tex Output/Gistic/BWA-CNVkit.ADC.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.ADC-AAH/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AIS/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "AAH" "AIS" "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC.Deletion.genes.tex Output/Gistic/BWA-CNVkit.ADC.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.ADC-AAH/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AIS/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "AAH" "AIS" "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-simple.Amplification.genes.tex Output/Gistic/BWA-CNVkit.ADC-simple.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.ADC-Precancer/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Precancer" "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC.Deletion-simple.genes.tex Output/Gistic/BWA-CNVkit.ADC-simple.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.ADC-Precancer/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Precancer" "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC.Amplification.genes.tex Output/Gistic/Bowtie2-CNVkit.SQC.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.SQC-Dysplasia/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Dysplasia" "CIS" "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC.Deletion.genes.tex Output/Gistic/Bowtie2-CNVkit.SQC.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.SQC-Dysplasia/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Dysplasia" "CIS" "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-simple.Amplification.genes.tex Output/Gistic/Bowtie2-CNVkit.SQC-simple.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.SQC-Precancer/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Precancer" "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-simple.Deletion.genes.tex Output/Gistic/Bowtie2-CNVkit.SQC-simple.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.SQC-Precancer/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Precancer" "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC.Amplification.genes.tex Output/Gistic/Bowtie2-CNVkit.ADC.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.ADC-AAH/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AIS/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "AAH" "AIS" "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC.Deletion.genes.tex Output/Gistic/Bowtie2-CNVkit.ADC.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.ADC-AAH/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AIS/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "AAH" "AIS" "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-simple.Amplification.genes.tex Output/Gistic/Bowtie2-CNVkit.ADC-simple.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.ADC-Precancer/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Precancer" "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-simple.Deletion.genes.tex Output/Gistic/Bowtie2-CNVkit.ADC-simple.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.ADC-Precancer/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Precancer" "Primary" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-CIS.NR-R.Amplification.genes.tex Output/Gistic/BWA-CNVkit.SQC-CIS.NR-R.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.SQC-CIS.Nonrecur/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.Recur/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-CIS.NR-R.Deletion.genes.tex Output/Gistic/BWA-CNVkit.SQC-CIS.NR-R.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.SQC-CIS.Nonrecur/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.Recur/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Precancer.NR-R.Amplification.genes.tex Output/Gistic/BWA-CNVkit.SQC-Precancer.NR-R.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.SQC-Precancer.Nonrecur/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Recur/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Precancer.NR-R.Deletion.genes.tex Output/Gistic/BWA-CNVkit.SQC-Precancer.NR-R.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.SQC-Precancer.Nonrecur/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Recur/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Primary.NR-R.Amplification.genes.tex Output/Gistic/BWA-CNVkit.SQC-Primary.NR-R.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.SQC-Primary.Nonrecur/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Recur/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Primary.NR-R.Deletion.genes.tex Output/Gistic/BWA-CNVkit.SQC-Primary.NR-R.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.SQC-Primary.Nonrecur/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Recur/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-CIS.NR-R.Amplification.genes.tex Output/Gistic/Bowtie2-CNVkit.SQC-CIS.NR-R.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Nonrecur/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Recur/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-CIS.NR-R.Deletion.genes.tex Output/Gistic/Bowtie2-CNVkit.SQC-CIS.NR-R.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Nonrecur/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Recur/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.NR-R.Amplification.genes.tex Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.NR-R.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Nonrecur/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Recur/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.NR-R.Deletion.genes.tex Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.NR-R.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Nonrecur/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Recur/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Primary.NR-R.Amplification.genes.tex Output/Gistic/Bowtie2-CNVkit.SQC-Primary.NR-R.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Nonrecur/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Recur/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Primary.NR-R.Deletion.genes.tex Output/Gistic/Bowtie2-CNVkit.SQC-Primary.NR-R.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Nonrecur/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Recur/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-AAH.NR-R.Amplification.genes.tex Output/Gistic/BWA-CNVkit.ADC-AAH.NR-R.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.ADC-AAH.Nonrecur/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AAH.Recur/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-AAH.NR-R.Deletion.genes.tex Output/Gistic/BWA-CNVkit.ADC-AAH.NR-R.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.ADC-AAH.Nonrecur/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AAH.Recur/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Precancer.NR-R.Amplification.genes.tex Output/Gistic/BWA-CNVkit.ADC-Precancer.NR-R.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.ADC-Precancer.Nonrecur/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Precancer.Recur/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Precancer.NR-R.Deletion.genes.tex Output/Gistic/BWA-CNVkit.ADC-Precancer.NR-R.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.ADC-Precancer.Nonrecur/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Precancer.Recur/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Primary.NR-R.Amplification.genes.tex Output/Gistic/BWA-CNVkit.ADC-Primary.NR-R.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.ADC-Primary.Nonrecur/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary.Recur/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Primary.NR-R.Deletion.genes.tex Output/Gistic/BWA-CNVkit.ADC-Primary.NR-R.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.ADC-Primary.Nonrecur/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary.Recur/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-AAH.NR-R.Amplification.genes.tex Output/Gistic/Bowtie2-CNVkit.ADC-AAH.NR-R.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Nonrecur/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Recur/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-AAH.NR-R.Deletion.genes.tex Output/Gistic/Bowtie2-CNVkit.ADC-AAH.NR-R.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Nonrecur/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Recur/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.NR-R.Amplification.genes.tex Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.NR-R.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Nonrecur/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Recur/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.NR-R.Deletion.genes.tex Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.NR-R.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Nonrecur/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Recur/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Primary.NR-R.Amplification.genes.tex Output/Gistic/Bowtie2-CNVkit.ADC-Primary.NR-R.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Nonrecur/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Recur/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Primary.NR-R.Deletion.genes.tex Output/Gistic/Bowtie2-CNVkit.ADC-Primary.NR-R.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Nonrecur/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Recur/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Non-recurrence" "Recurrence" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking.Amplification.genes.tex Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Never/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Ex/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Current/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking.Deletion.genes.tex Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Never/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Ex/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Current/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking.Amplification.genes.tex Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Never/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Ex/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Current/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking.Deletion.genes.tex Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Never/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Ex/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Current/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking.Amplification.genes.tex Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Never/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Ex/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Current/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking.Deletion.genes.tex Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Never/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Ex/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Current/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking.Amplification.genes.tex Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking_Never/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking_Ex/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking_Current/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking.Deletion.genes.tex Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking_Never/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking_Ex/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Smoking_Current/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking.Amplification.genes.tex Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking_Never/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking_Ex/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking_Current/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking.Deletion.genes.tex Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking_Never/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking_Ex/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Smoking_Current/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking.Amplification.genes.tex Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking_Never/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking_Ex/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking_Current/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking.Deletion.genes.tex Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking_Never/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking_Ex/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Smoking_Current/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking.Amplification.genes.tex Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking_Ex/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking_Current/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking.Deletion.genes.tex Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking_Ex/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking_Current/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking.Amplification.genes.tex Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking_Ex/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking_Current/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking.Deletion.genes.tex Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking_Ex/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking_Current/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking.Amplification.genes.tex Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking_Ex/amp_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking_Current/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking.Deletion.genes.tex Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking_Ex/del_genes.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking_Current/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking.Amplification.genes.tex Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking_Ex/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking_Current/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking.Deletion.genes.tex Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking_Ex/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking_Current/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking.Amplification.genes.tex Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking_Ex/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking_Current/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking.Deletion.genes.tex Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking_Ex/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking_Current/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking.Amplification.genes.tex Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking.Amplification.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking_Ex/amp_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking_Current/amp_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking.Deletion.genes.tex Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking.Deletion.genes.pdf: Python/compare_gistic_peaks.py Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking_Ex/del_genes.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking_Current/del_genes.conf_99.txt $(GENE_CENCUS_BOTH_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@)) --annotation "Ex" "Current" 1> $@.stdout 2> $@.stderr

step43-SQC-stage-BWA: Output/Gistic/BWA-CNVkit.SQC.tsv Output/Gistic/BWA-CNVkit.SQC/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Dysplasia.tsv Output/Gistic/BWA-CNVkit.SQC-Dysplasia/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-AAH.tsv Output/Gistic/BWA-CNVkit.SQC-AAH/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.tsv Output/Gistic/BWA-CNVkit.SQC-CIS/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.tsv Output/Gistic/BWA-CNVkit.SQC-Primary/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.tsv Output/Gistic/BWA-CNVkit.SQC-Precancer/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC.Amplification.pdf Output/Gistic/BWA-CNVkit.SQC.Deletion.pdf Output/Gistic/BWA-CNVkit.SQC-simple.Amplification.pdf Output/Gistic/BWA-CNVkit.SQC-simple.Deletion.pdf Output/Gistic/BWA-CNVkit.SQC.Amplification.genes.pdf Output/Gistic/BWA-CNVkit.SQC.Deletion.genes.pdf Output/Gistic/BWA-CNVkit.SQC-simple.Amplification.genes.pdf Output/Gistic/BWA-CNVkit.SQC-simple.Deletion.genes.pdf
step43-SQC-stage-Bowtie2: Output/Gistic/Bowtie2-CNVkit.SQC.tsv Output/Gistic/Bowtie2-CNVkit.SQC/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Dysplasia.tsv Output/Gistic/Bowtie2-CNVkit.SQC-Dysplasia/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-AAH.tsv Output/Gistic/Bowtie2-CNVkit.SQC-AAH/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS.tsv Output/Gistic/Bowtie2-CNVkit.SQC-CIS/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary.tsv Output/Gistic/Bowtie2-CNVkit.SQC-Primary/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.tsv Output/Gistic/Bowtie2-CNVkit.SQC-Precancer/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC.Amplification.pdf Output/Gistic/Bowtie2-CNVkit.SQC.Deletion.pdf Output/Gistic/Bowtie2-CNVkit.SQC-simple.Amplification.pdf Output/Gistic/Bowtie2-CNVkit.SQC-simple.Deletion.pdf Output/Gistic/Bowtie2-CNVkit.SQC.Amplification.genes.pdf Output/Gistic/Bowtie2-CNVkit.SQC.Deletion.genes.pdf Output/Gistic/Bowtie2-CNVkit.SQC-simple.Amplification.genes.pdf Output/Gistic/Bowtie2-CNVkit.SQC-simple.Deletion.genes.pdf
step43-ADC-stage-BWA: Output/Gistic/BWA-CNVkit.ADC.tsv Output/Gistic/BWA-CNVkit.ADC/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AAH.tsv Output/Gistic/BWA-CNVkit.ADC-AAH/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AIS.tsv Output/Gistic/BWA-CNVkit.ADC-AIS/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary.tsv Output/Gistic/BWA-CNVkit.ADC-Primary/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Precancer.tsv Output/Gistic/BWA-CNVkit.ADC-Precancer/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC.Amplification.pdf Output/Gistic/BWA-CNVkit.ADC.Deletion.pdf Output/Gistic/BWA-CNVkit.ADC-simple.Amplification.pdf Output/Gistic/BWA-CNVkit.ADC-simple.Deletion.pdf Output/Gistic/BWA-CNVkit.ADC.Amplification.genes.pdf Output/Gistic/BWA-CNVkit.ADC.Deletion.genes.pdf Output/Gistic/BWA-CNVkit.ADC-simple.Amplification.genes.pdf Output/Gistic/BWA-CNVkit.ADC-simple.Deletion.genes.pdf
step43-ADC-stage-Bowtie2: Output/Gistic/Bowtie2-CNVkit.ADC.tsv Output/Gistic/Bowtie2-CNVkit.ADC/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AAH.tsv Output/Gistic/Bowtie2-CNVkit.ADC-AAH/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AIS.tsv Output/Gistic/Bowtie2-CNVkit.ADC-AIS/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary.tsv Output/Gistic/Bowtie2-CNVkit.ADC-Primary/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.tsv Output/Gistic/Bowtie2-CNVkit.ADC-Precancer/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC.Amplification.pdf Output/Gistic/Bowtie2-CNVkit.ADC.Deletion.pdf Output/Gistic/Bowtie2-CNVkit.ADC-simple.Amplification.pdf Output/Gistic/Bowtie2-CNVkit.ADC-simple.Deletion.pdf Output/Gistic/Bowtie2-CNVkit.ADC.Amplification.genes.pdf Output/Gistic/Bowtie2-CNVkit.ADC.Deletion.genes.pdf Output/Gistic/Bowtie2-CNVkit.ADC-simple.Amplification.genes.pdf Output/Gistic/Bowtie2-CNVkit.ADC-simple.Deletion.genes.pdf
step43-SQC-Recur-BWA: Output/Gistic/BWA-CNVkit.SQC-CIS.Nonrecur.tsv Output/Gistic/BWA-CNVkit.SQC-CIS.Recur.tsv Output/Gistic/BWA-CNVkit.SQC-Precancer.Nonrecur.tsv Output/Gistic/BWA-CNVkit.SQC-Precancer.Recur.tsv Output/Gistic/BWA-CNVkit.SQC-Primary.Nonrecur.tsv Output/Gistic/BWA-CNVkit.SQC-Primary.Recur.tsv Output/Gistic/BWA-CNVkit.SQC-CIS.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.Recur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Recur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Recur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.NR-R.Amplification.pdf Output/Gistic/BWA-CNVkit.SQC-CIS.NR-R.Deletion.pdf Output/Gistic/BWA-CNVkit.SQC-Precancer.NR-R.Amplification.pdf Output/Gistic/BWA-CNVkit.SQC-Precancer.NR-R.Deletion.pdf Output/Gistic/BWA-CNVkit.SQC-Primary.NR-R.Amplification.pdf Output/Gistic/BWA-CNVkit.SQC-Primary.NR-R.Deletion.pdf Output/Gistic/BWA-CNVkit.SQC-CIS.NR-R.Amplification.genes.pdf Output/Gistic/BWA-CNVkit.SQC-CIS.NR-R.Deletion.genes.pdf Output/Gistic/BWA-CNVkit.SQC-Precancer.NR-R.Amplification.genes.pdf Output/Gistic/BWA-CNVkit.SQC-Precancer.NR-R.Deletion.genes.pdf Output/Gistic/BWA-CNVkit.SQC-Primary.NR-R.Amplification.genes.pdf Output/Gistic/BWA-CNVkit.SQC-Primary.NR-R.Deletion.genes.pdf
step43-SQC-Recur-Bowtie2: Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Nonrecur.tsv Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Recur.tsv Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Nonrecur.tsv Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Recur.tsv Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Nonrecur.tsv Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Recur.tsv Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS.Recur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.Recur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-Primary.Recur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.SQC-CIS.NR-R.Amplification.pdf Output/Gistic/Bowtie2-CNVkit.SQC-CIS.NR-R.Deletion.pdf Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.NR-R.Amplification.pdf Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.NR-R.Deletion.pdf Output/Gistic/Bowtie2-CNVkit.SQC-Primary.NR-R.Amplification.pdf Output/Gistic/Bowtie2-CNVkit.SQC-Primary.NR-R.Deletion.pdf Output/Gistic/Bowtie2-CNVkit.SQC-CIS.NR-R.Amplification.genes.pdf Output/Gistic/Bowtie2-CNVkit.SQC-CIS.NR-R.Deletion.genes.pdf Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.NR-R.Amplification.genes.pdf Output/Gistic/Bowtie2-CNVkit.SQC-Precancer.NR-R.Deletion.genes.pdf Output/Gistic/Bowtie2-CNVkit.SQC-Primary.NR-R.Amplification.genes.pdf Output/Gistic/Bowtie2-CNVkit.SQC-Primary.NR-R.Deletion.genes.pdf
step43-ADC-Recur-BWA: Output/Gistic/BWA-CNVkit.ADC-AAH.Nonrecur.tsv Output/Gistic/BWA-CNVkit.ADC-AAH.Recur.tsv Output/Gistic/BWA-CNVkit.ADC-Precancer.Nonrecur.tsv Output/Gistic/BWA-CNVkit.ADC-Precancer.Recur.tsv Output/Gistic/BWA-CNVkit.ADC-Primary.Nonrecur.tsv Output/Gistic/BWA-CNVkit.ADC-Primary.Recur.tsv Output/Gistic/BWA-CNVkit.ADC-AAH.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AAH.Recur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Precancer.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Precancer.Recur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary.Nonrecur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary.Recur/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AAH.NR-R.Amplification.pdf Output/Gistic/BWA-CNVkit.ADC-AAH.NR-R.Deletion.pdf Output/Gistic/BWA-CNVkit.ADC-Precancer.NR-R.Amplification.pdf Output/Gistic/BWA-CNVkit.ADC-Precancer.NR-R.Deletion.pdf Output/Gistic/BWA-CNVkit.ADC-Primary.NR-R.Amplification.pdf Output/Gistic/BWA-CNVkit.ADC-Primary.NR-R.Deletion.pdf Output/Gistic/BWA-CNVkit.ADC-AAH.NR-R.Amplification.genes.pdf Output/Gistic/BWA-CNVkit.ADC-AAH.NR-R.Deletion.genes.pdf Output/Gistic/BWA-CNVkit.ADC-Precancer.NR-R.Amplification.genes.pdf Output/Gistic/BWA-CNVkit.ADC-Precancer.NR-R.Deletion.genes.pdf Output/Gistic/BWA-CNVkit.ADC-Primary.NR-R.Amplification.genes.pdf Output/Gistic/BWA-CNVkit.ADC-Primary.NR-R.Deletion.genes.pdf
step43-ADC-Recur-Bowtie2: Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Nonrecur.tsv Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Recur.tsv Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Nonrecur.tsv Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Recur.tsv Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Nonrecur.tsv Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Recur.tsv Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Recur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Recur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Nonrecur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Recur/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AAH.NR-R.Amplification.pdf Output/Gistic/Bowtie2-CNVkit.ADC-AAH.NR-R.Deletion.pdf Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.NR-R.Amplification.pdf Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.NR-R.Deletion.pdf Output/Gistic/Bowtie2-CNVkit.ADC-Primary.NR-R.Amplification.pdf Output/Gistic/Bowtie2-CNVkit.ADC-Primary.NR-R.Deletion.pdf Output/Gistic/Bowtie2-CNVkit.ADC-AAH.NR-R.Amplification.genes.pdf Output/Gistic/Bowtie2-CNVkit.ADC-AAH.NR-R.Deletion.genes.pdf Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.NR-R.Amplification.genes.pdf Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.NR-R.Deletion.genes.pdf Output/Gistic/Bowtie2-CNVkit.ADC-Primary.NR-R.Amplification.genes.pdf Output/Gistic/Bowtie2-CNVkit.ADC-Primary.NR-R.Deletion.genes.pdf
step43-SQC-Smoking-BWA: Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Never.tsv Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Ex.tsv Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Current.tsv Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Never.tsv Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Ex.tsv Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Current.tsv Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Never.tsv Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Ex.tsv Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Current.tsv Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Never/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking_Current/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Never/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking_Current/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Never/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking_Current/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking.Amplification.pdf Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking.Deletion.pdf Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking.Amplification.pdf Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking.Deletion.pdf Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking.Amplification.pdf Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking.Deletion.pdf Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking.Amplification.genes.pdf Output/Gistic/BWA-CNVkit.SQC-CIS.Smoking.Deletion.genes.pdf Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking.Amplification.genes.pdf Output/Gistic/BWA-CNVkit.SQC-Precancer.Smoking.Deletion.genes.pdf Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking.Amplification.genes.pdf Output/Gistic/BWA-CNVkit.SQC-Primary.Smoking.Deletion.genes.pdf
step43-ADC-Smoking-BWA: Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking_Ex.tsv Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking_Current.tsv Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking_Ex.tsv Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking_Current.tsv Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking_Ex.tsv Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking_Current.tsv Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking_Current/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking_Current/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking_Current/all_lesions.conf_99.txt Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking.Amplification.pdf Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking.Deletion.pdf Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking.Amplification.pdf Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking.Deletion.pdf Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking.Amplification.pdf Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking.Deletion.pdf Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking.Amplification.genes.pdf Output/Gistic/BWA-CNVkit.ADC-AAH.Smoking.Deletion.genes.pdf Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking.Amplification.genes.pdf Output/Gistic/BWA-CNVkit.ADC-Precancer.Smoking.Deletion.genes.pdf Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking.Amplification.genes.pdf Output/Gistic/BWA-CNVkit.ADC-Primary.Smoking.Deletion.genes.pdf
step43-ADC-Smoking-Bowtie2: Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking_Ex.tsv Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking_Current.tsv Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking_Ex.tsv Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking_Current.tsv Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking_Ex.tsv Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking_Current.tsv Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking_Current/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking_Current/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking_Ex/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking_Current/all_lesions.conf_99.txt Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking.Amplification.pdf Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking.Deletion.pdf Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking.Amplification.pdf Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking.Deletion.pdf Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking.Amplification.pdf Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking.Deletion.pdf Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking.Amplification.genes.pdf Output/Gistic/Bowtie2-CNVkit.ADC-AAH.Smoking.Deletion.genes.pdf Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking.Amplification.genes.pdf Output/Gistic/Bowtie2-CNVkit.ADC-Precancer.Smoking.Deletion.genes.pdf Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking.Amplification.genes.pdf Output/Gistic/Bowtie2-CNVkit.ADC-Primary.Smoking.Deletion.genes.pdf
step43: step43-SQC-stage-BWA step43-SQC-stage-Bowtie2 step43-ADC-stage-BWA step43-ADC-stage-Bowtie2 step43-SQC-Recur-BWA step43-SQC-Recur-Bowtie2 step43-ADC-Recur-BWA step43-ADC-Recur-Bowtie2 step43-SQC-Smoking-BWA step43-SQC-Smoking-Bowtie2 step43-ADC-Smoking-BWA step43-ADC-Smoking-Bowtie2
.PHONY: step43-SQC-stage-BWA step43-SQC-stage-Bowtie2 step43-ADC-stage-BWA step43-ADC-stage-Bowtie2 step43-SQC-Recur-BWA step43-SQC-Recur-Bowtie2 step43-ADC-Recur-BWA step43-ADC-Recur-Bowtie2 step43-SQC-Smoking-BWA step43-SQC-Smoking-Bowtie2 step43-ADC-Smoking-BWA step43-ADC-Smoking-Bowtie2 step43

# Step 44 (Mutation intersection & SigProfiler Signature score)
Output/MutationIntersection_vs_Signature:
	mkdir -p $@

Output/MutationIntersection_vs_Signature/BWA-SQC-SBS.tar: Python/compare_mutation_signature.py Output/MutationSharedProportion/BWA.tsv Output/SigProfiler/BWA/COSMIC_SBS96_Activities_refit.txt | build.log Output/MutationIntersection_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection_vs_Signature/BWA-ADC-SBS.tar: Python/compare_mutation_signature.py Output/MutationSharedProportion/BWA.tsv Output/SigProfiler/BWA/COSMIC_SBS96_Activities_refit.txt | build.log Output/MutationIntersection_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection_vs_Signature/BWA-SQC-DBS.tar: Python/compare_mutation_signature.py Output/MutationSharedProportion/BWA.tsv Output/SigProfiler/BWA/COSMIC_DBS78_Activities_refit.txt  | build.log Output/MutationIntersection_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection_vs_Signature/BWA-ADC-DBS.tar: Python/compare_mutation_signature.py Output/MutationSharedProportion/BWA.tsv Output/SigProfiler/BWA/COSMIC_DBS78_Activities_refit.txt | build.log Output/MutationIntersection_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection_vs_Signature/BWA-SQC-Indel.tar: Python/compare_mutation_signature.py Output/MutationSharedProportion/BWA.tsv Output/SigProfiler/BWA/COSMIC_ID83_Activities_refit.txt | build.log Output/MutationIntersection_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection_vs_Signature/BWA-ADC-Indel.tar: Python/compare_mutation_signature.py Output/MutationSharedProportion/BWA.tsv Output/SigProfiler/BWA/COSMIC_ID83_Activities_refit.txt | build.log Output/MutationIntersection_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection_vs_Signature/Bowtie2-SQC-SBS.tar: Python/compare_mutation_signature.py Output/MutationSharedProportion/Bowtie2.tsv Output/SigProfiler/Bowtie2/COSMIC_SBS96_Activities_refit.txt | build.log Output/MutationIntersection_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection_vs_Signature/Bowtie2-ADC-SBS.tar: Python/compare_mutation_signature.py Output/MutationSharedProportion/Bowtie2.tsv Output/SigProfiler/Bowtie2/COSMIC_SBS96_Activities_refit.txt | build.log Output/MutationIntersection_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection_vs_Signature/Bowtie2-SQC-DBS.tar: Python/compare_mutation_signature.py Output/MutationSharedProportion/Bowtie2.tsv Output/SigProfiler/Bowtie2/COSMIC_DBS78_Activities_refit.txt | build.log Output/MutationIntersection_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection_vs_Signature/Bowtie2-ADC-DBS.tar: Python/compare_mutation_signature.py Output/MutationSharedProportion/Bowtie2.tsv Output/SigProfiler/Bowtie2/COSMIC_DBS78_Activities_refit.txt | build.log Output/MutationIntersection_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection_vs_Signature/Bowtie2-SQC-Indel.tar: Python/compare_mutation_signature.py Output/MutationSharedProportion/Bowtie2.tsv Output/SigProfiler/Bowtie2/COSMIC_ID83_Activities_refit.txt | build.log Output/MutationIntersection_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection_vs_Signature/Bowtie2-ADC-Indel.tar: Python/compare_mutation_signature.py Output/MutationSharedProportion/Bowtie2.tsv Output/SigProfiler/Bowtie2/COSMIC_ID83_Activities_refit.txt | build.log Output/MutationIntersection_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step44-BWA-SQC: Output/MutationIntersection_vs_Signature/BWA-SQC-SBS.tar Output/MutationIntersection_vs_Signature/BWA-SQC-DBS.tar Output/MutationIntersection_vs_Signature/BWA-SQC-Indel.tar
step44-BWA-ADC: Output/MutationIntersection_vs_Signature/BWA-ADC-SBS.tar Output/MutationIntersection_vs_Signature/BWA-ADC-DBS.tar Output/MutationIntersection_vs_Signature/BWA-ADC-Indel.tar
step44-Bowtie2-SQC: Output/MutationIntersection_vs_Signature/Bowtie2-SQC-SBS.tar Output/MutationIntersection_vs_Signature/Bowtie2-SQC-DBS.tar Output/MutationIntersection_vs_Signature/Bowtie2-SQC-Indel.tar
step44-Bowtie2-ADC: Output/MutationIntersection_vs_Signature/Bowtie2-ADC-SBS.tar Output/MutationIntersection_vs_Signature/Bowtie2-ADC-DBS.tar Output/MutationIntersection_vs_Signature/Bowtie2-ADC-Indel.tar
step44: step44-BWA-SQC step44-BWA-ADC step44-Bowtie2-SQC step44-Bowtie2-ADC
.PHONY: step44-BWA-SQC step44-BWA-ADC step44-Bowtie2-SQC step44-Bowtie2-ADC step44

# Step 45 (Clinical Values & Mutational Signature score)
Output/Clinical_vs_Signature:
	mkdir -p $@

Output/Clinical_vs_Signature/BWA-SQC-VDT_SBS.tar: Python/compare_VDT_signature.py Output/SigProfiler/BWA/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-ADC-VDT_SBS.tar: Python/compare_VDT_signature.py Output/SigProfiler/BWA/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-SQC-VDT_DBS.tar: Python/compare_VDT_signature.py Output/SigProfiler/BWA/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-ADC-VDT_DBS.tar: Python/compare_VDT_signature.py Output/SigProfiler/BWA/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-SQC-VDT_Indel.tar: Python/compare_VDT_signature.py Output/SigProfiler/BWA/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-ADC-VDT_Indel.tar: Python/compare_VDT_signature.py Output/SigProfiler/BWA/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-SQC-VDT_SBS.tar: Python/compare_VDT_signature.py Output/SigProfiler/Bowtie2/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-ADC-VDT_SBS.tar: Python/compare_VDT_signature.py Output/SigProfiler/Bowtie2/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-SQC-VDT_DBS.tar: Python/compare_VDT_signature.py Output/SigProfiler/Bowtie2/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-ADC-VDT_DBS.tar: Python/compare_VDT_signature.py Output/SigProfiler/Bowtie2/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-SQC-VDT_Indel.tar: Python/compare_VDT_signature.py Output/SigProfiler/Bowtie2/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-ADC-VDT_Indel.tar: Python/compare_VDT_signature.py Output/SigProfiler/Bowtie2/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-SQC-RFS_SBS.tar: Python/compare_RFS_signature.py Output/SigProfiler/BWA/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-ADC-RFS_SBS.tar: Python/compare_RFS_signature.py Output/SigProfiler/BWA/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-SQC-RFS_DBS.tar: Python/compare_RFS_signature.py Output/SigProfiler/BWA/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-ADC-RFS_DBS.tar: Python/compare_RFS_signature.py Output/SigProfiler/BWA/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-SQC-RFS_Indel.tar: Python/compare_RFS_signature.py Output/SigProfiler/BWA/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-ADC-RFS_Indel.tar: Python/compare_RFS_signature.py Output/SigProfiler/BWA/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-SQC-RFS_SBS.tar: Python/compare_RFS_signature.py Output/SigProfiler/Bowtie2/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-ADC-RFS_SBS.tar: Python/compare_RFS_signature.py Output/SigProfiler/Bowtie2/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-SQC-RFS_DBS.tar: Python/compare_RFS_signature.py Output/SigProfiler/Bowtie2/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-ADC-RFS_DBS.tar: Python/compare_RFS_signature.py Output/SigProfiler/Bowtie2/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-SQC-RFS_Indel.tar: Python/compare_RFS_signature.py Output/SigProfiler/Bowtie2/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-ADC-RFS_Indel.tar: Python/compare_RFS_signature.py Output/SigProfiler/Bowtie2/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-SQC-OS_SBS.tar: Python/compare_OS_signature.py Output/SigProfiler/BWA/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-ADC-OS_SBS.tar: Python/compare_OS_signature.py Output/SigProfiler/BWA/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-SQC-OS_DBS.tar: Python/compare_OS_signature.py Output/SigProfiler/BWA/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-ADC-OS_DBS.tar: Python/compare_OS_signature.py Output/SigProfiler/BWA/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-SQC-OS_Indel.tar: Python/compare_OS_signature.py Output/SigProfiler/BWA/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/BWA-ADC-OS_Indel.tar: Python/compare_OS_signature.py Output/SigProfiler/BWA/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-SQC-OS_SBS.tar: Python/compare_OS_signature.py Output/SigProfiler/Bowtie2/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-ADC-OS_SBS.tar: Python/compare_OS_signature.py Output/SigProfiler/Bowtie2/COSMIC_SBS96_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-SQC-OS_DBS.tar: Python/compare_OS_signature.py Output/SigProfiler/Bowtie2/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-ADC-OS_DBS.tar: Python/compare_OS_signature.py Output/SigProfiler/Bowtie2/COSMIC_DBS78_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-SQC-OS_Indel.tar: Python/compare_OS_signature.py Output/SigProfiler/Bowtie2/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Signature/Bowtie2-ADC-OS_Indel.tar: Python/compare_OS_signature.py Output/SigProfiler/Bowtie2/COSMIC_ID83_Activities_refit.txt $(CLINICAL_DATA_CSV) | build.log Output/Clinical_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

step45-VDT: Output/Clinical_vs_Signature/BWA-SQC-VDT_SBS.tar Output/Clinical_vs_Signature/BWA-ADC-VDT_SBS.tar Output/Clinical_vs_Signature/BWA-SQC-VDT_DBS.tar Output/Clinical_vs_Signature/BWA-ADC-VDT_DBS.tar Output/Clinical_vs_Signature/BWA-SQC-VDT_Indel.tar Output/Clinical_vs_Signature/BWA-ADC-VDT_Indel.tar Output/Clinical_vs_Signature/Bowtie2-SQC-VDT_SBS.tar Output/Clinical_vs_Signature/Bowtie2-ADC-VDT_SBS.tar Output/Clinical_vs_Signature/Bowtie2-SQC-VDT_DBS.tar Output/Clinical_vs_Signature/Bowtie2-ADC-VDT_DBS.tar Output/Clinical_vs_Signature/Bowtie2-SQC-VDT_Indel.tar Output/Clinical_vs_Signature/Bowtie2-ADC-VDT_Indel.tar
step45-RFS: Output/Clinical_vs_Signature/BWA-SQC-RFS_SBS.tar Output/Clinical_vs_Signature/BWA-ADC-RFS_SBS.tar Output/Clinical_vs_Signature/BWA-SQC-RFS_DBS.tar Output/Clinical_vs_Signature/BWA-ADC-RFS_DBS.tar Output/Clinical_vs_Signature/BWA-SQC-RFS_Indel.tar Output/Clinical_vs_Signature/BWA-ADC-RFS_Indel.tar Output/Clinical_vs_Signature/Bowtie2-SQC-RFS_SBS.tar Output/Clinical_vs_Signature/Bowtie2-ADC-RFS_SBS.tar Output/Clinical_vs_Signature/Bowtie2-SQC-RFS_DBS.tar Output/Clinical_vs_Signature/Bowtie2-ADC-RFS_DBS.tar Output/Clinical_vs_Signature/Bowtie2-SQC-RFS_Indel.tar Output/Clinical_vs_Signature/Bowtie2-ADC-RFS_Indel.tar
step45-OS: Output/Clinical_vs_Signature/BWA-SQC-OS_SBS.tar Output/Clinical_vs_Signature/BWA-ADC-OS_SBS.tar Output/Clinical_vs_Signature/BWA-SQC-OS_DBS.tar Output/Clinical_vs_Signature/BWA-ADC-OS_DBS.tar Output/Clinical_vs_Signature/BWA-SQC-OS_Indel.tar Output/Clinical_vs_Signature/BWA-ADC-OS_Indel.tar Output/Clinical_vs_Signature/Bowtie2-SQC-OS_SBS.tar Output/Clinical_vs_Signature/Bowtie2-ADC-OS_SBS.tar Output/Clinical_vs_Signature/Bowtie2-SQC-OS_DBS.tar Output/Clinical_vs_Signature/Bowtie2-ADC-OS_DBS.tar Output/Clinical_vs_Signature/Bowtie2-SQC-OS_Indel.tar Output/Clinical_vs_Signature/Bowtie2-ADC-OS_Indel.tar
step45: step45-VDT step45-RFS step45-OS
.PHONY: step45-VDT step45-RFS step45-OS step45

# Step 46 (Deconvolution & Copy Number Variations)
Output/Deconvolution_vs_CNV:
	mkdir -p $@

Output/Deconvolution_vs_CNV/SQC.MuSiC-Bowtie2-%.Sequenza-BWA.tar: Python/compare_MuSiC_sequenza.py Output/MuSiC/Bowtie2.SQC.%.tsv $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/SQC.MuSiC-Bowtie2-%.Sequenza-Bowtie2.tar: Python/compare_MuSiC_sequenza.py Output/MuSiC/Bowtie2.SQC.%.tsv $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/SQC.MuSiC-STAR-%.Sequenza-BWA.tar: Python/compare_MuSiC_sequenza.py Output/MuSiC/STAR.SQC.%.tsv $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/SQC.MuSiC-STAR-%.Sequenza-Bowtie2.tar: Python/compare_MuSiC_sequenza.py Output/MuSiC/STAR.SQC.%.tsv $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/ADC.MuSiC-Bowtie2-%.Sequenza-BWA.tar: Python/compare_MuSiC_sequenza.py Output/MuSiC/Bowtie2.ADC.%.tsv $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/ADC.MuSiC-Bowtie2-%.Sequenza-Bowtie2.tar: Python/compare_MuSiC_sequenza.py Output/MuSiC/Bowtie2.ADC.%.tsv $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/ADC.MuSiC-STAR-%.Sequenza-BWA.tar: Python/compare_MuSiC_sequenza.py Output/MuSiC/STAR.ADC.%.tsv $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/ADC.MuSiC-STAR-%.Sequenza-Bowtie2.tar: Python/compare_MuSiC_sequenza.py Output/MuSiC/STAR.ADC.%.tsv $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/SQC.BisqueRNA-Bowtie2-%.Sequenza-BWA.tar: Python/compare_BisqueRNA_sequenza.py Output/BisqueRNA/Bowtie2.SQC.%.tsv $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/SQC.BisqueRNA-Bowtie2-%.Sequenza-Bowtie2.tar: Python/compare_BisqueRNA_sequenza.py Output/BisqueRNA/Bowtie2.SQC.%.tsv $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/SQC.BisqueRNA-STAR-%.Sequenza-BWA.tar: Python/compare_BisqueRNA_sequenza.py Output/BisqueRNA/STAR.SQC.%.tsv $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/SQC.BisqueRNA-STAR-%.Sequenza-Bowtie2.tar: Python/compare_BisqueRNA_sequenza.py Output/BisqueRNA/STAR.SQC.%.tsv $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/ADC.BisqueRNA-Bowtie2-%.Sequenza-BWA.tar: Python/compare_BisqueRNA_sequenza.py Output/BisqueRNA/Bowtie2.ADC.%.tsv $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/ADC.BisqueRNA-Bowtie2-%.Sequenza-Bowtie2.tar: Python/compare_BisqueRNA_sequenza.py Output/BisqueRNA/Bowtie2.ADC.%.tsv $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/ADC.BisqueRNA-STAR-%.Sequenza-BWA.tar: Python/compare_BisqueRNA_sequenza.py Output/BisqueRNA/STAR.ADC.%.tsv $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/ADC.BisqueRNA-STAR-%.Sequenza-Bowtie2.tar: Python/compare_BisqueRNA_sequenza.py Output/BisqueRNA/STAR.ADC.%.tsv $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/SQC.SCDC-Bowtie2-%.Sequenza-BWA.tar: Python/compare_MuSiC_sequenza.py Output/MuSiC/Bowtie2.SQC.%.tsv $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/SQC.SCDC-Bowtie2-%.Sequenza-Bowtie2.tar: Python/compare_MuSiC_sequenza.py Output/MuSiC/Bowtie2.SQC.%.tsv $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/SQC.SCDC-STAR-%.Sequenza-BWA.tar: Python/compare_MuSiC_sequenza.py Output/MuSiC/STAR.SQC.%.tsv $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/SQC.SCDC-STAR-%.Sequenza-Bowtie2.tar: Python/compare_MuSiC_sequenza.py Output/MuSiC/STAR.SQC.%.tsv $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/ADC.SCDC-Bowtie2-%.Sequenza-BWA.tar: Python/compare_MuSiC_sequenza.py Output/MuSiC/Bowtie2.ADC.%.tsv $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/ADC.SCDC-Bowtie2-%.Sequenza-Bowtie2.tar: Python/compare_MuSiC_sequenza.py Output/MuSiC/Bowtie2.ADC.%.tsv $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/ADC.SCDC-STAR-%.Sequenza-BWA.tar: Python/compare_MuSiC_sequenza.py Output/MuSiC/STAR.ADC.%.tsv $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/ADC.SCDC-STAR-%.Sequenza-Bowtie2.tar: Python/compare_MuSiC_sequenza.py Output/MuSiC/STAR.ADC.%.tsv $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/SQC.MuSiC-Bowtie2-%.PureCN-BWA.tar: Python/compare_MuSiC_PureCN.py Output/MuSiC/Bowtie2.SQC.%.tsv $(addprefix Output/PureCN/BWA/,$(addsuffix .csv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/SQC.MuSiC-Bowtie2-%.PureCN-Bowtie2.tar: Python/compare_MuSiC_PureCN.py Output/MuSiC/Bowtie2.SQC.%.tsv $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .csv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/SQC.MuSiC-STAR-%.PureCN-BWA.tar: Python/compare_MuSiC_PureCN.py Output/MuSiC/STAR.SQC.%.tsv $(addprefix Output/PureCN/BWA/,$(addsuffix .csv,$(WES_REAL_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/Deconvolution_vs_CNV/SQC.MuSiC-STAR-%.PureCN-Bowtie2.tar: Python/compare_MuSiC_PureCN.py Output/MuSiC/STAR.SQC.%.tsv $(addprefix Output/PureCN/Bowtie2/,$(addsuffix .csv,$(WES_TUMOR_ID))) $(CLINICAL_DATA_CSV) | Output/Deconvolution_vs_CNV build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

step46-MuSiC-Sequenza: Output/Deconvolution_vs_CNV/SQC.MuSiC-Bowtie2-GSE131907.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/SQC.MuSiC-Bowtie2-GSE131907.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/SQC.MuSiC-STAR-GSE131907.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/SQC.MuSiC-STAR-GSE131907.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/SQC.MuSiC-Bowtie2-GSE162498.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/SQC.MuSiC-Bowtie2-GSE162498.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/SQC.MuSiC-STAR-GSE162498.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/SQC.MuSiC-STAR-GSE162498.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/ADC.MuSiC-Bowtie2-GSE131907.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/ADC.MuSiC-Bowtie2-GSE131907.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/ADC.MuSiC-STAR-GSE131907.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/ADC.MuSiC-STAR-GSE131907.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/ADC.MuSiC-Bowtie2-GSE162498.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/ADC.MuSiC-Bowtie2-GSE162498.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/ADC.MuSiC-STAR-GSE162498.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/ADC.MuSiC-STAR-GSE162498.Sequenza-Bowtie2.tar
step46-BisqueRNA-Sequenza: Output/Deconvolution_vs_CNV/SQC.BisqueRNA-Bowtie2-GSE131907.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/SQC.BisqueRNA-Bowtie2-GSE131907.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/SQC.BisqueRNA-STAR-GSE131907.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/SQC.BisqueRNA-STAR-GSE131907.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/SQC.BisqueRNA-Bowtie2-GSE162498.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/SQC.BisqueRNA-Bowtie2-GSE162498.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/SQC.BisqueRNA-STAR-GSE162498.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/SQC.BisqueRNA-STAR-GSE162498.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/ADC.BisqueRNA-Bowtie2-GSE131907.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/ADC.BisqueRNA-Bowtie2-GSE131907.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/ADC.BisqueRNA-STAR-GSE131907.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/ADC.BisqueRNA-STAR-GSE131907.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/ADC.BisqueRNA-Bowtie2-GSE162498.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/ADC.BisqueRNA-Bowtie2-GSE162498.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/ADC.BisqueRNA-STAR-GSE162498.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/ADC.BisqueRNA-STAR-GSE162498.Sequenza-Bowtie2.tar
step46-SCDC-Sequenza: Output/Deconvolution_vs_CNV/SQC.SCDC-Bowtie2-GSE131907.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/SQC.SCDC-Bowtie2-GSE131907.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/SQC.SCDC-STAR-GSE131907.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/SQC.SCDC-STAR-GSE131907.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/SQC.SCDC-Bowtie2-GSE162498.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/SQC.SCDC-Bowtie2-GSE162498.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/SQC.SCDC-STAR-GSE162498.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/SQC.SCDC-STAR-GSE162498.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/ADC.SCDC-Bowtie2-GSE131907.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/ADC.SCDC-Bowtie2-GSE131907.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/ADC.SCDC-STAR-GSE131907.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/ADC.SCDC-STAR-GSE131907.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/ADC.SCDC-Bowtie2-GSE162498.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/ADC.SCDC-Bowtie2-GSE162498.Sequenza-Bowtie2.tar Output/Deconvolution_vs_CNV/ADC.SCDC-STAR-GSE162498.Sequenza-BWA.tar Output/Deconvolution_vs_CNV/ADC.SCDC-STAR-GSE162498.Sequenza-Bowtie2.tar
step46: step46-MuSiC-Sequenza step46-BisqueRNA-Sequenza step46-SCDC-Sequenza
.PHONY: step46-MuSiC-Sequenza step46-BisqueRNA-Sequenza step46-SCDC-Sequenza step46

# Step 47 (Kaplan-Meier Plot with Shared Proportion)
Output/Clinical_vs_KM-Plot:
	mkdir -p $@

Output/Clinical_vs_KM-Plot/BWA-SQC-SharedProportion.tar: Python/draw_KM_SharedProportion.py Output/MutationSharedProportion/BWA.tsv | build.log Output/Clinical_vs_KM-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --median 1> $@.stdout 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_KM-Plot/BWA-SQC-SharedProportion.mean.tar: Python/draw_KM_SharedProportion.py Output/MutationSharedProportion/BWA.tsv | build.log Output/Clinical_vs_KM-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --mean 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_KM-Plot/BWA-ADC-SharedProportion.tar: Python/draw_KM_SharedProportion.py Output/MutationSharedProportion/BWA.tsv | build.log Output/Clinical_vs_KM-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --median 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_KM-Plot/BWA-ADC-SharedProportion.mean.tar: Python/draw_KM_SharedProportion.py Output/MutationSharedProportion/BWA.tsv | build.log Output/Clinical_vs_KM-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --mean 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_KM-Plot/Bowtie2-SQC-SharedProportion.tar: Python/draw_KM_SharedProportion.py Output/MutationSharedProportion/Bowtie2.tsv | build.log Output/Clinical_vs_KM-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --median 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_KM-Plot/Bowtie2-SQC-SharedProportion.mean.tar: Python/draw_KM_SharedProportion.py Output/MutationSharedProportion/Bowtie2.tsv | build.log Output/Clinical_vs_KM-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --mean 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_KM-Plot/Bowtie2-ADC-SharedProportion.tar: Python/draw_KM_SharedProportion.py Output/MutationSharedProportion/Bowtie2.tsv | build.log Output/Clinical_vs_KM-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --median 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_KM-Plot/Bowtie2-ADC-SharedProportion.mean.tar: Python/draw_KM_SharedProportion.py Output/MutationSharedProportion/Bowtie2.tsv | build.log Output/Clinical_vs_KM-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --mean 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_KM-Plot/BWA-SQC-SharedProportion.cutting.tar: Python/draw_KM_SharedProportion.py Output/MutationSharedProportion/BWA.tsv | build.log Output/Clinical_vs_KM-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cutting --median 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_KM-Plot/BWA-SQC-SharedProportion.mean.cutting.tar: Python/draw_KM_SharedProportion.py Output/MutationSharedProportion/BWA.tsv | build.log Output/Clinical_vs_KM-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cutting --mean 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_KM-Plot/BWA-ADC-SharedProportion.cutting.tar: Python/draw_KM_SharedProportion.py Output/MutationSharedProportion/BWA.tsv | build.log Output/Clinical_vs_KM-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@)  --ADC --cutting --median 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_KM-Plot/BWA-ADC-SharedProportion.mean.cutting.tar: Python/draw_KM_SharedProportion.py Output/MutationSharedProportion/BWA.tsv | build.log Output/Clinical_vs_KM-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cutting --mean 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_KM-Plot/Bowtie2-SQC-SharedProportion.cutting.tar: Python/draw_KM_SharedProportion.py Output/MutationSharedProportion/Bowtie2.tsv | build.log Output/Clinical_vs_KM-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cutting --median 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_KM-Plot/Bowtie2-SQC-SharedProportion.mean.cutting.tar: Python/draw_KM_SharedProportion.py Output/MutationSharedProportion/Bowtie2.tsv | build.log Output/Clinical_vs_KM-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cutting --mean 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_KM-Plot/Bowtie2-ADC-SharedProportion.cutting.tar: Python/draw_KM_SharedProportion.py Output/MutationSharedProportion/Bowtie2.tsv | build.log Output/Clinical_vs_KM-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cutting --median 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_KM-Plot/Bowtie2-ADC-SharedProportion.mean.cutting.tar: Python/draw_KM_SharedProportion.py Output/MutationSharedProportion/Bowtie2.tsv | build.log Output/Clinical_vs_KM-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cutting --mean 1> $@.stdout 2> $@.stderr

step47-SQC: Output/Clinical_vs_KM-Plot/BWA-SQC-SharedProportion.tar Output/Clinical_vs_KM-Plot/Bowtie2-SQC-SharedProportion.tar Output/Clinical_vs_KM-Plot/BWA-SQC-SharedProportion.cutting.tar Output/Clinical_vs_KM-Plot/Bowtie2-SQC-SharedProportion.cutting.tar Output/Clinical_vs_KM-Plot/BWA-SQC-SharedProportion.mean.tar Output/Clinical_vs_KM-Plot/Bowtie2-SQC-SharedProportion.mean.tar Output/Clinical_vs_KM-Plot/BWA-SQC-SharedProportion.mean.cutting.tar Output/Clinical_vs_KM-Plot/Bowtie2-SQC-SharedProportion.mean.cutting.tar
step47-ADC: Output/Clinical_vs_KM-Plot/BWA-ADC-SharedProportion.tar Output/Clinical_vs_KM-Plot/Bowtie2-ADC-SharedProportion.tar Output/Clinical_vs_KM-Plot/BWA-ADC-SharedProportion.cutting.tar Output/Clinical_vs_KM-Plot/Bowtie2-ADC-SharedProportion.cutting.tar Output/Clinical_vs_KM-Plot/BWA-ADC-SharedProportion.mean.tar Output/Clinical_vs_KM-Plot/Bowtie2-ADC-SharedProportion.mean.tar Output/Clinical_vs_KM-Plot/BWA-ADC-SharedProportion.mean.cutting.tar Output/Clinical_vs_KM-Plot/Bowtie2-ADC-SharedProportion.mean.cutting.tar
step47: step47-SQC
.PHONY: step47-SQC step47

# Step 48 (TIMER 2.0)
Output/TIMER Output/TIMER-CNV:
	mkdir -p $@

Output/TIMER/%.SQC.csv: Python/change_TSV_CSV.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log Output/TIMER
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/TIMER/%.ADC.csv: Python/change_TSV_CSV.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | build.log Output/TIMER
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/TIMER/%.estimation_matrix.csv: Output/TIMER/%.csv
	@echo "You have to manually execute TIMER on < http://timer.cistrome.org > with $^ CSV file in order to make $@ CSV file!!"

Output/TIMER/%.tar: Python/draw_violin_TIMER.py Output/TIMER/%.estimation_matrix.csv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/TIMER-CNV/%.Sequenza-BWA.tar: Python/compare_TIMER_sequenza.py Output/TIMER/%.estimation_matrix.csv $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) | build.log Output/TIMER-CNV
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/TIMER-CNV/%.Sequenza-Bowtie2.tar: Python/compare_TIMER_sequenza.py Output/TIMER/%.estimation_matrix.csv $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) | build.log Output/TIMER-CNV
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step48-CSV: Output/TIMER/STAR.SQC.csv Output/TIMER/STAR.TPM.SQC.csv Output/TIMER/STAR.FPKM.SQC.csv Output/TIMER/STAR.ADC.csv Output/TIMER/STAR.TPM.ADC.csv Output/TIMER/STAR.FPKM.ADC.csv Output/TIMER/Bowtie2.SQC.csv Output/TIMER/Bowtie2.TPM.SQC.csv Output/TIMER/Bowtie2.FPKM.SQC.csv Output/TIMER/Bowtie2.ADC.csv Output/TIMER/Bowtie2.TPM.ADC.csv Output/TIMER/Bowtie2.FPKM.ADC.csv
step48-violin: Output/TIMER/STAR.SQC.tar Output/TIMER/STAR.TPM.SQC.tar Output/TIMER/STAR.FPKM.SQC.tar Output/TIMER/STAR.ADC.tar Output/TIMER/STAR.TPM.ADC.tar Output/TIMER/STAR.FPKM.ADC.tar Output/TIMER/Bowtie2.SQC.tar Output/TIMER/Bowtie2.TPM.SQC.tar Output/TIMER/Bowtie2.FPKM.SQC.tar Output/TIMER/Bowtie2.ADC.tar Output/TIMER/Bowtie2.TPM.ADC.tar Output/TIMER/Bowtie2.FPKM.ADC.tar
step48-CNV: Output/TIMER-CNV/STAR.SQC.Sequenza-BWA.tar Output/TIMER-CNV/STAR.ADC.Sequenza-BWA.tar Output/TIMER-CNV/Bowtie2.SQC.Sequenza-BWA.tar Output/TIMER-CNV/Bowtie2.ADC.Sequenza-BWA.tar Output/TIMER-CNV/STAR.SQC.Sequenza-Bowtie2.tar Output/TIMER-CNV/STAR.ADC.Sequenza-Bowtie2.tar Output/TIMER-CNV/Bowtie2.SQC.Sequenza-Bowtie2.tar Output/TIMER-CNV/Bowtie2.ADC.Sequenza-Bowtie2.tar
step48: step48-CSV step48-violin
.PHONY: step48-CSV step48-violin step48

# Step 49 (MiXCR)
Output/MiXCR/WES Output/MiXCR/WTS:
	mkdir -p $@

Output/MiXCR/WES/%.clns Output/MiXCR/WES/%.clonotypes.ALL.txt Output/MiXCR/WES/%.clonotypes.TRB.txt Output/MiXCR/WES/%.vdjca Output/MiXCR/WES/%.report: Data/RawData/%_*.fastq.gz | Output/MiXCR/WES build.log
	$(DOCKER) mixcr analyze shotgun --verbose --force-overwrite --species "HomoSapiens" --only-productive --starting-material "dna" --threads $(CPUS) $(addprefix /,$(sort $^) $(@D)/$*) 1> $(@D)/$*.stdout 2> $(@D)/$*.stderr

Output/MiXCR/WTS/%.clns Output/MiXCR/WTS/%.clonotypes.ALL.txt Output/MiXCR/WTS/%.clonotypes.TRB.txt Output/MiXCR/WTS/%.vdjca Output/MiXCR/WTS/%.report: Data/WTS/%_*.fastq.gz | Output/MiXCR/WTS build.log
	$(DOCKER) mixcr analyze shotgun --verbose --force-overwrite --species "HomoSapiens" --only-productive --starting-material "rna" --threads $(CPUS) $(addprefix /,$(sort $^) $(@D)/$*) 1> $(@D)/$*.stdout 2> $(@D)/$*.stderr

Output/MiXCR/WTS.SQC.Density.pdf: Python/draw_MiXCR_density.py $(addprefix Output/MiXCR/WTS/,$(addsuffix .clonotypes.ALL.txt,$(WTS_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/MiXCR/WTS.ADC.Density.pdf: Python/draw_MiXCR_density.py $(addprefix Output/MiXCR/WTS/,$(addsuffix .clonotypes.ALL.txt,$(WTS_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/MiXCR/WTS.SQC.Diversity.tar: Python/draw_MiXCR_diversity.py $(addprefix Output/MiXCR/WTS/,$(addsuffix .clonotypes.ALL.txt,$(WTS_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/MiXCR/WTS.ADC.Diversity.tar: Python/draw_MiXCR_diversity.py $(addprefix Output/MiXCR/WTS/,$(addsuffix .clonotypes.ALL.txt,$(WTS_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/MiXCR/WTS.SQC.Clonality.pdf: Python/draw_MiXCR_clonality.py $(addprefix Output/MiXCR/WTS/,$(addsuffix .clonotypes.ALL.txt,$(WTS_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/MiXCR/WTS.ADC.Clonality.pdf: Python/draw_MiXCR_clonality.py $(addprefix Output/MiXCR/WTS/,$(addsuffix .clonotypes.ALL.txt,$(WTS_ID))) $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

step49-MiXCR: $(addprefix Output/MiXCR/WES/,$(addsuffix .clonotypes.ALL.txt,$(WES_ID))) $(addprefix Output/MiXCR/WTS/,$(addsuffix .clonotypes.ALL.txt,$(WTS_ID)))
step49-Density: Output/MiXCR/WTS.SQC.Density.pdf Output/MiXCR/WTS.ADC.Density.pdf
step49-Diversity: Output/MiXCR/WTS.SQC.Diversity.tar Output/MiXCR/WTS.ADC.Diversity.tar
step49-Clonality: Output/MiXCR/WTS.SQC.Clonality.pdf Output/MiXCR/WTS.ADC.Clonality.pdf
step49: step49-MiXCR step49-Density step49-Diversity step49-Clonality
.PHONY: step49-MiXCR step49-Density step49-Diversity step49-Clonality step49

# Step 50 (CoMut plot with Copy Number Variations)
Output/SNV_CNV/Mutect-Sequenza:
	mkdir -p $@

Output/SNV_CNV/Mutect-Sequenza/BWA-SQC.pdf: Python/aggregate_mutect_sequenza.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) | build.log Output/SNV_CNV/Mutect-Sequenza
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter-out %/sample_segments.txt,$^) $@) --cnv $(addprefix /,$(filter %/sample_segments.txt,$^)) --cpus $(CPUS) --SQC --include 1> $@.stdout 2> $@.stderr

Output/SNV_CNV/Mutect-Sequenza/BWA-ADC.pdf: Python/aggregate_mutect_sequenza.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) | build.log Output/SNV_CNV/Mutect-Sequenza
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter-out %/sample_segments.txt,$^) $@) --cnv $(addprefix /,$(filter %/sample_segments.txt,$^)) --cpus $(CPUS) --ADC --include 1> $@.stdout 2> $@.stderr

Output/SNV_CNV/Mutect-Sequenza/Bowtie2-SQC.pdf: Python/aggregate_mutect_sequenza.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) | build.log Output/SNV_CNV/Mutect-Sequenza
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter-out %/sample_segments.txt,$^) $@) --cnv $(addprefix /,$(filter %/sample_segments.txt,$^)) --cpus $(CPUS) --SQC --include 1> $@.stdout 2> $@.stderr

Output/SNV_CNV/Mutect-Sequenza/Bowtie2-ADC.pdf: Python/aggregate_mutect_sequenza.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) | build.log Output/SNV_CNV/Mutect-Sequenza
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter-out %/sample_segments.txt,$^) $@) --cnv $(addprefix /,$(filter %/sample_segments.txt,$^)) --cpus $(CPUS) --ADC --include 1> $@.stdout 2> $@.stderr

Output/SNV_CNV/Mutect-Sequenza/BWA-SQC.stage.pdf: Python/aggregate_mutect_sequenza.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) | build.log Output/SNV_CNV/Mutect-Sequenza
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter-out %/sample_segments.txt,$^) $@) --cnv $(addprefix /,$(filter %/sample_segments.txt,$^)) --cpus $(CPUS) --SQC --include --stage 1> $@.stdout 2> $@.stderr

Output/SNV_CNV/Mutect-Sequenza/BWA-ADC.stage.pdf: Python/aggregate_mutect_sequenza.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) | build.log Output/SNV_CNV/Mutect-Sequenza
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter-out %/sample_segments.txt,$^) $@) --cnv $(addprefix /,$(filter %/sample_segments.txt,$^)) --cpus $(CPUS) --ADC --include --stage 1> $@.stdout 2> $@.stderr

Output/SNV_CNV/Mutect-Sequenza/Bowtie2-SQC.stage.pdf: Python/aggregate_mutect_sequenza.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) | build.log Output/SNV_CNV/Mutect-Sequenza
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter-out %/sample_segments.txt,$^) $@) --cnv $(addprefix /,$(filter %/sample_segments.txt,$^)) --cpus $(CPUS) --SQC --include --stage 1> $@.stdout 2> $@.stderr

Output/SNV_CNV/Mutect-Sequenza/Bowtie2-ADC.stage.pdf: Python/aggregate_mutect_sequenza.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) | build.log Output/SNV_CNV/Mutect-Sequenza
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter-out %/sample_segments.txt,$^) $@) --cnv $(addprefix /,$(filter %/sample_segments.txt,$^)) --cpus $(CPUS) --ADC --include --stage 1> $@.stdout 2> $@.stderr

Output/SNV_CNV/Mutect-Sequenza/BWA-SQC.gene.pdf: Python/aggregate_mutect_sequenza.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) | build.log Output/SNV_CNV/Mutect-Sequenza
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter-out %/sample_segments.txt,$^) $@) --cnv $(addprefix /,$(filter %/sample_segments.txt,$^)) --cpus $(CPUS) --SQC --include --gene 1> $@.stdout 2> $@.stderr

Output/SNV_CNV/Mutect-Sequenza/BWA-ADC.gene.pdf: Python/aggregate_mutect_sequenza.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/BWA/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) | build.log Output/SNV_CNV/Mutect-Sequenza
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter-out %/sample_segments.txt,$^) $@) --cnv $(addprefix /,$(filter %/sample_segments.txt,$^)) --cpus $(CPUS) --ADC --include --gene 1> $@.stdout 2> $@.stderr

Output/SNV_CNV/Mutect-Sequenza/Bowtie2-SQC.gene.pdf: Python/aggregate_mutect_sequenza.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/SQC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) | build.log Output/SNV_CNV/Mutect-Sequenza
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter-out %/sample_segments.txt,$^) $@) --cnv $(addprefix /,$(filter %/sample_segments.txt,$^)) --cpus $(CPUS) --SQC --include --gene 1> $@.stdout 2> $@.stderr

Output/SNV_CNV/Mutect-Sequenza/Bowtie2-ADC.gene.pdf: Python/aggregate_mutect_sequenza.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) Output/MutEnricher/Bowtie2/ADC/MutEnricher_gene_hotspot_Fisher_enrichments.txt $(GENE_CENCUS_ONE_CSV) $(CLINICAL_DATA_CSV) $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_segments.txt,$(WES_REAL_TUMOR_ID))) | build.log Output/SNV_CNV/Mutect-Sequenza
	$(DOCKER) $(PYTHON) $(addprefix /,$(filter-out %/sample_segments.txt,$^) $@) --cnv $(addprefix /,$(filter %/sample_segments.txt,$^)) --cpus $(CPUS) --ADC --include --gene 1> $@.stdout 2> $@.stderr

step50-Sequenza: Output/SNV_CNV/Mutect-Sequenza/BWA-SQC.pdf Output/SNV_CNV/Mutect-Sequenza/BWA-ADC.pdf Output/SNV_CNV/Mutect-Sequenza/Bowtie2-SQC.pdf Output/SNV_CNV/Mutect-Sequenza/Bowtie2-ADC.pdf Output/SNV_CNV/Mutect-Sequenza/BWA-SQC.stage.pdf Output/SNV_CNV/Mutect-Sequenza/BWA-ADC.stage.pdf Output/SNV_CNV/Mutect-Sequenza/Bowtie2-SQC.stage.pdf Output/SNV_CNV/Mutect-Sequenza/Bowtie2-ADC.stage.pdf Output/SNV_CNV/Mutect-Sequenza/BWA-SQC.gene.pdf Output/SNV_CNV/Mutect-Sequenza/BWA-ADC.gene.pdf Output/SNV_CNV/Mutect-Sequenza/Bowtie2-SQC.gene.pdf Output/SNV_CNV/Mutect-Sequenza/Bowtie2-ADC.gene.pdf
step50: step50-Sequenza
.PHONY: step50-Sequenza step50

# Step 51 (Clinical Values & Gene Expression Values)
Output/Clinical_vs_DEG:
	mkdir -p $@

Output/Clinical_vs_DEG/%.SQC.Recurrence.tar: Python/compare_clinical_DEG_discrete.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | Output/Clinical_vs_DEG
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --column "Recurrence" --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_DEG/%.ADC.Recurrence.tar: Python/compare_clinical_DEG_discrete.py Output/RSEM/DEG/%.tsv $(CLINICAL_DATA_CSV) | Output/Clinical_vs_DEG
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --column "Recurrence" --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step51-Recur: Output/Clinical_vs_DEG/STAR.SQC.Recurrence.tar Output/Clinical_vs_DEG/STAR.TPM.SQC.Recurrence.tar Output/Clinical_vs_DEG/STAR.FPKM.SQC.Recurrence.tar Output/Clinical_vs_DEG/Bowtie2.SQC.Recurrence.tar Output/Clinical_vs_DEG/Bowtie2.TPM.SQC.Recurrence.tar Output/Clinical_vs_DEG/Bowtie2.FPKM.SQC.Recurrence.tar Output/Clinical_vs_DEG/STAR.ADC.Recurrence.tar Output/Clinical_vs_DEG/STAR.TPM.ADC.Recurrence.tar Output/Clinical_vs_DEG/STAR.FPKM.ADC.Recurrence.tar Output/Clinical_vs_DEG/Bowtie2.ADC.Recurrence.tar Output/Clinical_vs_DEG/Bowtie2.TPM.ADC.Recurrence.tar Output/Clinical_vs_DEG/Bowtie2.FPKM.ADC.Recurrence.tar
step51: step51-Recur
.PHONY: step51-Recur step51

# Step 52 (RTCR)
Output/RTCR:
	mkdir -p $@

step52:
.PHONY: step52

# Step 53 (deconstructSigs)
Output/deconstructSigs Output/deconstructSigs/Violin Output/deconstructSigs/Bar Output/deconstructSigs/Distribution:
	mkdir -p $@

Output/deconstructSigs/BWA.tsv: Python/aggregate_mutect_deconstructSigs.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_REAL_TUMOR_ID))) | Output/deconstructSigs build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bowtie2.tsv: Python/aggregate_mutect_deconstructSigs.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_REAL_TUMOR_ID))) | Output/deconstructSigs build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/%.signature.tsv: R/run_deconstructSigs.R Output/deconstructSigs/%.tsv | build.log
	$(DOCKER) $(R) $(addprefix /,$(filter %.R,$^)) --input $(addprefix /,$(filter %.tsv,$^)) --output $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Violin/%.SQC.Stage.Absolute.tar: Python/plot_signature_deconstructSigs.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Violin/%.SQC.Stage.Relative.tar: Python/plot_signature_deconstructSigs_relative.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Violin/%.ADC.Stage.Absolute.tar: Python/plot_signature_deconstructSigs.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Violin/%.ADC.Stage.Relative.tar: Python/plot_signature_deconstructSigs_relative.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Violin/%.SQC.Recurrence.Absolute.tar: Python/plot_signature_deconstructSigs_clinical.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Violin/%.SQC.Recurrence.Relative.tar: Python/plot_signature_deconstructSigs_clinical_relative.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Violin/%.SQC.Smoking.Absolute.tar: Python/plot_signature_deconstructSigs_clinical.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Violin/%.SQC.Smoking.Relative.tar: Python/plot_signature_deconstructSigs_clinical_relative.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Violin/%.ADC.Recurrence.Absolute.tar: Python/plot_signature_deconstructSigs_clinical.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Violin/%.ADC.Recurrence.Relative.tar: Python/plot_signature_deconstructSigs_clinical_relative.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Violin/%.ADC.Smoking.Absolute.tar: Python/plot_signature_deconstructSigs_clinical.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Violin/%.ADC.Smoking.Relative.tar: Python/plot_signature_deconstructSigs_clinical_relative.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Violin build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bar/%.SQC.Absolute.pdf: Python/plot_signature_deconstructSigs_bar.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Bar build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bar/%.SQC.Relative.pdf: Python/plot_signature_deconstructSigs_bar.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Bar build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bar/%.ADC.Absolute.pdf: Python/plot_signature_deconstructSigs_bar.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Bar build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bar/%.ADC.Relative.pdf: Python/plot_signature_deconstructSigs_bar.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Bar build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bar/%.SQC.Subtype.Absolute.pdf: Python/plot_signature_deconstructSigs_bar_subtype.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Bar build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bar/%.SQC.Subtype.Relative.pdf: Python/plot_signature_deconstructSigs_bar_subtype.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Bar build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bar/%.ADC.Subtype.Absolute.pdf: Python/plot_signature_deconstructSigs_bar_subtype.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Bar build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bar/%.ADC.Subtype.Relative.pdf: Python/plot_signature_deconstructSigs_bar_subtype.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Bar build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bar/%.SQC.Recurrence.Absolute.pdf: Python/plot_signature_deconstructSigs_bar_subtype_clinical.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Bar build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bar/%.SQC.Recurrence.Relative.pdf: Python/plot_signature_deconstructSigs_bar_subtype_clinical.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Bar build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bar/%.ADC.Recurrence.Absolute.pdf: Python/plot_signature_deconstructSigs_bar_subtype_clinical.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Bar build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bar/%.ADC.Recurrence.Relative.pdf: Python/plot_signature_deconstructSigs_bar_subtype_clinical.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Bar build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bar/%.SQC.Smoking.Absolute.pdf: Python/plot_signature_deconstructSigs_bar_subtype_clinical.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Bar build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bar/%.SQC.Smoking.Relative.pdf: Python/plot_signature_deconstructSigs_bar_subtype_clinical.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Bar build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bar/%.ADC.Smoking.Absolute.pdf: Python/plot_signature_deconstructSigs_bar_subtype_clinical.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Bar build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute --compare "Smoking-Detail" "Never" "Ex" "Current" 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Bar/%.ADC.Smoking.Relative.pdf: Python/plot_signature_deconstructSigs_bar_subtype_clinical.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Bar build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative --compare "Smoking-Detail" "Never" "Ex" "Current"  1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Distribution/%.SQC.Absolute.pdf: Python/plot_signature_distribution_deconstructSigs.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Distribution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --absolute 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Distribution/%.SQC.Relative.pdf: Python/plot_signature_distribution_deconstructSigs.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Distribution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --relative 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Distribution/%.ADC.Absolute.pdf: Python/plot_signature_distribution_deconstructSigs.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Distribution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --absolute 1> $@.stdout 2> $@.stderr

Output/deconstructSigs/Distribution/%.ADC.Relative.pdf: Python/plot_signature_distribution_deconstructSigs.py Output/deconstructSigs/%.signature.tsv $(CLINICAL_DATA_CSV) | Output/deconstructSigs/Distribution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --relative 1> $@.stdout 2> $@.stderr

step53-Preparation: Output/deconstructSigs/BWA.tsv Output/deconstructSigs/Bowtie2.tsv
step53-deconstructSigs: Output/deconstructSigs/BWA.signature.tsv Output/deconstructSigs/Bowtie2.signature.tsv
step53-Violin-Stage: Output/deconstructSigs/Violin/BWA.SQC.Stage.Absolute.tar Output/deconstructSigs/Violin/Bowtie2.SQC.Stage.Absolute.tar Output/deconstructSigs/Violin/BWA.SQC.Stage.Relative.tar Output/deconstructSigs/Violin/Bowtie2.SQC.Stage.Relative.tar Output/deconstructSigs/Violin/BWA.ADC.Stage.Absolute.tar Output/deconstructSigs/Violin/Bowtie2.ADC.Stage.Absolute.tar Output/deconstructSigs/Violin/BWA.ADC.Stage.Relative.tar Output/deconstructSigs/Violin/Bowtie2.ADC.Stage.Relative.tar
step53-Violin-Recurrence: Output/deconstructSigs/Violin/BWA.SQC.Recurrence.Absolute.tar Output/deconstructSigs/Violin/Bowtie2.SQC.Recurrence.Absolute.tar Output/deconstructSigs/Violin/BWA.SQC.Recurrence.Relative.tar Output/deconstructSigs/Violin/Bowtie2.SQC.Recurrence.Relative.tar Output/deconstructSigs/Violin/BWA.ADC.Recurrence.Absolute.tar Output/deconstructSigs/Violin/Bowtie2.ADC.Recurrence.Absolute.tar Output/deconstructSigs/Violin/BWA.ADC.Recurrence.Relative.tar Output/deconstructSigs/Violin/Bowtie2.ADC.Recurrence.Relative.tar
step53-Violin-Smoking: Output/deconstructSigs/Violin/BWA.SQC.Smoking.Absolute.tar Output/deconstructSigs/Violin/Bowtie2.SQC.Smoking.Absolute.tar Output/deconstructSigs/Violin/BWA.SQC.Smoking.Relative.tar Output/deconstructSigs/Violin/Bowtie2.SQC.Smoking.Relative.tar Output/deconstructSigs/Violin/BWA.ADC.Smoking.Absolute.tar Output/deconstructSigs/Violin/Bowtie2.ADC.Smoking.Absolute.tar Output/deconstructSigs/Violin/BWA.ADC.Smoking.Relative.tar Output/deconstructSigs/Violin/Bowtie2.ADC.Smoking.Relative.tar
step53-Bar-All: Output/deconstructSigs/Bar/BWA.SQC.Absolute.pdf Output/deconstructSigs/Bar/BWA.SQC.Relative.pdf Output/deconstructSigs/Bar/BWA.ADC.Absolute.pdf Output/deconstructSigs/Bar/BWA.ADC.Relative.pdf Output/deconstructSigs/Bar/Bowtie2.SQC.Absolute.pdf Output/deconstructSigs/Bar/Bowtie2.SQC.Relative.pdf Output/deconstructSigs/Bar/Bowtie2.ADC.Absolute.pdf Output/deconstructSigs/Bar/Bowtie2.ADC.Relative.pdf
step53-Bar-Subtype: Output/deconstructSigs/Bar/BWA.SQC.Subtype.Absolute.pdf Output/deconstructSigs/Bar/BWA.SQC.Subtype.Relative.pdf Output/deconstructSigs/Bar/BWA.ADC.Subtype.Absolute.pdf Output/deconstructSigs/Bar/BWA.ADC.Subtype.Relative.pdf Output/deconstructSigs/Bar/Bowtie2.SQC.Subtype.Absolute.pdf Output/deconstructSigs/Bar/Bowtie2.SQC.Subtype.Relative.pdf Output/deconstructSigs/Bar/Bowtie2.ADC.Subtype.Absolute.pdf Output/deconstructSigs/Bar/Bowtie2.ADC.Subtype.Relative.pdf
step53-Bar-Recurrence: Output/deconstructSigs/Bar/BWA.SQC.Recurrence.Absolute.pdf Output/deconstructSigs/Bar/BWA.SQC.Recurrence.Relative.pdf Output/deconstructSigs/Bar/BWA.ADC.Recurrence.Absolute.pdf Output/deconstructSigs/Bar/BWA.ADC.Recurrence.Relative.pdf Output/deconstructSigs/Bar/Bowtie2.SQC.Recurrence.Absolute.pdf Output/deconstructSigs/Bar/Bowtie2.SQC.Recurrence.Relative.pdf Output/deconstructSigs/Bar/Bowtie2.ADC.Recurrence.Absolute.pdf Output/deconstructSigs/Bar/Bowtie2.ADC.Recurrence.Relative.pdf
step53-Bar-Smoking: Output/deconstructSigs/Bar/BWA.SQC.Smoking.Absolute.pdf Output/deconstructSigs/Bar/BWA.SQC.Smoking.Relative.pdf Output/deconstructSigs/Bar/BWA.ADC.Smoking.Absolute.pdf Output/deconstructSigs/Bar/BWA.ADC.Smoking.Relative.pdf Output/deconstructSigs/Bar/Bowtie2.SQC.Smoking.Absolute.pdf Output/deconstructSigs/Bar/Bowtie2.SQC.Smoking.Relative.pdf Output/deconstructSigs/Bar/Bowtie2.ADC.Smoking.Absolute.pdf Output/deconstructSigs/Bar/Bowtie2.ADC.Smoking.Relative.pdf
step53-Violin: step53-Violin-Stage step53-Violin-Recurrence step53-Violin-Smoking
step53-Bar: step53-Bar-All step53-Bar-Subtype step53-Bar-Recurrence step53-Bar-Smoking
step53-Distribution: Output/deconstructSigs/Distribution/BWA.SQC.Absolute.pdf Output/deconstructSigs/Distribution/BWA.SQC.Relative.pdf Output/deconstructSigs/Distribution/BWA.ADC.Absolute.pdf Output/deconstructSigs/Distribution/BWA.ADC.Relative.pdf Output/deconstructSigs/Distribution/Bowtie2.SQC.Absolute.pdf Output/deconstructSigs/Distribution/Bowtie2.SQC.Relative.pdf Output/deconstructSigs/Distribution/Bowtie2.ADC.Absolute.pdf Output/deconstructSigs/Distribution/Bowtie2.ADC.Relative.pdf
step53: step53-Preparation step53-deconstructSigs step53-Violin step53-Bar step53-Distribution
.PHONY: step53-Preparation step53-deconstructSigs step53-Violin-Stage step53-Violin-Recurrence step53-Violin-Smoking step53-Violin step53-Bar-All step53-Bar-Subtype step53-Bar-Recurrence step53-Bar-Smoking step53-Bar step53-Distribution step53

# Step 54 (Mutational Shared Porportion distribution)
Output/Clinical_Distribution:
	mkdir -p $@

Output/Clinical_Distribution/BWA.SQC.SharedProportion.pdf: Python/draw_distribution_SharedProportion.py Output/MutationSharedProportion/BWA.tsv | Output/Clinical_Distribution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_Distribution/BWA.ADC.SharedProportion.pdf: Python/draw_distribution_SharedProportion.py Output/MutationSharedProportion/BWA.tsv | Output/Clinical_Distribution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_Distribution/Bowtie2.SQC.SharedProportion.pdf: Python/draw_distribution_SharedProportion.py Output/MutationSharedProportion/Bowtie2.tsv | Output/Clinical_Distribution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_Distribution/Bowtie2.ADC.SharedProportion.pdf: Python/draw_distribution_SharedProportion.py Output/MutationSharedProportion/Bowtie2.tsv | Output/Clinical_Distribution build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

step54-SQC: Output/Clinical_Distribution/BWA.SQC.SharedProportion.pdf Output/Clinical_Distribution/Bowtie2.SQC.SharedProportion.pdf
step54-ADC: Output/Clinical_Distribution/BWA.ADC.SharedProportion.pdf Output/Clinical_Distribution/Bowtie2.ADC.SharedProportion.pdf
step54: step54-SQC step54-ADC
.PHONY: step54-SQC step54-ADC step43

# Step 55 (Mutation intersection & deconstructSigs Signature score)
Output/MutationIntersection_vs_Signature/BWA-SQC-deconstructSigs.tar: Python/compare_mutation_signature_deconstructSigs.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) Output/deconstructSigs/BWA.signature.tsv $(CLINICAL_DATA_CSV) | build.log Output/MutationIntersection_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection_vs_Signature/BWA-ADC-deconstructSigs.tar: Python/compare_mutation_signature_deconstructSigs.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) Output/deconstructSigs/BWA.signature.tsv $(CLINICAL_DATA_CSV) | build.log Output/MutationIntersection_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection_vs_Signature/Bowtie2-SQC-deconstructSigs.tar: Python/compare_mutation_signature_deconstructSigs.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) Output/deconstructSigs/Bowtie2.signature.tsv $(CLINICAL_DATA_CSV) | build.log Output/MutationIntersection_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection_vs_Signature/Bowtie2-ADC-deconstructSigs.tar: Python/compare_mutation_signature_deconstructSigs.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) Output/deconstructSigs/Bowtie2.signature.tsv $(CLINICAL_DATA_CSV) | build.log Output/MutationIntersection_vs_Signature
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step55: Output/MutationIntersection_vs_Signature/BWA-SQC-deconstructSigs.tar Output/MutationIntersection_vs_Signature/BWA-ADC-deconstructSigs.tar Output/MutationIntersection_vs_Signature/Bowtie2-SQC-deconstructSigs.tar Output/MutationIntersection_vs_Signature/Bowtie2-ADC-deconstructSigs.tar
.PHONY: step55

# Step 56 (Mutation intersection & Gene)
Output/MutationIntersection:
	mkdir -p $@

Output/MutationIntersection/BWA.SQC.Median.pdf: Python/compare_mutation_shared_proportion.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID)))  Output/MutationSharedProportion/BWA.tsv $(GENE_CENCUS_ONE_CSV) | build.log Output/MutationIntersection
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --median --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection/BWA.SQC.Mean.pdf: Python/compare_mutation_shared_proportion.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID)))  Output/MutationSharedProportion/BWA.tsv $(GENE_CENCUS_ONE_CSV) | build.log Output/MutationIntersection
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --mean --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection/BWA.ADC.Median.pdf: Python/compare_mutation_shared_proportion.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID)))  Output/MutationSharedProportion/BWA.tsv $(GENE_CENCUS_ONE_CSV) | build.log Output/MutationIntersection
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --median --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection/BWA.ADC.Mean.pdf: Python/compare_mutation_shared_proportion.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID)))  Output/MutationSharedProportion/BWA.tsv $(GENE_CENCUS_ONE_CSV) | build.log Output/MutationIntersection
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --mean --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection/Bowtie2.SQC.Median.pdf: Python/compare_mutation_shared_proportion.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/Bowtie2.tsv  $(GENE_CENCUS_ONE_CSV) | build.log Output/MutationIntersection
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --median --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection/Bowtie2.SQC.Mean.pdf: Python/compare_mutation_shared_proportion.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/Bowtie2.tsv  $(GENE_CENCUS_ONE_CSV) | build.log Output/MutationIntersection
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --mean --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection/Bowtie2.ADC.Median.pdf: Python/compare_mutation_shared_proportion.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/Bowtie2.tsv  $(GENE_CENCUS_ONE_CSV) | build.log Output/MutationIntersection
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --median --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/MutationIntersection/Bowtie2.ADC.Mean.pdf: Python/compare_mutation_shared_proportion.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/Bowtie2.tsv $(GENE_CENCUS_ONE_CSV) | build.log Output/MutationIntersection
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --mean --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step56-BWA: Output/MutationIntersection/BWA.SQC.Median.pdf Output/MutationIntersection/BWA.SQC.Mean.pdf Output/MutationIntersection/BWA.ADC.Median.pdf Output/MutationIntersection/BWA.ADC.Mean.pdf
step56-Bowtie2: Output/MutationIntersection/Bowtie2.SQC.Median.pdf Output/MutationIntersection/Bowtie2.SQC.Mean.pdf Output/MutationIntersection/Bowtie2.ADC.Median.pdf Output/MutationIntersection/Bowtie2.ADC.Mean.pdf
step56: step56-BWA step56-Bowtie2
.PHONY: step56-BWA step56-Bowtie2 step56

# Step 57 (PathSeq - WES)
Output/Pathseq/Reference Output/Pathseq/WES/BWA Output/Pathseq/WES/Bowtie2:
	mkdir -p $@

Output/Pathseq/Reference/pathseq_catalog.gz: | Output/Pathseq/Reference
	wget "https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/RefSeq-release212.catalog.gz" -O $@

Output/Pathseq/Reference/pathseq_taxdump.tar.gz: | Output/Pathseq/Reference
	wget "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" -O $@

Output/Pathseq/Reference/pathseq_host.tar.gz: | Output/Pathseq/Reference
	wget "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/pathseq/pathseq_host.tar.gz" -O $@

Output/Pathseq/Reference/pathseq_microbe.tar.gz: | Output/Pathseq/Reference
	wget "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/pathseq/pathseq_microbe.tar.gz" -O $@

Output/Pathseq/Reference/pathseq_host.fa Output/Pathseq/Reference/pathseq_host.bfi Output/Pathseq/Reference/pathseq_host.fa.img: Output/Pathseq/Reference/pathseq_host.tar.gz
	cd $(dir $^) && tar -zxvf $(notdir $^ $@)
	touch $@

Output/Pathseq/Reference/pathseq_microbe.fa Output/Pathseq/Reference/pathseq_microbe.fa.img Output/Pathseq/Reference/pathseq_microbe.dict: Output/Pathseq/Reference/pathseq_microbe.tar.gz
	cd $(dir $^) && tar -zxvf $(notdir $^ $@)
	touch $@

Output/Pathseq/Reference/pathseq_taxonomy.db: Output/Pathseq/Reference/pathseq_taxonomy.tar.gz
	cd $(dir $^) && tar -zxvf $(notdir $^ $@)
	touch $@

Output/Pathseq/WES/BWA/%.tsv Output/Pathseq/WES/BWA/%.bam: Output/WES/%.BWA.Sort.MarkDuplicate.BQSR.bam Output/Pathseq/Reference/pathseq_microbe.fa.img Output/Pathseq/Reference/pathseq_microbe.dict Output/Pathseq/Reference/pathseq_taxonomy.db $(REFERENCE_FASTA) Output/Pathseq/Reference/pathseq_host.bfi Output/Pathseq/Reference/pathseq_host.fa.img | Output/Pathseq/WES/BWA
	$(DOCKER) $(GATK) PathSeqPipelineSpark --input $(addprefix /,$(word 1,$^)) --microbe-bwa-image $(addprefix /,$(word 2,$^)) --microbe-dict $(addprefix /,$(word 3,$^)) --taxonomy-file $(addprefix /,$(word 4,$^)) --reference $(addprefix /,$(word 5,$^)) --kmer-file $(addprefix /,$(word 6,$^)) --filter-bwa-image $(addprefix /,$(word 7,$^)) --scores-output $(addprefix /,$(@D)/$*.tsv) --output $(addprefix /,$(@D)/$*.bam) --is-host-aligned true 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2/%.tsv Output/Pathseq/WES/Bowtie2/%.bam: Output/WES/%.Bowtie2.Sort.MarkDuplicate.BQSR.bam Output/Pathseq/Reference/pathseq_microbe.fa.img Output/Pathseq/Reference/pathseq_microbe.dict Output/Pathseq/Reference/pathseq_taxonomy.db $(REFERENCE_FASTA) Output/Pathseq/Reference/pathseq_host.bfi Output/Pathseq/Reference/pathseq_host.fa.img | Output/Pathseq/WES/Bowtie2
	$(DOCKER) $(GATK) PathSeqPipelineSpark --input $(addprefix /,$(word 1,$^)) --microbe-bwa-image $(addprefix /,$(word 2,$^)) --microbe-dict $(addprefix /,$(word 3,$^)) --taxonomy-file $(addprefix /,$(word 4,$^)) --reference $(addprefix /,$(word 5,$^)) --kmer-file $(addprefix /,$(word 6,$^)) --filter-bwa-image $(addprefix /,$(word 7,$^)) --scores-output $(addprefix /,$(@D)/$*.tsv) --output $(addprefix /,$(@D)/$*.bam) --is-host-aligned true 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA.%.tsv: Python/aggregate_PathSeq_tsv.py $(addprefix Output/Pathseq/WES/BWA/,$(addsuffix .tsv,$(WES_REAL_ID))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2.%.tsv: Python/aggregate_PathSeq_tsv.py $(addprefix Output/Pathseq/WES/Bowtie2/,$(addsuffix .tsv,$(WES_REAL_ID))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-BWA-SQC.%.pdf: Python/aggregate_PathSeq.py Output/Pathseq/WES/BWA.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-BWA-ADC.%.pdf: Python/aggregate_PathSeq.py Output/Pathseq/WES/BWA.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-Bowtie2-SQC.%.pdf: Python/aggregate_PathSeq.py Output/Pathseq/WES/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-Bowtie2-ADC.%.pdf: Python/aggregate_PathSeq.py Output/Pathseq/WES/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-BWA-SQC.%.tar: Python/aggregate_PathSeq_violin.py Output/Pathseq/WES/BWA.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-BWA-ADC.%.tar: Python/aggregate_PathSeq_violin.py Output/Pathseq/WES/BWA.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-Bowtie2-SQC.%.tar: Python/aggregate_PathSeq_violin.py Output/Pathseq/WES/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-Bowtie2-ADC.%.tar: Python/aggregate_PathSeq_violin.py Output/Pathseq/WES/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-BWA-SQC.Smoking.%.tar: Python/aggregate_PathSeq_violin_clinical.py Output/Pathseq/WES/BWA.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-BWA-SQC.Recurrence.%.tar: Python/aggregate_PathSeq_violin_clinical.py Output/Pathseq/WES/BWA.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --level "$*" --compare "Recurrence" "NO" "YES" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-BWA-ADC.Smoking.%.tar: Python/aggregate_PathSeq_violin_clinical.py Output/Pathseq/WES/BWA.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-BWA-ADC.Recurrence.%.tar: Python/aggregate_PathSeq_violin_clinical.py Output/Pathseq/WES/BWA.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --level "$*" --compare "Recurrence" "NO" "YES" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-Bowtie2-SQC.Smoking.%.tar: Python/aggregate_PathSeq_violin_clinical.py Output/Pathseq/WES/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-Bowtie2-SQC.Recurrence.%.tar: Python/aggregate_PathSeq_violin_clinical.py Output/Pathseq/WES/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --level "$*" --compare "Recurrence" "NO" "YES" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-Bowtie2-ADC.Smoking.%.tar: Python/aggregate_PathSeq_violin_clinical.py Output/Pathseq/WES/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-Bowtie2-ADC.Recurrence.%.tar: Python/aggregate_PathSeq_violin_clinical.py Output/Pathseq/WES/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --level "$*" --compare "Recurrence" "NO" "YES" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-BWA-SQC.MutationSharedProportion-Median.%.tar: Python/aggregate_PathSeq_violin_MutationSharedProportion.py Output/Pathseq/WES/BWA.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --median --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-BWA-SQC.MutationSharedProportion-Mean.%.tar: Python/aggregate_PathSeq_violin_MutationSharedProportion.py Output/Pathseq/WES/BWA.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --mean --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-BWA-ADC.MutationSharedProportion-Median.%.tar: Python/aggregate_PathSeq_violin_MutationSharedProportion.py Output/Pathseq/WES/BWA.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --median --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-BWA-ADC.MutationSharedProportion-Mean.%.tar: Python/aggregate_PathSeq_violin_MutationSharedProportion.py Output/Pathseq/WES/BWA.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --mean --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-Bowtie2-SQC.MutationSharedProportion-Median.%.tar: Python/aggregate_PathSeq_violin_MutationSharedProportion.py Output/Pathseq/WES/Bowtie2.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --median --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-Bowtie2-SQC.MutationSharedProportion-Mean.%.tar: Python/aggregate_PathSeq_violin_MutationSharedProportion.py Output/Pathseq/WES/Bowtie2.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --mean --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-Bowtie2-ADC.MutationSharedProportion-Median.%.tar: Python/aggregate_PathSeq_violin_MutationSharedProportion.py Output/Pathseq/WES/Bowtie2.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --median --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Violin-Bowtie2-ADC.MutationSharedProportion-Mean.%.tar: Python/aggregate_PathSeq_violin_MutationSharedProportion.py Output/Pathseq/WES/Bowtie2.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --mean --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Cluster-BWA-SQC.%.pdf: Python/aggregate_PathSeq_clustermap.py Output/Pathseq/WES/BWA.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Cluster-BWA-ADC.%.pdf: Python/aggregate_PathSeq_clustermap.py Output/Pathseq/WES/BWA.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Cluster-Bowtie2-SQC.%.pdf: Python/aggregate_PathSeq_clustermap.py Output/Pathseq/WES/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Cluster-Bowtie2-ADC.%.pdf: Python/aggregate_PathSeq_clustermap.py Output/Pathseq/WES/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-BWA-SQC.Smoking.%.pdf: Python/aggregate_PathSeq_clinical.py Output/Pathseq/WES/BWA.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-BWA-ADC.Smoking.%.pdf: Python/aggregate_PathSeq_clinical.py Output/Pathseq/WES/BWA.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-Bowtie2-SQC.Smoking.%.pdf: Python/aggregate_PathSeq_clinical.py Output/Pathseq/WES/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-Bowtie2-ADC.Smoking.%.pdf: Python/aggregate_PathSeq_clinical.py Output/Pathseq/WES/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-BWA-SQC.Recurrence.%.pdf: Python/aggregate_PathSeq_clinical.py Output/Pathseq/WES/BWA.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-BWA-ADC.Recurrence.%.pdf: Python/aggregate_PathSeq_clinical.py Output/Pathseq/WES/BWA.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-Bowtie2-SQC.Recurrence.%.pdf: Python/aggregate_PathSeq_clinical.py Output/Pathseq/WES/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-Bowtie2-ADC.Recurrence.%.pdf: Python/aggregate_PathSeq_clinical.py Output/Pathseq/WES/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-BWA-SQC.MutationSharedProportion-Median.%.pdf: Python/aggregate_PathSeq_MutationSharedProportion.py Output/Pathseq/WES/BWA.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --median --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-BWA-SQC.MutationSharedProportion-Mean.%.pdf: Python/aggregate_PathSeq_MutationSharedProportion.py Output/Pathseq/WES/BWA.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --mean --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-BWA-ADC.MutationSharedProportion-Median.%.pdf: Python/aggregate_PathSeq_MutationSharedProportion.py Output/Pathseq/WES/BWA.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --median --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-BWA-ADC.MutationSharedProportion-Mean.%.pdf: Python/aggregate_PathSeq_MutationSharedProportion.py Output/Pathseq/WES/BWA.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --mean --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-Bowtie2-SQC.MutationSharedProportion-Median.%.pdf: Python/aggregate_PathSeq_MutationSharedProportion.py Output/Pathseq/WES/Bowtie2.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --median --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-Bowtie2-SQC.MutationSharedProportion-Mean.%.pdf: Python/aggregate_PathSeq_MutationSharedProportion.py Output/Pathseq/WES/Bowtie2.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --mean --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-Bowtie2-ADC.MutationSharedProportion-Median.%.pdf: Python/aggregate_PathSeq_MutationSharedProportion.py Output/Pathseq/WES/Bowtie2.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --median --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bar-Bowtie2-ADC.MutationSharedProportion-Mean.%.pdf: Python/aggregate_PathSeq_MutationSharedProportion.py Output/Pathseq/WES/Bowtie2.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --mean --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA-alpha.%.tsv: Python/calculate_alpha_PathSeq.py Output/Pathseq/WES/BWA.%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2-alpha.%.tsv: Python/calculate_alpha_PathSeq.py Output/Pathseq/WES/Bowtie2.%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA.tree.nwk: Python/make_tree_PathSeq.py $(addprefix Output/Pathseq/WES/BWA/,$(addsuffix .tsv,$(WES_REAL_ID))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2.tree.nwk: Python/make_tree_PathSeq.py $(addprefix Output/Pathseq/WES/Bowtie2/,$(addsuffix .tsv,$(WES_REAL_ID))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA-beta.%.Unweighted.tsv: Python/calculate_beta_PathSeq.py Output/Pathseq/WES/BWA.%.tsv Output/Pathseq/WES/BWA.tree.nwk | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --beta "unweighted_unifrac" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA-beta.%.Weighted.tsv: Python/calculate_beta_PathSeq.py Output/Pathseq/WES/BWA.%.tsv Output/Pathseq/WES/BWA.tree.nwk | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --beta "weighted_unifrac" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2-beta.%.Unweighted.tsv: Python/calculate_beta_PathSeq.py Output/Pathseq/WES/Bowtie2.%.tsv Output/Pathseq/WES/Bowtie2.tree.nwk | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --beta "unweighted_unifrac" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2-beta.%.Weighted.tsv: Python/calculate_beta_PathSeq.py Output/Pathseq/WES/Bowtie2.%.tsv Output/Pathseq/WES/Bowtie2.tree.nwk | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --beta "weighted_unifrac" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA-SQC-alpha.%.tar: Python/draw_alpha_PathSeq.py Output/Pathseq/WES/BWA-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2-SQC-alpha.%.tar: Python/draw_alpha_PathSeq.py Output/Pathseq/WES/Bowtie2-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA-SQC-alpha.Smoking.%.tar: Python/draw_alpha_PathSeq_clinical.py Output/Pathseq/WES/BWA-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA-SQC-alpha.Recurrence.%.tar: Python/draw_alpha_PathSeq_clinical.py Output/Pathseq/WES/BWA-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA-ADC-alpha.Smoking.%.tar: Python/draw_alpha_PathSeq_clinical.py Output/Pathseq/WES/BWA-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA-ADC-alpha.Recurrence.%.tar: Python/draw_alpha_PathSeq_clinical.py Output/Pathseq/WES/BWA-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA-SQC-alpha.MutationSharedProportion-Median.%.tar: Python/draw_alpha_PathSeq_MutationSharedProportion.py Output/Pathseq/WES/BWA-alpha.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --median 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA-SQC-alpha.MutationSharedProportion-Mean.%.tar: Python/draw_alpha_PathSeq_MutationSharedProportion.py Output/Pathseq/WES/BWA-alpha.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --mean 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA-ADC-alpha.MutationSharedProportion-Median.%.tar: Python/draw_alpha_PathSeq_MutationSharedProportion.py Output/Pathseq/WES/BWA-alpha.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --median 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA-ADC-alpha.MutationSharedProportion-Mean.%.tar: Python/draw_alpha_PathSeq_MutationSharedProportion.py Output/Pathseq/WES/BWA-alpha.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --mean 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2-SQC-alpha.MutationSharedProportion-Median.%.tar: Python/draw_alpha_PathSeq_MutationSharedProportion.py Output/Pathseq/WES/Bowtie2-alpha.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --median 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2-SQC-alpha.MutationSharedProportion-Mean.%.tar: Python/draw_alpha_PathSeq_MutationSharedProportion.py Output/Pathseq/WES/Bowtie2-alpha.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --mean 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2-ADC-alpha.MutationSharedProportion-Median.%.tar: Python/draw_alpha_PathSeq_MutationSharedProportion.py Output/Pathseq/WES/Bowtie2-alpha.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --median 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2-ADC-alpha.MutationSharedProportion-Mean.%.tar: Python/draw_alpha_PathSeq_MutationSharedProportion.py Output/Pathseq/WES/Bowtie2-alpha.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --mean 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2-SQC-alpha.Smoking.%.tar: Python/draw_alpha_PathSeq_clinical.py Output/Pathseq/WES/Bowtie2-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2-SQC-alpha.Recurrence.%.tar: Python/draw_alpha_PathSeq_clinical.py Output/Pathseq/WES/Bowtie2-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2-ADC-alpha.Smoking.%.tar: Python/draw_alpha_PathSeq_clinical.py Output/Pathseq/WES/Bowtie2-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2-ADC-alpha.Recurrence.%.tar: Python/draw_alpha_PathSeq_clinical.py Output/Pathseq/WES/Bowtie2-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA-SQC-beta.%.Unweighted.pdf: Python/draw_beta_PathSeq.py Output/Pathseq/WES/BWA-beta.%.Unweighted.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA-SQC-beta.%.Weighted.pdf: Python/draw_beta_PathSeq.py Output/Pathseq/WES/BWA-beta.%.Weighted.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA-ADC-beta.%.Unweighted.pdf: Python/draw_beta_PathSeq.py Output/Pathseq/WES/BWA-beta.%.Unweighted.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/BWA-ADC-beta.%.Weighted.pdf: Python/draw_beta_PathSeq.py Output/Pathseq/WES/BWA-beta.%.Weighted.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2-SQC-beta.%.Unweighted.pdf: Python/draw_beta_PathSeq.py Output/Pathseq/WES/Bowtie2-beta.%.Unweighted.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2-SQC-beta.%.Weighted.pdf: Python/draw_beta_PathSeq.py Output/Pathseq/WES/Bowtie2-beta.%.Weighted.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2-ADC-beta.%.Unweighted.pdf: Python/draw_beta_PathSeq.py Output/Pathseq/WES/Bowtie2-beta.%.Unweighted.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WES/Bowtie2-ADC-beta.%.Weighted.pdf: Python/draw_beta_PathSeq.py Output/Pathseq/WES/Bowtie2-beta.%.Weighted.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step57-download: Output/Pathseq/Reference/pathseq_catalog.gz Output/Pathseq/Reference/pathseq_taxdump.tar.gz Output/Pathseq/Reference/pathseq_host.tar.gz Output/Pathseq/Reference/pathseq_microbe.tar.gz
step57-reference: Output/Pathseq/Reference/pathseq_microbe.fa.img Output/Pathseq/Reference/pathseq_microbe.dict Output/Pathseq/Reference/pathseq_taxonomy.db Output/Pathseq/Reference/pathseq_host.bfi Output/Pathseq/Reference/pathseq_host.fa.img
step57-BWA-Merge: $(addprefix Output/Pathseq/WES/BWA/,$(addsuffix .tsv,$(WES_ID))) $(addprefix Output/Pathseq/WES/BWA.,$(addsuffix .tsv,$(TAXONOMY_TYPES)))
step57-BWA-Bar: $(addprefix Output/Pathseq/WES/Bar-BWA-SQC.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-BWA-ADC.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-BWA-SQC.Smoking.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-BWA-ADC.Smoking.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-BWA-SQC.Recurrence.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-BWA-ADC.Recurrence.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-BWA-SQC.MutationSharedProportion-Median.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-BWA-SQC.MutationSharedProportion-Mean.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-BWA-ADC.MutationSharedProportion-Median.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-BWA-ADC.MutationSharedProportion-Mean.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-SQC.MutationSharedProportion-Median.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-SQC.MutationSharedProportion-Mean.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-ADC.MutationSharedProportion-Median.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-ADC.MutationSharedProportion-Mean.,$(addsuffix .pdf,$(TAXONOMY_TYPES)))
step57-BWA-Violin: $(addprefix Output/Pathseq/WES/Violin-BWA-SQC.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-BWA-ADC.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-BWA-SQC.Smoking.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-BWA-ADC.Smoking.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-BWA-SQC.Recurrence.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-BWA-ADC.Recurrence.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-BWA-SQC.MutationSharedProportion-Median.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-BWA-SQC.MutationSharedProportion-Mean.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-BWA-ADC.MutationSharedProportion-Median.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-BWA-ADC.MutationSharedProportion-Mean.,$(addsuffix .tar,$(TAXONOMY_TYPES)))
step57-BWA-Cluster: $(addprefix Output/Pathseq/WES/Cluster-BWA-SQC.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Cluster-BWA-ADC.,$(addsuffix .pdf,$(TAXONOMY_TYPES)))
step57-BWA-Diversity: $(addprefix Output/Pathseq/WES/BWA-alpha.,$(addsuffix .tsv,$(TAXONOMY_TYPES))) Output/Pathseq/WES/BWA.tree.nwk $(addprefix Output/Pathseq/WES/BWA-beta.,$(addsuffix .Unweighted.tsv,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/BWA-beta.,$(addsuffix .Weighted.tsv,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/BWA-SQC-alpha.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/BWA-SQC-alpha.Smoking.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/BWA-SQC-alpha.Recurrence.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/BWA-ADC-alpha.Smoking.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/BWA-ADC-alpha.Recurrence.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/BWA-SQC-beta.,$(addsuffix .Unweighted.pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/BWA-SQC-beta.,$(addsuffix .Weighted.pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/BWA-ADC-beta.,$(addsuffix .Unweighted.pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/BWA-ADC-beta.,$(addsuffix .Weighted.pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/BWA-SQC-alpha.MutationSharedProportion-Median.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/BWA-SQC-alpha.MutationSharedProportion-Mean.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/BWA-ADC-alpha.MutationSharedProportion-Median.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/BWA-ADC-alpha.MutationSharedProportion-Mean.,$(addsuffix .tar,$(TAXONOMY_TYPES)))
step57-BWA: step57-BWA-Merge step57-BWA-Bar step57-BWA-Violin step57-BWA-cluster step57-BWA-Diversity
step57-Bowtie2-Merge: $(addprefix Output/Pathseq/WES/Bowtie2/,$(addsuffix .tsv,$(WES_ID))) $(addprefix Output/Pathseq/WES/Bowtie2.,$(addsuffix .tsv,$(TAXONOMY_TYPES)))
step57-Bowtie2-Bar: $(addprefix Output/Pathseq/WES/Bar-Bowtie2-SQC.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-ADC.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-SQC.Smoking.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-ADC.Smoking.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-SQC.Recurrence.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-ADC.Recurrence.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-SQC.MutationSharedProportion-Median.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-SQC.MutationSharedProportion-Mean.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-ADC.MutationSharedProportion-Median.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-ADC.MutationSharedProportion-Mean.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-SQC.MutationSharedProportion-Median.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-SQC.MutationSharedProportion-Mean.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-ADC.MutationSharedProportion-Median.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bar-Bowtie2-ADC.MutationSharedProportion-Mean.,$(addsuffix .pdf,$(TAXONOMY_TYPES)))
step57-Bowtie2-Violin: $(addprefix Output/Pathseq/WES/Violin-Bowtie2-SQC.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-Bowtie2-ADC.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-Bowtie2-SQC.Smoking.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-Bowtie2-ADC.Smoking.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-Bowtie2-SQC.Recurrence.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-Bowtie2-ADC.Recurrence.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-Bowtie2-SQC.MutationSharedProportion-Median.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-Bowtie2-SQC.MutationSharedProportion-Mean.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-Bowtie2-ADC.MutationSharedProportion-Median.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Violin-Bowtie2-ADC.MutationSharedProportion-Mean.,$(addsuffix .tar,$(TAXONOMY_TYPES)))
step57-Bowtie2-Cluster: $(addprefix Output/Pathseq/WES/Cluster-Bowtie2-SQC.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Cluster-Bowtie2-ADC.,$(addsuffix .pdf,$(TAXONOMY_TYPES)))
step57-Bowtie2-Diversity:  $(addprefix Output/Pathseq/WES/Bowtie2-alpha.,$(addsuffix .tsv,$(TAXONOMY_TYPES))) Output/Pathseq/WES/Bowtie2.tree.nwk $(addprefix Output/Pathseq/WES/Bowtie2-beta.,$(addsuffix .Unweighted.tsv,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bowtie2-beta.,$(addsuffix .Weighted.tsv,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bowtie2-SQC-alpha.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bowtie2-SQC-alpha.Smoking.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bowtie2-SQC-alpha.Recurrence.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bowtie2-ADC-alpha.Smoking.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bowtie2-ADC-alpha.Recurrence.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bowtie2-SQC-beta.,$(addsuffix .Unweighted.pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bowtie2-SQC-beta.,$(addsuffix .Weighted.pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bowtie2-ADC-beta.,$(addsuffix .Unweighted.pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bowtie2-ADC-beta.,$(addsuffix .Weighted.pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bowtie2-SQC-alpha.MutationSharedProportion-Median.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bowtie2-SQC-alpha.MutationSharedProportion-Mean.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bowtie2-ADC-alpha.MutationSharedProportion-Median.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WES/Bowtie2-ADC-alpha.MutationSharedProportion-Mean.,$(addsuffix .tar,$(TAXONOMY_TYPES)))
step57-Bowtie2: step57-Bowtie2-Merge step57-Bowtie2-Bar step57-Bowtie2-Violin step57-Bowtie2-Cluster step57-Bowtie2-Diversity
step57: step57-download step57-reference step57-BWA-Merge step57-BWA-Bar step57-BWA-Violin step57-BWA-Cluster step57-BWA-Diversity step57-BWA step57-Bowtie2-Merge step57-Bowtie2-Bar step57-Bowtie2-Violin step57-Bowtie2-Cluster step57-Bowtie2-Diversity step57-Bowtie2
.PHONY: step57-download step57-reference step57-BWA-Merge step57-BWA-Bar step57-BWA-Violin step57-BWA-Cluster step57-BWA-Diversity step57-BWA step57-Bowtie2-Merge step57-Bowtie2-Bar step57-Bowtie2-Violin step57-Bowtie2-Cluster step57-Bowtie2-Diversity step57-Bowtie2

# Step 58 (WTS STAR alignment with unmapped)
Output/WTS/STAR/Reference:
	mkdir -p $@

Output/WTS/STAR/Reference/Genome: $(REFERENCE_FASTA) | Output/WTS/STAR/Reference build.log
	$(DOCKER) STAR --runMode genomeGenerate --genomeDir $(addprefix /,$(@D)) --genomeFastaFiles $(addprefix /,$^) --runThreadN $(CPUS) 1> $@.stdout 2> $@.stderr

Output/WTS/STAR/%/Aligned.out.sam: Data/WTS/%_*.fastq.gz Output/WTS/STAR/Reference/Genome | build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D)) && mkdir $(@D)
	$(WORKDIR_DOCKER) STAR --genomeDir $(addprefix /,$(dir $(filter %/Genome,$^))) --readFilesIn $(addprefix /,$(sort $(filter %.gz,$^))) --readFilesCommand zcat --outSAMunmapped "Within" "KeepPairs" --runThreadN $(CPUS) 1> $@.stdout 2> $@.stderr

Output/WTS/STAR/%.bam: $(REFERENCE_FASTA) Output/WTS/STAR/%/Aligned.out.sam | build.log
	$(DOCKER) samtools view -b -h --threads $(CPUS) --reference $(addprefix /,$(filter %.fasta,$^)) $(addprefix /,$(filter %.sam,$^)) -o $(addprefix /,$@)

Output/WTS/STAR/%.Sort.bam: $(REFERENCE_FASTA) Output/WTS/STAR/%.bam | build.log
	$(DOCKER) $(PICARD) AddOrReplaceReadGroups --INPUT $(addprefix /,$(filter %.bam,$^)) --OUTPUT $(addprefix /,$@) --RGLB "Library" --RGPL "Illumina" --RGPU "PlatformUnit" --RGSM "$*" --RGID "$*" --COMPRESSION_LEVEL 9 --REFERENCE_SEQUENCE $(addprefix /,$(filter %.fasta,$^)) --SORT_ORDER "coordinate" 1> $@.stdout 2> $@.stderr

Output/WTS/STAR/%.Sort.bai: Output/WTS/STAR/%.Sort.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/WTS/STAR/%.Sort.dm-tmp.bam: $(REFERENCE_FASTA) Output/WTS/STAR/%.Sort.bam Output/WTS/STAR/%.Sort.bai | build.log
	$(DOCKER) $(PICARD) MarkDuplicates --INPUT $(addprefix /,$(filter %.bam,$^)) --REFERENCE_SEQUENCE $(addprefix /,$(filter %.fasta,$^)) --OUTPUT $(addprefix /,$@) --METRICS_FILE $(addprefix /,$@.metrics) --ASSUME_SORTED true --VALIDATION_STRINGENCY SILENT --COMPRESSION_LEVEL 9 1> $@.stdout 2> $@.stderr

Output/WTS/STAR/%.Sort.MarkDuplicate.bam: Bash/remove_CR.bash Output/WTS/STAR/%.Sort.dm-tmp.bam | build.log
	$(DOCKER) $(BASH) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/WTS/STAR/%.Sort.MarkDuplicate.bai: Output/WTS/STAR/%.Sort.MarkDuplicate.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/WTS/STAR/%.Sort.MarkDuplicate.BQSR.table: $(REFERENCE_FASTA) $(REFERENCE_KNOWN_VCF) Output/WTS/STAR/%.Sort.MarkDuplicate.bam Output/WTS/STAR/%.Sort.MarkDuplicate.bai | build.log
	$(DOCKER) $(GATK) BaseRecalibrator --input $(addprefix /,$(filter %.bam,$^)) $(foreach site,$(addprefix /,$(filter %.vcf,$^)),--known-sites $(site)) --output $(addprefix /,$@) --reference $(addprefix /,$(filter %.fasta,$^)) 1> $@.stdout 2> $@.stderr

Output/WTS/STAR/%.Sort.MarkDuplicate.BQSR.bam: $(REFERENCE_FASTA) Output/WTS/STAR/%.Sort.MarkDuplicate.bam Output/WTS/STAR/%.Sort.MarkDuplicate.BQSR.table | build.log
	$(DOCKER) $(GATK) ApplyBQSR --bqsr-recal-file $(addprefix /,$(filter %.table,$^)) --input $(addprefix /,$(filter %.bam,$^)) --output $(addprefix /,$@) --reference $(addprefix /,$(filter %.fasta,$^)) 1> $@.stdout 2> $@.stderr

Output/WTS/STAR/%.Sort.MarkDuplicate.BQSR.bai: Output/WTS/STAR/%.Sort.MarkDuplicate.BQSR.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step58: Output/WTS/STAR/Reference/Genome $(addprefix Output/WTS/STAR/,$(addsuffix /Aligned.out.sam,$(WTS_ID))) $(addprefix Output/WTS/STAR/,$(addsuffix .bam,$(WTS_ID))) $(addprefix Output/WTS/STAR/,$(addsuffix .Sort.bam,$(WTS_ID))) $(addprefix Output/WTS/STAR/,$(addsuffix .Sort.bai,$(WTS_ID))) $(addprefix Output/WTS/STAR/,$(addsuffix .Sort.dm-tmp.bam,$(WTS_ID))) $(addprefix Output/WTS/STAR/,$(addsuffix .Sort.MarkDuplicate.bam,$(WTS_ID))) $(addprefix Output/WTS/STAR/,$(addsuffix .Sort.MarkDuplicate.bai,$(WTS_ID))) $(addprefix Output/WTS/STAR/,$(addsuffix .Sort.MarkDuplicate.BQSR.table,$(WTS_ID))) $(addprefix Output/WTS/STAR/,$(addsuffix .Sort.MarkDuplicate.BQSR.bam,$(WTS_ID))) $(addprefix Output/WTS/STAR/,$(addsuffix .Sort.MarkDuplicate.BQSR.bai,$(WTS_ID)))
.PHONY: step58

# Step 59 (WTS Bowtie2 alignment with unmapped)
Output/WTS/Bowtie2/Reference:
	mkdir -p $@

Output/WTS/Bowtie2/Reference/Bowtie2.1.bt2: $(REFERENCE_FASTA) | Output/WTS/Bowtie2/Reference build.log
	$(DOCKER) bowtie2-build --threads $(CPUS) --seed "12345" $(addprefix /,$(filter %.fasta,$^)) $(addprefix /,$(basename $@)) 1> $@.stdout 2> $@.stderr

Output/WTS/Bowtie2/%.sam: Output/WTS/Bowtie2/Reference/Bowtie2.1.bt2 Data/WTS/%_*.fastq.gz | build.log
	$(DOCKER) bowtie2 -q --threads $(CPUS) --qc-filter --seed "12345" --rg-id $* --rg "ID:$*\tPL:ILLUMINA\tLB:$*\tSM:$*\tCN:UNIST" --time --met-stderr -x $(addprefix /,$(basename $(basename $(filter %.bt2,$^)))) -1 $(addprefix /,$(word 1,$(sort $(filter %.fastq.gz,$^)))) -2 $(addprefix /,$(word 2,$(sort $(filter %.fastq.gz,$^)))) -S $(addprefix /,$@) 1> $@.stdout 2> $@.stderr

Output/WTS/Bowtie2/%.bam: $(REFERENCE_FASTA) Output/WTS/Bowtie2/%.sam | build.log
	$(DOCKER) samtools view -b -h --threads $(CPUS) --reference $(addprefix /,$(filter %.fasta,$^)) $(addprefix /,$(filter %.sam,$^)) -o $(addprefix /,$@)

Output/WTS/Bowtie2/%.Sort.bam: $(REFERENCE_FASTA) Output/WTS/Bowtie2/%.bam | build.log
	$(DOCKER) $(PICARD) AddOrReplaceReadGroups --INPUT $(addprefix /,$(filter %.bam,$^)) --OUTPUT $(addprefix /,$@) --RGLB "Library" --RGPL "Illumina" --RGPU "PlatformUnit" --RGSM "$*" --RGID "$*" --COMPRESSION_LEVEL 9 --REFERENCE_SEQUENCE $(addprefix /,$(filter %.fasta,$^)) --SORT_ORDER "coordinate" --VALIDATION_STRINGENCY "SILENT" 1> $@.stdout 2> $@.stderr

Output/WTS/Bowtie2/%.Sort.bai: Output/WTS/Bowtie2/%.Sort.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/WTS/Bowtie2/%.Sort.dm-tmp.bam: $(REFERENCE_FASTA) Output/WTS/Bowtie2/%.Sort.bam Output/WTS/Bowtie2/%.Sort.bai | build.log
	$(DOCKER) $(PICARD) MarkDuplicates --INPUT $(addprefix /,$(filter %.bam,$^)) --REFERENCE_SEQUENCE $(addprefix /,$(filter %.fasta,$^)) --OUTPUT $(addprefix /,$@) --METRICS_FILE $(addprefix /,$@.metrics) --ASSUME_SORTED true --VALIDATION_STRINGENCY SILENT --COMPRESSION_LEVEL 9 1> $@.stdout 2> $@.stderr

Output/WTS/Bowtie2/%.Sort.MarkDuplicate.bam: Bash/remove_CR.bash Output/WTS/Bowtie2/%.Sort.dm-tmp.bam | build.log
	$(DOCKER) $(BASH) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/WTS/Bowtie2/%.Sort.MarkDuplicate.bai: Output/WTS/Bowtie2/%.Sort.MarkDuplicate.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/WTS/Bowtie2/%.Sort.MarkDuplicate.BQSR.table: $(REFERENCE_FASTA) $(REFERENCE_KNOWN_VCF) Output/WTS/Bowtie2/%.Sort.MarkDuplicate.bam Output/WTS/Bowtie2/%.Sort.MarkDuplicate.bai | build.log
	$(DOCKER) $(GATK) BaseRecalibrator --input $(addprefix /,$(filter %.bam,$^)) $(foreach site,$(addprefix /,$(filter %.vcf,$^)),--known-sites $(site)) --output $(addprefix /,$@) --reference $(addprefix /,$(filter %.fasta,$^)) 1> $@.stdout 2> $@.stderr

Output/WTS/Bowtie2/%.Sort.MarkDuplicate.BQSR.bam: $(REFERENCE_FASTA) Output/WTS/Bowtie2/%.Sort.MarkDuplicate.bam Output/WTS/Bowtie2/%.Sort.MarkDuplicate.BQSR.table | build.log
	$(DOCKER) $(GATK) ApplyBQSR --bqsr-recal-file $(addprefix /,$(filter %.table,$^)) --input $(addprefix /,$(filter %.bam,$^)) --output $(addprefix /,$@) --reference $(addprefix /,$(filter %.fasta,$^)) 1> $@.stdout 2> $@.stderr

Output/WTS/Bowtie2/%.Sort.MarkDuplicate.BQSR.bai: Output/WTS/Bowtie2/%.Sort.MarkDuplicate.BQSR.bam | build.log
	$(DOCKER) samtools index -b -@ $(CPUS) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

step59: Output/WTS/Bowtie2/Reference/Bowtie2.1.bt2 $(addprefix Output/WTS/Bowtie2/,$(addsuffix .sam,$(WTS_ID))) $(addprefix Output/WTS/Bowtie2/,$(addsuffix .bam,$(WTS_ID))) $(addprefix Output/WTS/Bowtie2/,$(addsuffix .Sort.bam,$(WTS_ID))) $(addprefix Output/WTS/Bowtie2/,$(addsuffix .Sort.bai,$(WTS_ID))) $(addprefix Output/WTS/Bowtie2/,$(addsuffix .Sort.dm-tmp.bam,$(WTS_ID))) $(addprefix Output/WTS/Bowtie2/,$(addsuffix .Sort.MarkDuplicate.bam,$(WTS_ID))) $(addprefix Output/WTS/Bowtie2/,$(addsuffix .Sort.MarkDuplicate.bai,$(WTS_ID))) $(addprefix Output/WTS/Bowtie2/,$(addsuffix .Sort.MarkDuplicate.BQSR.table,$(WTS_ID))) $(addprefix Output/WTS/Bowtie2/,$(addsuffix .Sort.MarkDuplicate.BQSR.bam,$(WTS_ID))) $(addprefix Output/WTS/Bowtie2/,$(addsuffix .Sort.MarkDuplicate.BQSR.bai,$(WTS_ID)))
.PHONY: step59

# Step 60 (PathSeq - WTS)
Output/Pathseq/WTS/Bowtie2 Output/Pathseq/WTS/STAR:
	mkdir -p $@

Output/Pathseq/WTS/Bowtie2/%.tsv Output/Pathseq/WTS/Bowtie2/%.bam: Output/WTS/Bowtie2/%.Sort.MarkDuplicate.BQSR.bam Output/Pathseq/Reference/pathseq_microbe.fa.img Output/Pathseq/Reference/pathseq_microbe.dict Output/Pathseq/Reference/pathseq_taxonomy.db $(REFERENCE_FASTA) Output/Pathseq/Reference/pathseq_host.bfi Output/Pathseq/Reference/pathseq_host.fa.img | Output/Pathseq/WTS/Bowtie2
	$(DOCKER) $(GATK) PathSeqPipelineSpark --input $(addprefix /,$(word 1,$^)) --microbe-bwa-image $(addprefix /,$(word 2,$^)) --microbe-dict $(addprefix /,$(word 3,$^)) --taxonomy-file $(addprefix /,$(word 4,$^)) --reference $(addprefix /,$(word 5,$^)) --kmer-file $(addprefix /,$(word 6,$^)) --filter-bwa-image $(addprefix /,$(word 7,$^)) --scores-output $(addprefix /,$(@D)/$*.tsv) --output $(addprefix /,$(@D)/$*.bam) --is-host-aligned true 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR/%.tsv Output/Pathseq/WTS/STAR/%.bam: Output/WTS/STAR/%.Sort.MarkDuplicate.BQSR.bam Output/Pathseq/Reference/pathseq_microbe.fa.img Output/Pathseq/Reference/pathseq_microbe.dict Output/Pathseq/Reference/pathseq_taxonomy.db $(REFERENCE_FASTA) Output/Pathseq/Reference/pathseq_host.bfi Output/Pathseq/Reference/pathseq_host.fa.img | Output/Pathseq/WTS/STAR
	$(DOCKER) $(GATK) PathSeqPipelineSpark --input $(addprefix /,$(word 1,$^)) --microbe-bwa-image $(addprefix /,$(word 2,$^)) --microbe-dict $(addprefix /,$(word 3,$^)) --taxonomy-file $(addprefix /,$(word 4,$^)) --reference $(addprefix /,$(word 5,$^)) --kmer-file $(addprefix /,$(word 6,$^)) --filter-bwa-image $(addprefix /,$(word 7,$^)) --scores-output $(addprefix /,$(@D)/$*.tsv) --output $(addprefix /,$(@D)/$*.bam) --is-host-aligned true 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR.%.tsv: Python/aggregate_PathSeq_tsv.py $(addprefix Output/Pathseq/WTS/STAR/,$(addsuffix .tsv,$(WTS_ID))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2.%.tsv: Python/aggregate_PathSeq_tsv.py $(addprefix Output/Pathseq/WTS/Bowtie2/,$(addsuffix .tsv,$(WTS_ID))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-STAR-SQC.%.pdf: Python/aggregate_PathSeq.py Output/Pathseq/WTS/STAR.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-STAR-ADC.%.pdf: Python/aggregate_PathSeq.py Output/Pathseq/WTS/STAR.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-Bowtie2-SQC.%.pdf: Python/aggregate_PathSeq.py Output/Pathseq/WTS/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-Bowtie2-ADC.%.pdf: Python/aggregate_PathSeq.py Output/Pathseq/WTS/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-STAR-SQC.%.tar: Python/aggregate_PathSeq_violin.py Output/Pathseq/WTS/STAR.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-STAR-ADC.%.tar: Python/aggregate_PathSeq_violin.py Output/Pathseq/WTS/STAR.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-Bowtie2-SQC.%.tar: Python/aggregate_PathSeq_violin.py Output/Pathseq/WTS/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-Bowtie2-ADC.%.tar: Python/aggregate_PathSeq_violin.py Output/Pathseq/WTS/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-STAR-SQC.Smoking.%.tar: Python/aggregate_PathSeq_violin_clinical.py Output/Pathseq/WTS/STAR.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-STAR-SQC.Recurrence.%.tar: Python/aggregate_PathSeq_violin_clinical.py Output/Pathseq/WTS/STAR.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --level "$*" --compare "Recurrence" "NO" "YES" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-STAR-ADC.Smoking.%.tar: Python/aggregate_PathSeq_violin_clinical.py Output/Pathseq/WTS/STAR.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-STAR-ADC.Recurrence.%.tar: Python/aggregate_PathSeq_violin_clinical.py Output/Pathseq/WTS/STAR.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --level "$*" --compare "Recurrence" "NO" "YES" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-Bowtie2-SQC.Smoking.%.tar: Python/aggregate_PathSeq_violin_clinical.py Output/Pathseq/WTS/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-Bowtie2-SQC.Recurrence.%.tar: Python/aggregate_PathSeq_violin_clinical.py Output/Pathseq/WTS/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --level "$*" --compare "Recurrence" "NO" "YES" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-Bowtie2-ADC.Smoking.%.tar: Python/aggregate_PathSeq_violin_clinical.py Output/Pathseq/WTS/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-Bowtie2-ADC.Recurrence.%.tar: Python/aggregate_PathSeq_violin_clinical.py Output/Pathseq/WTS/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --level "$*" --compare "Recurrence" "NO" "YES" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-STAR-SQC.MutationSharedProportion-Median.%.tar: Python/aggregate_PathSeq_violin_MutationSharedProportion.py Output/Pathseq/WTS/STAR.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --median --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-STAR-SQC.MutationSharedProportion-Mean.%.tar: Python/aggregate_PathSeq_violin_MutationSharedProportion.py Output/Pathseq/WTS/STAR.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --mean --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-Bowtie2-SQC.MutationSharedProportion-Median.%.tar: Python/aggregate_PathSeq_violin_MutationSharedProportion.py Output/Pathseq/WTS/Bowtie2.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --median --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-Bowtie2-SQC.MutationSharedProportion-Mean.%.tar: Python/aggregate_PathSeq_violin_MutationSharedProportion.py Output/Pathseq/WTS/Bowtie2.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --mean --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-STAR-ADC.MutationSharedProportion-Median.%.tar: Python/aggregate_PathSeq_violin_MutationSharedProportion.py Output/Pathseq/WTS/STAR.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --median --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-STAR-ADC.MutationSharedProportion-Mean.%.tar: Python/aggregate_PathSeq_violin_MutationSharedProportion.py Output/Pathseq/WTS/STAR.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --mean --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-Bowtie2-ADC.MutationSharedProportion-Median.%.tar: Python/aggregate_PathSeq_violin_MutationSharedProportion.py Output/Pathseq/WTS/Bowtie2.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --median --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Violin-Bowtie2-ADC.MutationSharedProportion-Mean.%.tar: Python/aggregate_PathSeq_violin_MutationSharedProportion.py Output/Pathseq/WTS/Bowtie2.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --mean --level "$*" --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Cluster-STAR-SQC.%.pdf: Python/aggregate_PathSeq_clustermap.py Output/Pathseq/WTS/STAR.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Cluster-STAR-ADC.%.pdf: Python/aggregate_PathSeq_clustermap.py Output/Pathseq/WTS/STAR.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Cluster-Bowtie2-SQC.%.pdf: Python/aggregate_PathSeq_clustermap.py Output/Pathseq/WTS/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Cluster-Bowtie2-ADC.%.pdf: Python/aggregate_PathSeq_clustermap.py Output/Pathseq/WTS/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-STAR-SQC.Smoking.%.pdf: Python/aggregate_PathSeq_clinical.py Output/Pathseq/WTS/STAR.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-STAR-ADC.Smoking.%.pdf: Python/aggregate_PathSeq_clinical.py Output/Pathseq/WTS/STAR.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-Bowtie2-SQC.Smoking.%.pdf: Python/aggregate_PathSeq_clinical.py Output/Pathseq/WTS/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-Bowtie2-ADC.Smoking.%.pdf: Python/aggregate_PathSeq_clinical.py Output/Pathseq/WTS/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-STAR-SQC.Recurrence.%.pdf: Python/aggregate_PathSeq_clinical.py Output/Pathseq/WTS/STAR.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-STAR-ADC.Recurrence.%.pdf: Python/aggregate_PathSeq_clinical.py Output/Pathseq/WTS/STAR.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-Bowtie2-SQC.Recurrence.%.pdf: Python/aggregate_PathSeq_clinical.py Output/Pathseq/WTS/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --type --level "$*" --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-Bowtie2-ADC.Recurrence.%.pdf: Python/aggregate_PathSeq_clinical.py Output/Pathseq/WTS/Bowtie2.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --type --level "$*" --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-STAR-SQC.MutationSharedProportion-Median.%.pdf: Python/aggregate_PathSeq_MutationSharedProportion.py Output/Pathseq/WTS/STAR.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --level "$*" --median 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-STAR-SQC.MutationSharedProportion-Mean.%.pdf: Python/aggregate_PathSeq_MutationSharedProportion.py Output/Pathseq/WTS/STAR.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --level "$*" --mean 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-Bowtie2-SQC.MutationSharedProportion-Median.%.pdf: Python/aggregate_PathSeq_MutationSharedProportion.py Output/Pathseq/WTS/Bowtie2.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --level "$*" --median 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-Bowtie2-SQC.MutationSharedProportion-Mean.%.pdf: Python/aggregate_PathSeq_MutationSharedProportion.py Output/Pathseq/WTS/Bowtie2.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --level "$*" --mean 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-STAR-ADC.MutationSharedProportion-Median.%.pdf: Python/aggregate_PathSeq_MutationSharedProportion.py Output/Pathseq/WTS/STAR.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --level "$*" --median 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-STAR-ADC.MutationSharedProportion-Mean.%.pdf: Python/aggregate_PathSeq_MutationSharedProportion.py Output/Pathseq/WTS/STAR.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --level "$*" --mean 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-Bowtie2-ADC.MutationSharedProportion-Median.%.pdf: Python/aggregate_PathSeq_MutationSharedProportion.py Output/Pathseq/WTS/Bowtie2.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --level "$*" --median 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bar-Bowtie2-ADC.MutationSharedProportion-Mean.%.pdf: Python/aggregate_PathSeq_MutationSharedProportion.py Output/Pathseq/WTS/Bowtie2.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --level "$*" --mean 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR-alpha.%.tsv: Python/calculate_alpha_PathSeq.py Output/Pathseq/WTS/STAR.%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2-alpha.%.tsv: Python/calculate_alpha_PathSeq.py Output/Pathseq/WTS/Bowtie2.%.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR.tree.nwk: Python/make_tree_PathSeq.py $(addprefix Output/Pathseq/WTS/STAR/,$(addsuffix .tsv,$(WTS_ID))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2.tree.nwk: Python/make_tree_PathSeq.py $(addprefix Output/Pathseq/WTS/Bowtie2/,$(addsuffix .tsv,$(WTS_ID))) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR-beta.%.Unweighted.tsv: Python/calculate_beta_PathSeq.py Output/Pathseq/WTS/STAR.%.tsv Output/Pathseq/WTS/STAR.tree.nwk | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --beta "unweighted_unifrac" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR-beta.%.Weighted.tsv: Python/calculate_beta_PathSeq.py Output/Pathseq/WTS/STAR.%.tsv Output/Pathseq/WTS/STAR.tree.nwk | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --beta "weighted_unifrac" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2-beta.%.Unweighted.tsv: Python/calculate_beta_PathSeq.py Output/Pathseq/WTS/Bowtie2.%.tsv Output/Pathseq/WTS/Bowtie2.tree.nwk | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --beta "unweighted_unifrac" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2-beta.%.Weighted.tsv: Python/calculate_beta_PathSeq.py Output/Pathseq/WTS/Bowtie2.%.tsv Output/Pathseq/WTS/Bowtie2.tree.nwk | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --beta "weighted_unifrac" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR-SQC-alpha.%.tar: Python/draw_alpha_PathSeq.py Output/Pathseq/WTS/STAR-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2-SQC-alpha.%.tar: Python/draw_alpha_PathSeq.py Output/Pathseq/WTS/Bowtie2-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR-SQC-alpha.Smoking.%.tar: Python/draw_alpha_PathSeq_clinical.py Output/Pathseq/WTS/STAR-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR-SQC-alpha.Recurrence.%.tar: Python/draw_alpha_PathSeq_clinical.py Output/Pathseq/WTS/STAR-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR-ADC-alpha.Smoking.%.tar: Python/draw_alpha_PathSeq_clinical.py Output/Pathseq/WTS/STAR-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR-ADC-alpha.Recurrence.%.tar: Python/draw_alpha_PathSeq_clinical.py Output/Pathseq/WTS/STAR-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2-SQC-alpha.Smoking.%.tar: Python/draw_alpha_PathSeq_clinical.py Output/Pathseq/WTS/Bowtie2-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2-SQC-alpha.Recurrence.%.tar: Python/draw_alpha_PathSeq_clinical.py Output/Pathseq/WTS/Bowtie2-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2-ADC-alpha.Smoking.%.tar: Python/draw_alpha_PathSeq_clinical.py Output/Pathseq/WTS/Bowtie2-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2-ADC-alpha.Recurrence.%.tar: Python/draw_alpha_PathSeq_clinical.py Output/Pathseq/WTS/Bowtie2-alpha.%.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --compare "Recurrence" "NO" "YES" 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR-SQC-alpha.MutationSharedProportion-Median.%.tar: Python/draw_alpha_PathSeq_MutationSharedProportion.py Output/Pathseq/WTS/STAR-alpha.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --median 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR-SQC-alpha.MutationSharedProportion-Mean.%.tar: Python/draw_alpha_PathSeq_MutationSharedProportion.py Output/Pathseq/WTS/STAR-alpha.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --mean 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2-SQC-alpha.MutationSharedProportion-Median.%.tar: Python/draw_alpha_PathSeq_MutationSharedProportion.py Output/Pathseq/WTS/Bowtie2-alpha.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --median 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2-SQC-alpha.MutationSharedProportion-Mean.%.tar: Python/draw_alpha_PathSeq_MutationSharedProportion.py Output/Pathseq/WTS/Bowtie2-alpha.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --mean 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR-ADC-alpha.MutationSharedProportion-Median.%.tar: Python/draw_alpha_PathSeq_MutationSharedProportion.py Output/Pathseq/WTS/STAR-alpha.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --median 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR-ADC-alpha.MutationSharedProportion-Mean.%.tar: Python/draw_alpha_PathSeq_MutationSharedProportion.py Output/Pathseq/WTS/STAR-alpha.%.tsv Output/MutationSharedProportion/BWA.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --mean 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2-ADC-alpha.MutationSharedProportion-Median.%.tar: Python/draw_alpha_PathSeq_MutationSharedProportion.py Output/Pathseq/WTS/Bowtie2-alpha.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --median 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2-ADC-alpha.MutationSharedProportion-Mean.%.tar: Python/draw_alpha_PathSeq_MutationSharedProportion.py Output/Pathseq/WTS/Bowtie2-alpha.%.tsv Output/MutationSharedProportion/Bowtie2.tsv | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --mean 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR-SQC-beta.%.Unweighted.pdf: Python/draw_beta_PathSeq.py Output/Pathseq/WTS/STAR-beta.%.Unweighted.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR-SQC-beta.%.Weighted.pdf: Python/draw_beta_PathSeq.py Output/Pathseq/WTS/STAR-beta.%.Weighted.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR-ADC-beta.%.Unweighted.pdf: Python/draw_beta_PathSeq.py Output/Pathseq/WTS/STAR-beta.%.Unweighted.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/STAR-ADC-beta.%.Weighted.pdf: Python/draw_beta_PathSeq.py Output/Pathseq/WTS/STAR-beta.%.Weighted.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2-SQC-beta.%.Unweighted.pdf: Python/draw_beta_PathSeq.py Output/Pathseq/WTS/Bowtie2-beta.%.Unweighted.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2-SQC-beta.%.Weighted.pdf: Python/draw_beta_PathSeq.py Output/Pathseq/WTS/Bowtie2-beta.%.Weighted.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2-ADC-beta.%.Unweighted.pdf: Python/draw_beta_PathSeq.py Output/Pathseq/WTS/Bowtie2-beta.%.Unweighted.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Pathseq/WTS/Bowtie2-ADC-beta.%.Weighted.pdf: Python/draw_beta_PathSeq.py Output/Pathseq/WTS/Bowtie2-beta.%.Weighted.tsv $(CLINICAL_DATA_CSV) | build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step60-Bowtie2-Merge: $(addprefix Output/Pathseq/WTS/Bowtie2/,$(addsuffix .tsv,$(WTS_ID))) $(addprefix Output/Pathseq/WTS/Bowtie2.,$(addsuffix .tsv,$(TAXONOMY_TYPES)))
step60-Bowite2-Bar: $(addprefix Output/Pathseq/WTS/Bar-Bowtie2-SQC.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-Bowtie2-ADC.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-Bowtie2-ADC.Smoking.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-Bowtie2-SQC.Recurrence.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-Bowtie2-ADC.Recurrence.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-Bowtie2-SQC.Smoking.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-Bowtie2-SQC.MutationSharedProportion-Median.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-Bowtie2-SQC.MutationSharedProportion-Mean.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-Bowtie2-ADC.MutationSharedProportion-Median.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-Bowtie2-ADC.MutationSharedProportion-Mean.,$(addsuffix .pdf,$(TAXONOMY_TYPES)))
step60-Bowtie2-Violin: $(addprefix Output/Pathseq/WTS/Violin-Bowtie2-SQC.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-Bowtie2-ADC.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-Bowtie2-SQC.Smoking.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-Bowtie2-ADC.Smoking.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-Bowtie2-SQC.Recurrence.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-Bowtie2-ADC.Recurrence.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-Bowtie2-SQC.MutationSharedProportion-Median.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-Bowtie2-SQC.MutationSharedProportion-Mean.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-Bowtie2-ADC.MutationSharedProportion-Median.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-Bowtie2-ADC.MutationSharedProportion-Mean.,$(addsuffix .tar,$(TAXONOMY_TYPES)))
step60-Bowtie2-Cluster: $(addprefix Output/Pathseq/WTS/Cluster-Bowtie2-SQC.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Cluster-Bowtie2-ADC.,$(addsuffix .pdf,$(TAXONOMY_TYPES)))
step60-Bowtie2-Diversity: $(addprefix Output/Pathseq/WTS/Bowtie2-alpha.,$(addsuffix .tsv,$(TAXONOMY_TYPES))) Output/Pathseq/WTS/Bowtie2.tree.nwk $(addprefix Output/Pathseq/WTS/Bowtie2-beta.,$(addsuffix .Unweighted.tsv,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bowtie2-beta.,$(addsuffix .Weighted.tsv,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bowtie2-SQC-alpha.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bowtie2-SQC-alpha.Smoking.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bowtie2-SQC-alpha.Recurrence.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bowtie2-ADC-alpha.Smoking.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bowtie2-ADC-alpha.Recurrence.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bowtie2-SQC-beta.,$(addsuffix .Unweighted.pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bowtie2-SQC-beta.,$(addsuffix .Weighted.pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bowtie2-ADC-beta.,$(addsuffix .Unweighted.pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bowtie2-ADC-beta.,$(addsuffix .Weighted.pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bowtie2-SQC-alpha.MutationSharedProportion-Median.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bowtie2-SQC-alpha.MutationSharedProportion-Mean.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bowtie2-ADC-alpha.MutationSharedProportion-Median.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bowtie2-ADC-alpha.MutationSharedProportion-Mean.,$(addsuffix .tar,$(TAXONOMY_TYPES)))
step60-Bowtie2: step60-Bowtie2-Merge step60-Bowite2-Bar step60-Bowtie2-Violin step60-Bowtie2-Cluster step60-Bowtie2-Diversity
step60-STAR-Merge: $(addprefix Output/Pathseq/WTS/STAR/,$(addsuffix .tsv,$(WTS_ID))) $(addprefix Output/Pathseq/WTS/STAR.,$(addsuffix .tsv,$(TAXONOMY_TYPES)))
step60-Bowite2-Bar: $(addprefix Output/Pathseq/WTS/Bar-STAR-SQC.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-STAR-ADC.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-STAR-ADC.Smoking.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-STAR-SQC.Recurrence.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-STAR-ADC.Recurrence.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-STAR-SQC.Smoking.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-STAR-SQC.MutationSharedProportion-Median.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-STAR-SQC.MutationSharedProportion-Mean.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-STAR-ADC.MutationSharedProportion-Median.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Bar-STAR-ADC.MutationSharedProportion-Mean.,$(addsuffix .pdf,$(TAXONOMY_TYPES)))
step60-STAR-Violin: $(addprefix Output/Pathseq/WTS/Violin-STAR-SQC.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-STAR-ADC.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-STAR-SQC.Smoking.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-STAR-ADC.Smoking.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-STAR-SQC.Recurrence.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-STAR-ADC.Recurrence.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-STAR-SQC.MutationSharedProportion-Median.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-STAR-SQC.MutationSharedProportion-Mean.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-STAR-ADC.MutationSharedProportion-Median.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Violin-STAR-ADC.MutationSharedProportion-Mean.,$(addsuffix .tar,$(TAXONOMY_TYPES)))
step60-STAR-Cluster: $(addprefix Output/Pathseq/WTS/Cluster-STAR-SQC.,$(addsuffix .pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/Cluster-STAR-ADC.,$(addsuffix .pdf,$(TAXONOMY_TYPES)))
step60-STAR-Diversity: $(addprefix Output/Pathseq/WTS/STAR-alpha.,$(addsuffix .tsv,$(TAXONOMY_TYPES))) Output/Pathseq/WTS/STAR.tree.nwk $(addprefix Output/Pathseq/WTS/STAR-beta.,$(addsuffix .Unweighted.tsv,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/STAR-beta.,$(addsuffix .Weighted.tsv,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/STAR-SQC-alpha.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/STAR-SQC-alpha.Smoking.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/STAR-SQC-alpha.Recurrence.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/STAR-ADC-alpha.Smoking.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/STAR-ADC-alpha.Recurrence.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/STAR-SQC-beta.,$(addsuffix .Unweighted.pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/STAR-SQC-beta.,$(addsuffix .Weighted.pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/STAR-ADC-beta.,$(addsuffix .Unweighted.pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/STAR-ADC-beta.,$(addsuffix .Weighted.pdf,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/STAR-SQC-alpha.MutationSharedProportion-Median.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/STAR-SQC-alpha.MutationSharedProportion-Mean.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/STAR-ADC-alpha.MutationSharedProportion-Median.,$(addsuffix .tar,$(TAXONOMY_TYPES))) $(addprefix Output/Pathseq/WTS/STAR-ADC-alpha.MutationSharedProportion-Mean.,$(addsuffix .tar,$(TAXONOMY_TYPES)))
step60-STAR: step60-STAR-Merge step60-Bowite2-Bar step60-STAR-Violin step60-STAR-Cluster step60-STAR-Diversity
step60: step60-Bowtie2-Merge step60-Bowite2-Bar step60-Bowtie2-Violin step60-Bowtie2-Cluster step60-Bowtie2-Diversity step60-Bowtie2 step60-STAR-Merge step60-Bowite2-Bar step60-STAR-Violin step60-STAR-Cluster step60-STAR-Diversity step60-STAR
.PHONY: step60-Bowtie2-Merge step60-Bowite2-Bar step60-Bowtie2-Violin step60-Bowtie2-Cluster step60-Bowtie2-Diversity step60-Bowtie2 step60-STAR-Merge step60-Bowite2-Bar step60-STAR-Violin step60-STAR-Cluster step60-STAR-Diversity step60-STAR step60

# Step 61 (Short variants analysis upon clinical data with Mutation Sharing)
Output/Gene_Clinical_Sharing:
	mkdir -p $@

Output/Gene_Clinical_Sharing/BWA-SQC.Mean.tsv Output/Gene_Clinical_Sharing/BWA-SQC.Mean.pdf: Python/plot_gene_sharing.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/BWA.tsv $(GENE_CENCUS_ONE_CSV) | build.log Output/Gene_Clinical_Sharing
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --SQC --mean 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Sharing/BWA-SQC.Median.tsv Output/Gene_Clinical_Sharing/BWA-SQC.Median.pdf: Python/plot_gene_sharing.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/BWA.tsv  $(GENE_CENCUS_ONE_CSV) | build.log Output/Gene_Clinical_Sharing
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --SQC --median 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Sharing/BWA-ADC.Mean.tsv Output/Gene_Clinical_Sharing/BWA-ADC.Mean.pdf: Python/plot_gene_sharing.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/BWA.tsv $(GENE_CENCUS_ONE_CSV) | build.log Output/Gene_Clinical_Sharing
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --ADC --mean 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Sharing/BWA-ADC.Median.tsv Output/Gene_Clinical_Sharing/BWA-ADC.Median.pdf: Python/plot_gene_sharing.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/BWA.tsv $(GENE_CENCUS_ONE_CSV) | build.log Output/Gene_Clinical_Sharing
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --ADC --median 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Sharing/Bowtie2-SQC.Mean.tsv Output/Gene_Clinical_Sharing/Bowtie2-SQC.Mean.pdf: Python/plot_gene_sharing.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/Bowtie2.tsv $(GENE_CENCUS_ONE_CSV) | build.log Output/Gene_Clinical_Sharing
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --SQC --mean 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Sharing/Bowtie2-SQC.Median.tsv Output/Gene_Clinical_Sharing/Bowtie2-SQC.Median.pdf: Python/plot_gene_sharing.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/Bowtie2.tsv $(GENE_CENCUS_ONE_CSV) | build.log Output/Gene_Clinical_Sharing
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --SQC --median 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Sharing/Bowtie2-ADC.Mean.tsv Output/Gene_Clinical_Sharing/Bowtie2-ADC.Mean.pdf: Python/plot_gene_sharing.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/Bowtie2.tsv $(GENE_CENCUS_ONE_CSV) | build.log Output/Gene_Clinical_Sharing
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --ADC --mean 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Sharing/Bowtie2-ADC.Median.tsv Output/Gene_Clinical_Sharing/Bowtie2-ADC.Median.pdf: Python/plot_gene_sharing.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/Bowtie2.tsv $(GENE_CENCUS_ONE_CSV) | build.log Output/Gene_Clinical_Sharing
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --ADC --median 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Sharing/BWA-SQC.Mean-Both.tsv Output/Gene_Clinical_Sharing/BWA-SQC.Mean-Both.pdf: Python/plot_gene_sharing.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/BWA.tsv $(GENE_CENCUS_BOTH_CSV) | build.log Output/Gene_Clinical_Sharing
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --SQC --mean 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Sharing/BWA-SQC.Median-Both.tsv Output/Gene_Clinical_Sharing/BWA-SQC.Median-Both.pdf: Python/plot_gene_sharing.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/BWA.tsv  $(GENE_CENCUS_BOTH_CSV) | build.log Output/Gene_Clinical_Sharing
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --SQC --median 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Sharing/BWA-ADC.Mean-Both.tsv Output/Gene_Clinical_Sharing/BWA-ADC.Mean-Both.pdf: Python/plot_gene_sharing.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/BWA.tsv $(GENE_CENCUS_BOTH_CSV) | build.log Output/Gene_Clinical_Sharing
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --ADC --mean 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Sharing/BWA-ADC.Median-Both.tsv Output/Gene_Clinical_Sharing/BWA-ADC.Median-Both.pdf: Python/plot_gene_sharing.py $(addprefix Output/Mutect2/BWA/,$(addsuffix .BWA.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/BWA.tsv $(GENE_CENCUS_BOTH_CSV) | build.log Output/Gene_Clinical_Sharing
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --ADC --median 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Sharing/Bowtie2-SQC.Mean-Both.tsv Output/Gene_Clinical_Sharing/Bowtie2-SQC.Mean-Both.pdf: Python/plot_gene_sharing.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/Bowtie2.tsv $(GENE_CENCUS_BOTH_CSV) | build.log Output/Gene_Clinical_Sharing
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --SQC --mean 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Sharing/Bowtie2-SQC.Median-Both.tsv Output/Gene_Clinical_Sharing/Bowtie2-SQC.Median-Both.pdf: Python/plot_gene_sharing.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/Bowtie2.tsv $(GENE_CENCUS_BOTH_CSV) | build.log Output/Gene_Clinical_Sharing
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --SQC --median 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Sharing/Bowtie2-ADC.Mean-Both.tsv Output/Gene_Clinical_Sharing/Bowtie2-ADC.Mean-Both.pdf: Python/plot_gene_sharing.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/Bowtie2.tsv $(GENE_CENCUS_BOTH_CSV) | build.log Output/Gene_Clinical_Sharing
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --ADC --mean 1> $@.stdout 2> $@.stderr

Output/Gene_Clinical_Sharing/Bowtie2-ADC.Median-Both.tsv Output/Gene_Clinical_Sharing/Bowtie2-ADC.Median-Both.pdf: Python/plot_gene_sharing.py $(addprefix Output/Mutect2/Bowtie2/,$(addsuffix .Bowtie2.PASS.maf,$(WES_TUMOR_ID))) Output/MutationSharedProportion/Bowtie2.tsv $(GENE_CENCUS_BOTH_CSV) | build.log Output/Gene_Clinical_Sharing
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $(basename $@).tsv $(basename $@).pdf) --cpus $(CPUS) --ADC --median 1> $@.stdout 2> $@.stderr

step61-BWA: Output/Gene_Clinical_Sharing/BWA-SQC.Mean.pdf Output/Gene_Clinical_Sharing/BWA-SQC.Median.pdf Output/Gene_Clinical_Sharing/BWA-SQC.Mean-Both.pdf Output/Gene_Clinical_Sharing/BWA-SQC.Median-Both.pdf
step61-Bowtie2: Output/Gene_Clinical_Sharing/Bowtie2-SQC.Mean.pdf Output/Gene_Clinical_Sharing/Bowtie2-SQC.Median.pdf Output/Gene_Clinical_Sharing/Bowtie2-SQC.Mean-Both.pdf Output/Gene_Clinical_Sharing/Bowtie2-SQC.Median-Both.pdf
step61: step61-BWA step61-Bowtie2
.PHONY: step61-BWA step61-Bowtie2 step61

# Step 62 (Joint plot with Mutation Shared Proportion)
Output/Clinical_vs_Joint-Plot:
	mkdir -p $@

Output/Clinical_vs_Joint-Plot/BWA-SQC-SharedProportion.tar: Python/draw_Joint_SharedProportion.py Output/MutationSharedProportion/BWA.tsv | build.log Output/Clinical_vs_Joint-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Joint-Plot/BWA-ADC-SharedProportion.tar: Python/draw_Joint_SharedProportion.py Output/MutationSharedProportion/BWA.tsv | build.log Output/Clinical_vs_Joint-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Joint-Plot/Bowtie2-SQC-SharedProportion.tar: Python/draw_Joint_SharedProportion.py Output/MutationSharedProportion/Bowtie2.tsv | build.log Output/Clinical_vs_Joint-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --SQC 1> $@.stdout 2> $@.stderr

Output/Clinical_vs_Joint-Plot/Bowtie2-ADC-SharedProportion.tar: Python/draw_Joint_SharedProportion.py Output/MutationSharedProportion/Bowtie2.tsv | build.log Output/Clinical_vs_Joint-Plot
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --ADC 1> $@.stdout 2> $@.stderr

step62: Output/Clinical_vs_Joint-Plot/BWA-SQC-SharedProportion.tar Output/Clinical_vs_Joint-Plot/BWA-ADC-SharedProportion.tar Output/Clinical_vs_Joint-Plot/Bowtie2-SQC-SharedProportion.tar Output/Clinical_vs_Joint-Plot/Bowtie2-ADC-SharedProportion.tar
.PHONY: step62

# Step 63 (Helmsman)
Output/Helmsman/BWA Output/Helmsman/Bowtie2:
	mkdir -p $@

Output/Helmsman/BWA/%/subtype_count_matrix_spectra.txt Output/Helmsman/BWA/%/subtype_count_matrix.txt: $(REFERENCE_FASTA) Output/Mutect2/BWA/%.BWA.PASS.maf | Output/Helmsman/BWA build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D)) && mkdir -p $(@)
	$(DOCKER) $(HELMSMAN) --cpus $(CPUS) --mode "maf" --seed 12345 --verbose --input $(addprefix /,$(filter %.maf,$^)) --fastafile $(addprefix /,$(filter %.fasta,$^)) --impute --projectdir $(@D) 1> $(@D).stdout 2> $(@D).stderr

Output/Helmsman/Bowtie2/%/subtype_count_matrix_spectra.txt Output/Helmsman/Bowtie2/%/subtype_count_matrix.txt: $(REFERENCE_FASTA) Output/Mutect2/Bowtie2/%.Bowtie2.PASS.maf | Output/Helmsman/Bowtie2 build.log
	$(DOCKER) rm -rf $(addprefix /,$(@D)) && mkdir -p $(@)
	$(DOCKER) $(HELMSMAN) --cpus $(CPUS) --mode "maf" --seed 12345 --verbose --input $(addprefix /,$(filter %.maf,$^)) --fastafile $(addprefix /,$(filter %.fasta,$^)) --impute --projectdir $(@D) 1> $(@D).stdout 2> $(@D).stderr

step63-BWA: $(addprefix Output/Helmsman/BWA/,$(addsuffix /subtype_count_matrix_spectra.txt,$(WES_TUMOR_ID)))
step63-Bowtie2: $(addprefix Output/Helmsman/Bowtie2/,$(addsuffix /subtype_count_matrix_spectra.txt,$(WES_TUMOR_ID)))
step63: step63-BWA step63-Bowtie2
.PHONY: step63-BWA step63-Bowtie2 step63

# Step 64 (Histogram of file size)
Output/Size:
	mkdir -p $@

Output/Size/WES.pdf: Python/analysis_filesize_distribution.py Data/RawData/*.fastq.gz | Output/Size build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

Output/Size/WTS.pdf: Python/analysis_filesize_distribution.py Data/WTS/*.fastq.gz | Output/Size build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) 1> $@.stdout 2> $@.stderr

step64: Output/Size/WES.pdf Output/Size/WTS.pdf
.PHONY: step64

# Step 65 (Mutation Shared Proportion vs. Sequenza)
Output/MutationSharedProportion_vs_CNV/Sequenza:
	mkdir -p $@

Output/MutationSharedProportion_vs_CNV/Sequenza/BWA-SQC.tar: Python/compare_MutationSharedProportion_Sequenza.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_REAL_TUMOR_ID))) Output/MutationSharedProportion/BWA.tsv | build.log Output/MutationSharedProportion_vs_CNV/Sequenza
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/MutationSharedProportion_vs_CNV/Sequenza/BWA-ADC.tar: Python/compare_MutationSharedProportion_Sequenza.py $(addprefix Output/Sequenza/BWA/,$(addsuffix /sample_alternative_solutions.txt,$(WES_REAL_TUMOR_ID))) Output/MutationSharedProportion/BWA.tsv | build.log Output/MutationSharedProportion_vs_CNV/Sequenza
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

Output/MutationSharedProportion_vs_CNV/Sequenza/Bowtie2-SQC.tar: Python/compare_MutationSharedProportion_Sequenza.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_REAL_TUMOR_ID))) Output/MutationSharedProportion/Bowtie2.tsv | build.log Output/MutationSharedProportion_vs_CNV/Sequenza
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --SQC 1> $@.stdout 2> $@.stderr

Output/MutationSharedProportion_vs_CNV/Sequenza/Bowtie2-ADC.tar: Python/compare_MutationSharedProportion_Sequenza.py $(addprefix Output/Sequenza/Bowtie2/,$(addsuffix /sample_alternative_solutions.txt,$(WES_REAL_TUMOR_ID))) Output/MutationSharedProportion/Bowtie2.tsv | build.log Output/MutationSharedProportion_vs_CNV/Sequenza
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --cpus $(CPUS) --ADC 1> $@.stdout 2> $@.stderr

step65: Output/MutationSharedProportion_vs_CNV/Sequenza/BWA-SQC.tar Output/MutationSharedProportion_vs_CNV/Sequenza/Bowtie2-SQC.tar
.PHONY: step65

# Step 66 (Compare PathSeq WES vs. WTS)
Output/Pathseq/WES_vs_WTS:
	mkdir -p $@

Output/Pathseq/WES_vs_WTS/venn.%.pdf: Python/compare_PathSeq_venn.py Output/Pathseq/WES/BWA.%.tsv Output/Pathseq/WES/Bowtie2.%.tsv Output/Pathseq/WTS/STAR.%.tsv Output/Pathseq/WTS/Bowtie2.%.tsv | Output/Pathseq/WES_vs_WTS build.log
	$(DOCKER) $(PYTHON) $(addprefix /,$^ $@) --annotation "WES-BWA" "WES-Bowtie2" "WTS-STAR" "WTS-Bowtie2" --SQC 1> $@.stdout 2> $@.stderr

step66-venn: $(addprefix Output/Pathseq/WES_vs_WTS/venn.,$(addsuffix .pdf,$(TAXONOMY_TYPES)))
step66: step66-venn
.PHONY: step66-venn step66

# Step ## (QAPA)
# https://github.com/morrislab/qapa

# Step ## (Roar)
# http://bioconductor.org/packages/release/bioc/html/roar.html

# Step ## (MesKit)
# https://www.bioconductor.org/packages/release/bioc/html/MesKit.html

# Step ## (Spruce)
# https://github.com/raphael-group/spruce/

# Step ## (CRUST)
