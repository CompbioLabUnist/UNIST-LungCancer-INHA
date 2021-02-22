FROM ubuntu:20.04
LABEL maintainer="jwlee230@unist.ac.kr"

RUN ln -sf /usr/share/zoneinfo/Asia/Seoul /etc/localtime

RUN apt-get update && apt-get upgrade -y && apt-get install -y make gcc g++ zlib1g-dev libbz2-dev liblzma-dev autoconf libncurses5-dev unzip openjdk-8-jdk r-base r-base-dev libfontconfig1 libssl-dev openssl libffi-dev tabix libcurl4-openssl-dev curl perl

# TOOL directory
RUN mkdir /Tools

# Python3
ADD Tools/Python.tar.xz /Tools
RUN cd /Tools/Python-3.9.1 && ./configure && make -j && make -j install
RUN ln -sf /usr/local/bin/python3 /usr/bin/python

# PIP
ADD requirements.txt /Tools/requirements.txt
RUN pip3 install --requirement /Tools/requirements.txt

# Samtools
ADD Tools/samtools.tar.bz2 /Tools
RUN cd /Tools/samtools-1.11 && autoheader && autoconf -Wno-syntax && ./configure && make -j && make -j install
ENV PATH=/Tools/samtools-1.10:$PATH

# GATK
ADD Tools/gatk.zip /Tools
RUN cd Tools && unzip gatk.zip
ENV PATH=/Tools/gatk-4.1.9.0:$PATH

# BWA
ADD Tools/bwa.tar.bz2 /Tools
RUN cd /Tools/bwa-0.7.17 && make -j
ENV PATH=/Tools/bwa-0.7.17:$PATH

# Picard
ADD Tools/picard.jar /Tools
ENV PICARD=/Tools/picard.jar

# FastQC
ADD Tools/FastQC.zip /Tools
RUN cd /Tools && unzip FastQC.zip && cd FastQC && chmod +x fastqc
ENV PATH=/Tools/FastQC:$PATH

# STAR
ADD Tools/STAR.tar.gz /Tools
ENV PATH=/Tools/STAR-2.7.7a/bin/Linux_x86_64:$PATH

# Bowtie2
ADD Tools/Bowtie2.zip /Tools
RUN cd /Tools && unzip Bowtie2.zip
ENV PATH=/Tools/bowtie2-2.4.2-linux-x86_64:$PATH

# VCF2MAF
ADD Tools/vcf2maf.tar.gz /Tools

# R packages
RUN Rscript -e "install.packages(c('optparse', 'BiocManager'))" && Rscript -e "BiocManager::install('copynumber')" && Rscript -e "install.packages('sequenza')"

VOLUME /Output
VOLUME /Data
VOLUME /Python

CMD ["/bin/bash"]
