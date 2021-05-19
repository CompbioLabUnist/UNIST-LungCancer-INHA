FROM ubuntu:20.04
LABEL maintainer="jwlee230@unist.ac.kr"

RUN ln -sf /usr/share/zoneinfo/Asia/Seoul /etc/localtime && apt-get update && apt-get upgrade -y && apt-get install -y make gcc g++ zlib1g-dev libbz2-dev liblzma-dev autoconf libncurses5-dev unzip openjdk-8-jdk r-base r-base-dev libfontconfig1 libssl-dev openssl libffi-dev tabix libcurl4-openssl-dev curl perl cpanminus git libpng-dev libmysqlclient-dev libxml2-dev cmake libboost-all-dev

# TOOL directory
RUN mkdir /Tools
ADD Tools/* /Tools/

# Python3
RUN cd /Tools/Python-3.9.1 && ./configure && make -j && make -j install && ln -sf /usr/local/bin/python3 /usr/bin/python
ADD requirements.txt /Tools/requirements.txt
RUN pip3 install --no-cache-dir --upgrade wheel pip && pip3 install --requirement /Tools/requirements.txt --no-cache-dir

# Conda
RUN cd /Tools && bash Miniconda.sh -b
ENV PATH=/root/miniconda/bin:$PATH
RUN conda install pyclone -c bioconda -c conda-forge

# Samtools
RUN cd /Tools/samtools-1.11 && autoheader && autoconf -Wno-syntax && ./configure && make -j && make -j install
ENV PATH=/Tools/samtools-1.10:$PATH

# GATK
RUN cd Tools && unzip gatk.zip
ENV PATH=/Tools/gatk-4.1.9.0:$PATH

# BWA
RUN cd /Tools/bwa-0.7.17 && make -j
ENV PATH=/Tools/bwa-0.7.17:$PATH

# FastQC
RUN cd /Tools && unzip FastQC.zip && cd FastQC && chmod +x fastqc
ENV PATH=/Tools/FastQC:$PATH

# STAR
ENV PATH=/Tools/STAR-2.7.7a/bin/Linux_x86_64:$PATH

# Bowtie2
RUN cd /Tools && unzip Bowtie2.zip
ENV PATH=/Tools/bowtie2-2.4.2-linux-x86_64:$PATH BOWTIE2_INDEXES=/Data/hg38

# R packages
RUN mkdir -p $HOME/.R && echo "MAKEFLAGS = -j" > $HOME/.R/Makevars && Rscript -e "install.packages(c('optparse', 'BiocManager', 'data.table'))" && Rscript -e "BiocManager::install(c('copynumber', 'DESeq2', 'MesKit'))" && Rscript -e "install.packages(c('sequenza'))"

# VEP
ENV KENT_SRC=/Tools/kent-335_base/src MACHTYPE=x86_64 CFLAGS="-fPIC"
RUN cd /Tools/kent-335_base/src/lib && make -j && cd ../jkOwnLib && make -j && cd /Tools && unzip VEP.zip && cd ensembl-vep-release-103 && cpanm --with-recommends --notest --installdeps . && perl INSTALL.pl --NO_TEST --NO_UPDATE && cd /Tools/ensembl-vep-release-103/Bio/EnsEMBL/Variation && sed -i -e '900,906d' TranscriptVariationAllele.pm
ENV PATH=/Tools/ensembl-vep-release-103:$PATH

# BCFtools
RUN cd /Tools/bcftools-1.12 && ./configure && make -j && make -j install

# MutEnricher
RUN cd /Tools/MutEnricher-1.3.1/math_funcs && python3 setup.py build_ext --inplace

# RSEM
RUN cd /Tools/RSEM-1.3.3 && make -j && make -j ebseq && make -j install

# Lemon
RUN cd /Tools/lemon-1.3 && mkdir build && cd build && cmake .. && make -j && make -j install

# Spruce
RUN cd /Tools && git clone "https://github.com/raphael-group/spruce.git" && cd spruce && mkdir build && cd build && cmake .. && make -j
ENV PATH=/Tools/spruce/build:$PATH

VOLUME /Output
VOLUME /Data
VOLUME /Python

CMD ["/bin/bash"]
