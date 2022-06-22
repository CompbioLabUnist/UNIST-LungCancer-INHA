FROM mathworks/matlab:r2021b
USER root
WORKDIR /
LABEL maintainer="jwlee230@unist.ac.kr"

# apt install & update
RUN ln -sf /usr/share/zoneinfo/Asia/Seoul /etc/localtime
RUN apt-get update && apt-get upgrade -y && apt-get install -y software-properties-common && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 089EBE08314DF160 && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 089EBE08314DF160 51716619E084DAB9 && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && apt-get update && apt-get upgrade -y && apt-get install -y make gcc g++ zlib1g-dev libbz2-dev liblzma-dev autoconf libncurses5-dev unzip openjdk-8-jdk r-base r-base-core r-recommended r-base-dev libfontconfig1 libssl-dev openssl libffi-dev tabix libcurl4-openssl-dev curl perl cpanminus git libpng-dev libmysqlclient-dev libxml2-dev cmake libboost-all-dev gcc-7 g++-7 gfortran-7 wget libdb-dev rsync git-lfs libncurses5 cabextract libgeos-dev libsqlite3-dev && apt-get install libsqlite3-dev -y && wget "https://www.freedesktop.org/software/fontconfig/webfonts/webfonts.tar.gz" && tar -zxvf webfonts.tar.gz && cd msfonts && cabextract *.exe && mkdir -p ~/.local/share/fonts && cp *.ttf *.TTF ~/.local/share/fonts/

# TOOL directory
RUN mkdir /Tools

# Python3
ADD Tools/Python.tar.xz /Tools/
RUN cd /Tools/Python-3.9.1 && ./configure --enable-loadable-sqlite-extensions --enable-optimizations && make -j && make -j install && make -j altinstall
ADD requirements.txt /Tools/requirements.txt
RUN pip3 install --upgrade wheel pip && pip3 install --requirement /Tools/requirements.txt

# Conda
ADD Tools/Miniconda.sh /Tools/
RUN cd /Tools && bash Miniconda.sh -b
ENV PATH=/root/anaconda2/bin:$PATH
RUN conda create -n pyclone -c bioconda -c conda-forge pyclone

# Samtools
ADD Tools/samtools.tar.bz2 /Tools/
RUN cd /Tools/samtools-1.11 && autoheader && autoconf -Wno-syntax && ./configure && make -j && make -j install
ENV PATH=/Tools/samtools-1.10:$PATH

# GATK
ADD Tools/gatk.zip /Tools/
RUN cd Tools && unzip gatk.zip
ENV PATH=/Tools/gatk-4.1.9.0:$PATH

# PICARD
ADD Tools/picard.jar /Tools/

# VCF2MAF
ADD Tools/vcf2maf.tar.gz /Tools/

# BWA
ADD Tools/bwa.tar.bz2 /Tools/
RUN cd /Tools/bwa-0.7.17 && make -j
ENV PATH=/Tools/bwa-0.7.17:$PATH

# FastQC
ADD Tools/FastQC.zip /Tools/
RUN cd /Tools && unzip FastQC.zip && cd FastQC && chmod +x fastqc
ENV PATH=/Tools/FastQC:$PATH

# STAR
ADD Tools/STAR.tar.gz /Tools/
ENV PATH=/Tools/STAR-2.7.8a/bin/Linux_x86_64:$PATH

# Bowtie2
ADD Tools/Bowtie2.zip /Tools/
RUN cd /Tools && unzip Bowtie2.zip
ENV PATH=/Tools/bowtie2-2.4.2-linux-x86_64:$PATH BOWTIE2_INDEXES=/Data/hg38

# R packages
ADD GitHub.pat /root/.Renviron
RUN mkdir -p $HOME/.R && echo "MAKEFLAGS = -j" > $HOME/.R/Makevars && Rscript -e "install.packages(c('optparse', 'BiocManager', 'data.table', 'e1071', 'colorRamps', 'devtools', 'remotes', 'R.utils'))" && Rscript -e "remotes::install_github('renozao/xbioc')" && Rscript -e "BiocManager::install(c('copynumber', 'DESeq2', 'MesKit', 'preprocessCore', 'Biobase'))" && Rscript -e "install.packages(c('sequenza'))" && Rscript -e "install.packages('Rserve', repos='http://www.rforge.net/')" && Rscript -e "devtools::install_github(c('xuranw/MuSiC', 'meichendong/SCDC'))" && Rscript -e "install.packages(c('igraph', 'tidygraph', 'ggraph'))" && Rscript -e "devtools::install_github(c('caravagn/ctree', 'caravagn/mtree', 'caravagnalab/revolver'))" && Rscript -e "install.packages('remotes')" && Rscript -e "remotes::install_github('ShixiangWang/copynumber')" && Rscript -e "remotes::install_github('Subhayan18/CRUST')" && Rscript -e "BiocManager::install(c('roar', 'EBSeq', 'NOISeq', 'PureCN', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38'))" && Rscript -e "install.packages(c('BisqueRNA', 'Seurat', 'deconstructSigs'))"
ENV R_LIBS=/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library

# VEP
ADD Tools/VEP.zip /Tools/
ADD Tools/BigFile.tar.gz /Tools/
ENV KENT_SRC=/Tools/kent-335_base/src MACHTYPE=x86_64 CFLAGS="-fPIC"
RUN cd /Tools/kent-335_base/src/lib && make -j && cd ../jkOwnLib && make -j && cd /Tools && unzip VEP.zip && cd ensembl-vep-release-103 && cpanm --with-recommends --notest --installdeps . && perl INSTALL.pl --NO_TEST --NO_UPDATE && cd /Tools/ensembl-vep-release-103/Bio/EnsEMBL/Variation && sed -i -e '900,906d' TranscriptVariationAllele.pm
ENV PATH=/Tools/ensembl-vep-release-103:$PATH

# BCFtools
ADD Tools/bcftools.tar.bz2 /Tools/
RUN cd /Tools/bcftools-1.12 && ./configure && make -j && make -j install

# MutEnricher
ADD Tools/MutEnricher.tar.gz /Tools/
RUN cd /Tools/MutEnricher-1.3.1/math_funcs && python3 setup.py build_ext --inplace

# RSEM
ADD Tools/RSEM.tar.gz /Tools/
RUN cd /Tools/RSEM-1.3.3 && make -j && make -j ebseq && make -j install

# Lemon
ADD Tools/Lemon.tar.gz /Tools/
RUN cd /Tools/lemon-1.3 && mkdir build && cd build && cmake .. && make -j && make -j install

# Spruce
ADD Tools/Spruce.zip /Tools/
RUN cd /Tools && unzip Spruce.zip && cd spruce-master && mkdir build && cd build && cmake .. && make -j
ENV PATH=/Tools/spruce-master/build:$PATH

# ANNOVAR
ADD Tools/annovar.latest.tar.gz /Tools/
ENV PATH=/Tools/annovar:$PATH

# Arriba
ADD Tools/arriba.tar.gz /Tools/
RUN cd /Tools/arriba_v2.1.0 && sed -i "30s/zlib.net/zlib.net\/fossils/" Makefile && make -j
ENV PATH=/Tools/arriba_v2.1.0:$PATH

# FusionCatcher
ADD Tools/FusionCatcher.py /Tools/
RUN python2 /Tools/FusionCatcher.py --install-all-tools --yes && pip install biopython==1.76
ENV PATH=/opt/fusioncatcher/bin:$PATH

# STAR-Fusion
ADD Tools/Star-Fusion.tar.gz /Tools/
RUN cpanm --notest DB_File URI::Escape Set::IntervalTree Carp::Assert JSON::XS PerlIO::gzip
ENV PATH=/Tools/STAR-Fusion.v1.10.1:$PATH

# QAPA
ADD Tools/QAPA.zip /Tools/
RUN cd /Tools && unzip QAPA.zip && cd qapa-master && conda install -c bioconda -c conda-forge mamba && mamba env create -f environment.yml

# SigProfiler
RUN python3 -c "from SigProfilerMatrixGenerator import install as genInstall; genInstall.install('GRCh38')"

# Mobster
RUN cd /Tools && git clone "https://github.com/jyhehir/mobster.git" && cd mobster && git-lfs pull
RUN conda install -c bioconda -c conda-forge mobster

# Gistic
ADD Tools/Gistic.tar.gz /Tools/
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6
RUN cp /Tools/MCR_Installer/MCRInstaller.zip /usr && cd /usr && unzip MCRInstaller.zip && ./install -mode silent -agreeToLicense yes -destinationFolder /Tools/MATLAB_Compiler_Runtime && echo -e "/opt/matlab/R2021b/bin/glnxa64\n/opt/matlab/R2021b/sys/os/glnxa64\n/Tools/MATLAB_Compiler_Runtime/v83/runtime/glnxa64\n/Tools/MATLAB_Compiler_Runtime/v83/bin/glnxa64\n/Tools/MATLAB_Compiler_Runtime/v83/sys/os/glnxa64" > /etc/ld.so.conf.d/matlab.conf && ldconfig && unlink /Tools/MATLAB_Compiler_Runtime/v83/bin/glnxa64/libcurl.so.4

# MiXCR
ADD Tools/MiXCR.zip /Tools/
RUN cd /Tools && unzip MiXCR.zip
ENV PATH=/Tools/mixcr-3.0.13:$PATH

# RTCR
ADD Tools/RTCR.zip /Tools/
RUN cd /Tools && unzip RTCR.zip && cd uubram-RTCR-5835928 && python2 setup.py install

# Helmsman
ADD Tools/helmsman.tar.gz /Tools/
RUN cd /Tools/helmsman-1.5.3 && conda env create -n helmsman -f env.yml

VOLUME /Output
VOLUME /Data
VOLUME /Python

CMD ["/bin/bash"]
