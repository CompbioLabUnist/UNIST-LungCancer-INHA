FROM mathworks/matlab:r2021b
USER root
WORKDIR /
LABEL maintainer="jwlee230@unist.ac.kr"

# apt install & update
RUN ln -sf /usr/share/zoneinfo/Asia/Seoul /etc/localtime
RUN apt-get update && apt-get upgrade -y && apt-get install -y software-properties-common && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 089EBE08314DF160 && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 089EBE08314DF160 51716619E084DAB9 && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && apt-get update && apt-get upgrade -y && apt-get install -y make gcc g++ zlib1g-dev libbz2-dev liblzma-dev autoconf libncurses5-dev unzip openjdk-8-jdk r-base r-base-core r-recommended r-base-dev libfontconfig1 libssl-dev openssl libffi-dev tabix libcurl4-openssl-dev curl perl cpanminus git libpng-dev libmysqlclient-dev libxml2-dev cmake libboost-all-dev gcc-7 g++-7 gfortran-7 wget libdb-dev rsync git-lfs libncurses5 libgeos-dev libsqlite3-dev libharfbuzz-dev libfribidi-dev libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libcairo2-dev && mkdir -p /root/.local/share/fonts && wget 'https://raw.githubusercontent.com/matomo-org/travis-scripts/master/fonts/Arial.ttf' -O /root/.local/share/fonts/Arial.ttf

# TOOL directory
RUN mkdir /Tools

# Python3
ADD Tools/Python.tar.xz /Tools/
RUN cd /Tools/Python-3.9.1 && ./configure --enable-loadable-sqlite-extensions --enable-optimizations --enable-big-digits=30 --with-ensurepip && make -j && make -j install && make -j altinstall
ADD requirements.txt /Tools/requirements.txt
RUN pip3 install wheel==0.40.0 pip==23.1.2 && pip3 install --requirement /Tools/requirements.txt

# Conda
ADD Tools/Miniconda.sh /Tools/
RUN cd /Tools && bash Miniconda.sh -b
ENV PATH=/root/anaconda2/bin:$PATH
RUN conda create -n pyclone -c bioconda -c conda-forge pyclone==0.13.1

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
RUN mkdir -p $HOME/.R && echo "MAKEFLAGS = -j" > $HOME/.R/Makevars && Rscript -e "install.packages(c('optparse', 'BiocManager', 'data.table', 'e1071', 'colorRamps', 'devtools', 'R.utils'))" && Rscript -e "install.packages('remotes')" && Rscript -e "remotes::install_github('renozao/xbioc@1354168bd7e64be4ee1f9f74e971a61556d75003')" && Rscript -e "BiocManager::install(c('copynumber', 'DESeq2', 'MesKit', 'preprocessCore', 'Biobase'))" && Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/sequenza/sequenza_3.0.0.tar.gz', repos = NULL, type='source')" && Rscript -e "install.packages('Rserve', repos='http://www.rforge.net/')" && Rscript -e "devtools::install_github(c('xuranw/MuSiC@fc0ca72792dfc65003e10162a3619c02fd8742d1', 'meichendong/SCDC@2b0b2abcf5d06361003961bdec4f56223ce7bde2'))" && Rscript -e "install.packages(c('igraph', 'tidygraph', 'ggraph'))" && Rscript -e "devtools::install_github(c('caravagn/ctree@df8dc8395725f396e9f7cefa3a9e46d3b7b113ec', 'caravagn/mtree@8459f273264f425e60d8694258d94465c05217dd', 'caravagnalab/revolver@7ee825ae419d38ceb2b1613bb19cb409952ab00b'))" && Rscript -e "install.packages('remotes')" && Rscript -e "remotes::install_github('ShixiangWang/copynumber@b3de1b00165b112bdcbc44be61e236a88a634c9a')" && Rscript -e "remotes::install_github('Subhayan18/CRUST@105a60f211b411f9d00a6ab8e17c947f4a521227')" && Rscript -e "BiocManager::install(c('roar', 'EBSeq', 'NOISeq', 'PureCN', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38'))" && Rscript -e "install.packages(c('BisqueRNA', 'Seurat', 'deconstructSigs'))" && Rscript -e "BiocManager::install(c('affy', 'oligo', 'pd.clariom.d.human', 'affycoretools', 'clariomdhumantranscriptcluster.db'))"
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
RUN cd /Tools && unzip QAPA.zip && cd qapa-1.3.3 && conda install -c bioconda -c conda-forge mamba && mamba env create -f environment.yml

# SigProfiler
RUN python3 -c "from SigProfilerMatrixGenerator import install as genInstall; genInstall.install('GRCh38')"

# Mobster
ADD Tools/mobster.tar.gz /Tools/
RUN cd /Tools && git clone https://github.com/jyhehir/mobster && cd mobster && git checkout 5b7bae9b7eab057cc969c07c895c17401072eed0 && git-lfs pull && conda install -c bioconda -c conda-forge mobster

# Gistic
ADD Tools/Gistic.tar.gz /Tools/
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6
RUN cp /Tools/MCR_Installer/MCRInstaller.zip /usr && cd /usr && unzip MCRInstaller.zip && ./install -mode silent -agreeToLicense yes -destinationFolder /Tools/MATLAB_Compiler_Runtime && echo -e "/opt/matlab/R2021b/bin/glnxa64\n/opt/matlab/R2021b/sys/os/glnxa64\n/Tools/MATLAB_Compiler_Runtime/v83/runtime/glnxa64\n/Tools/MATLAB_Compiler_Runtime/v83/bin/glnxa64\n/Tools/MATLAB_Compiler_Runtime/v83/sys/os/glnxa64" > /etc/ld.so.conf.d/matlab.conf && ldconfig && unlink /Tools/MATLAB_Compiler_Runtime/v83/bin/glnxa64/libcurl.so.4 && ln -s /Tools/gistic2 /usr/bin

# MiXCR
ADD Tools/MiXCR.zip /Tools/
RUN cd /Tools && unzip MiXCR.zip
ENV PATH=/Tools/mixcr-3.0.13:$PATH

# Bedops
ADD Tools/Bedops.tar.bz2 /Tools/
ENV PATH=/Tools/bin:$PATH

# iREAD
ADD Tools/iRead.zip /Tools/
RUN cd /Tools && unzip iRead.zip && cd iread-0.8.9 && chmod +x * && ln -s /root/anaconda2/bin/python2 /usr/bin/python2 && cpanm --with-recommends --notest Parallel::ForkManager
ENV PATH=/Tools/iread-0.8.9:$PATH

# BeadArrayFiles
ADD Tools/BeadArrayFiles.tar.gz /Tools
RUN cd /Tools/BeadArrayFiles-1.3.4 && python3 setup.py install

# Musical
RUN cd /Tools & git clone https://github.com/parklab/MuSiCal && cd MuSiCal && git checkout 85ae0a1a9cbd5daed7800f49b138b9bffa5070a2 && pip3 install .

# Purge
RUN cd /Tools && rm -v *.zip

RUN cd /Tools
VOLUME /Output

CMD ["/bin/bash"]
