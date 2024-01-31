FROM python:3.9.1
#Image in buster flavor
# metadata
LABEL base.image="python:3.9.1"
LABEL version="1.0.3"
LABEL software="JLOH"
LABEL description="A tool to extract, filter, and manage blocks of loss of heterozygosity (LOH) based on single-nucleotide polymorphisms (SNPs), read mapping, and a reference genome"
LABEL website="https://github.com/Gabaldonlab/jloh"
LABEL license="GNU General Public License 3.0"
LABEL maintainer="Matteo Schiavinato (BSC): JLOH developer, Diego Fuentes (BSC): Docker image developer"
##Set up bash and install basic dependencies
SHELL ["/bin/bash", "-c"]
RUN apt-get update -qq \
    && apt-get install -y perl default-jre python3-pip git make nano automake wget g++ zlib1g-dev libbz2-dev libncurses5-dev libncursesw5-dev liblzma-dev curl mummer=3.23+dfsg-4

##Install from source htslib, samtools and mummer
RUN cd /usr/bin \
    && wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 \
    && tar -vxjf htslib-1.16.tar.bz2 \
    && rm htslib-1.16.tar.bz2 \
    && cd htslib-1.16 \
    && make install

RUN cd /usr/bin \
    && wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 \
    && tar -vxjf samtools-1.13.tar.bz2  \
    && rm samtools-1.13.tar.bz2 \
    && cd samtools-1.13 \
    && make install

RUN cd /usr/bin \
    && wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz \
    && tar -zxvf bedtools-2.30.0.tar.gz \
    && rm bedtools-2.30.0.tar.gz \
    && cd bedtools2 \
    && make install

RUN cd /usr/bin \
    && wget https://github.com/nextflow-io/nextflow/releases/download/v21.04.1/nextflow-21.04.1-all \
    && chmod -R 777 nextflow-21.04.1-all \
    && mv nextflow-21.04.1-all nextflow

##Install and upgrade python dependencies with pip
RUN python3 -m pip install --upgrade pip \
    && pip3 install --upgrade numpy pysam biopython pandas pandarallel pybedtools

##Git clone hisat2 and all2vcf
RUN mkdir -p /root/src \
    && cd /root/src \
    && git clone -b v0.7.8 https://github.com/MatteoSchiavinato/all2vcf \
    && ln -s /root/src/all2vcf/all2vcf /usr/bin

RUN cd /root/src \
    && git clone https://github.com/DaehwanKimLab/hisat2.git \
    && cd hisat2 \
    && make \
    && ln -s /root/src/hisat2/hisat2-* /usr/bin \
    && ln -s /root/src/hisat2/hisat2 /usr/bin

##Download R from source, compile and set up dependencies. Tends to be a long process
RUN apt-get upgrade -y gfortran libreadline6-dev libx11-dev libxt-dev libpng-dev libjpeg-dev libcairo2-dev xvfb libzstd-dev libcurl4-openssl-dev texinfo texlive texlive-fonts-extra screen libpcre2-dev \
    && cd /usr/local/src \
    && wget https://cran.r-project.org/src/base/R-3/R-3.6.0.tar.gz \
    && tar zxvf R-3.6.0.tar.gz && rm R-3.6.0.tar.gz && cd R-3.6.0/ \
    && ./configure --enable-R-shlib \
    && make \
    && make install 

RUN Rscript -e 'install.packages(c("ggplot2", "reshape2", "hash", "png"), repos = "http://cran.us.r-project.org")'

# clone JLOH 
RUN cd /root/src \
    && git clone -b v1.0.3 https://github.com/Gabaldonlab/jloh.git \
    && ln -s /root/src/jloh/jloh /usr/bin

#Purge unnecessary dependencies
RUN apt purge -y git make g++ zlib1g-dev python3-pip automake wget curl make zlib1g-dev libbz2-dev libncurses5-dev libncursesw5-dev liblzma-dev libzstd-dev libreadline6-dev libxt-dev \
    && rm -rf /var/lib/apt/lists/*

#Check R version
RUN R --version

#Running the tests"
RUN jloh stats --vcf /root/src/jloh/test_data/out.ff.vcf
RUN jloh extract --vcf /root/src/jloh/test_data/out.ff.vcf --bam /root/src/jloh/test_data/out.fs.bam --ref /root/src/jloh/test_data/S_para.chrXII.fa --min-snps-kbp 2,5 --output-dir /root/src/jloh/test_data/jloh_out
RUN jloh plot --one-ref --loh /root/src/jloh/test_data/jloh_out/jloh.LOH_blocks.tsv --het /root/src/jloh/test_data/jloh_out/jloh.exp.het_blocks.bed --contrast max
