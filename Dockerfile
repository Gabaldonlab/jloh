FROM python:3.9.1
#Image in buster flavor
# metadata
LABEL base.image="python:3.9.1"
LABEL version="0.21.0"
LABEL software="JLOH"
LABEL software.version="Latest"
LABEL description="A tool to extract, filter, and manage blocks of loss of heterozygosity (LOH) based on single-nucleotide polymorphisms (SNPs), read mapping, and a reference genome"
LABEL website="https://github.com/Gabaldonlab/jloh"
LABEL license="GNU General Public License 3.0"
LABEL maintainer="Matteo Schivinatto (BSC): JLOH developer, Diego Fuentes (BSC): Docker image developer"
##Set up bash and install basic dependencies
SHELL ["/bin/bash", "-c"]
RUN apt-get update -qq && apt-get install -y perl default-jre python3-pip git make nano automake wget g++ zlib1g-dev libbz2-dev libncurses5-dev libncursesw5-dev liblzma-dev curl mummer=3.23+dfsg-4
##Install from source htslib, samtools and mummer
RUN cd /usr/bin && wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && tar -vxjf htslib-1.16.tar.bz2 && rm htslib-1.16.tar.bz2 && cd htslib-1.16 && make install
RUN cd /usr/bin && wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && tar -vxjf samtools-1.13.tar.bz2  && rm samtools-1.13.tar.bz2 && cd samtools-1.13 && make install
RUN cd /usr/bin && wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz && tar -zxvf bedtools-2.30.0.tar.gz  && rm bedtools-2.30.0.tar.gz && cd bedtools2 && make install
RUN cd /usr/bin && wget https://github.com/nextflow-io/nextflow/releases/download/v21.04.1/nextflow-21.04.1-all && chmod -R 777 nextflow-21.04.1-all && mv nextflow-21.04.1-all nextflow
##Install and upgrade python dependencies with pip
RUN python3 -m pip install --upgrade pip && pip3 install --upgrade numpy pysam biopython pandas pandarallel pybedtools
##Git clone JLOH, hisat2 and all2vcf
RUN mkdir -p /root/src && cd /root/src && git clone https://github.com/MatteoSchiavinato/all2vcf && ln -s /root/src/all2vcf/all2vcf /usr/bin
RUN cd /root/src && git clone https://github.com/DaehwanKimLab/hisat2.git && cd hisat2 && make && ln -s /root/src/hisat2/hisat2-* /usr/bin && ln -s /root/src/hisat2/hisat2 /usr/bin
RUN cd /root/src && git clone https://github.com/Gabaldonlab/jloh.git && ln -s /root/src/jloh/jloh /usr/bin
#Purge unnecessary dependencies
RUN apt purge -y git make g++ zlib1g-dev python3-pip automake wget curl make zlib1g-dev libbz2-dev libncurses5-dev libncursesw5-dev liblzma-dev
#Test that nextflow runs as proof of concept
RUN sed -i 's/\${projectDir}/\/root\/src\/jloh\/workflows\/simulate_test_dataset/g' /root/src/jloh/workflows/simulate_test_dataset/nextflow.config
RUN cat /root/src/jloh/workflows/simulate_test_dataset/nextflow.config
RUN echo "Running the simulate_test_dataset.nf for 10 seconds!"
RUN nextflow -C /root/src/jloh/workflows/simulate_test_dataset/nextflow.config run /root/src/jloh/workflows/simulate_test_dataset/simulate_test_dataset.nf -work-dir /root/src/jloh/work & sleep 10
