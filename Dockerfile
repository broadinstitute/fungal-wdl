FROM openjdk:8-jre
MAINTAINER Aina Martinez Zurita <amartine@broadinstitute.org>

ENV TERM=xterm-256color

LABEL SAMTOOLS_VER=1.10
LABEL BWA_VER=0.7.12
LABEL BWA_VER=4.1.4.1


#Install Basic Utilities, Python and R

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    libbz2-dev \
    libcurl3-dev \
    libgsl-dev \
    liblzma-dev \
    libncurses5-dev \
    python \
    python3-pip \
    r-base \
    unzip \
    wget \
    zlib1g-dev

#Install R dependencies
RUN echo 'install.packages(c("ggplot2"), repos="http://cran.us.r-project.org", dependencies=TRUE)' > /tmp/packages.R && \
     Rscript /tmp/packages.R


#Install samtools
RUN cd /opt && \
    wget -O samtools-1.10.tar.bz2 'https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2'  && \
    tar xf samtools-1.10.tar.bz2 && \
    rm samtools-1.10.tar.bz2 && \
    cd samtools-1.10 && \
    make && \
    make install && \
    make clean

#Install BWA
RUN cd /opt && \
    wget 'https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2' -O bwa-0.7.17.tar.bz2 && \
    tar xf bwa-0.7.17.tar.bz2 && \
    rm bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make
ENV PATH /opt/bwa-0.7.17:$PATH

#Install GATK4
RUN cd /opt && \
    wget 'https://github.com/broadinstitute/gatk/releases/download/4.1.4.1/gatk-4.1.4.1.zip' -O gatk-4.1.4.1.zip && \
    unzip gatk-4.1.4.1.zip && \
    rm gatk-4.1.4.1.zip

ENV PATH /opt/gatk-4.1.4.1:$PATH
