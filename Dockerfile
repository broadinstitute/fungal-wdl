FROM openjdk:8-jre
MAINTAINER Aina Martinez Zurita <amartine@broadinstitute.org>

ENV TERM=xterm-256color

LABEL SAMTOOLS_VER=1.10
LABEL HTSLIB_VER=1.10
LABEL BWA_VER=0.7.12


#Install Basic Utilities, Python and R

RUN apt-get update && apt-get install -y \
    apt-utils \
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
    vim-common \
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

#Install GATK4

#Clean
