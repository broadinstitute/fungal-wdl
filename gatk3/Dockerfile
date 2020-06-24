FROM java:openjdk-8-jre
MAINTAINER Xiao Li <xiaoli@broadinstitute.org>

ENV TERM=xterm-256color
ENV DOCKER_FIX=''

LABEL GOTC_PICARD_VER=1.782
LABEL GOTC_GATK37_VER=3.7-93-ge9d8068
LABEL GOTC_SAMTOOLS_VER=1.3.1
LABEL GOTC_BWA_VER=0.7.12
LABEL GOTC_TABIX_VER=0.2.5_r1005
LABEL GOTC_BGZIP_VER=1.3

RUN apt-get -qq update && apt-get install -qqy \
    build-essential \
    curl \
    apt-utils \
    libbz2-dev \
    libcurl3-dev \
    libgsl-dev \
    liblzma-dev \
    libncurses5-dev \
    unzip \
    vim-common \
    wget \
    zlib1g-dev \
    python \
    python3-pip \
    r-base

# Install ggplot2
RUN echo 'install.packages(c("ggplot2"), repos="http://cran.us.r-project.org", dependencies=TRUE)' > /tmp/packages.R && \
     Rscript /tmp/packages.R

# Install samtools
RUN cd /opt && \
    wget 'https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2' -O samtools-1.3.1.tar.bz && \
    tar xf samtools-1.3.1.tar.bz && \
    rm samtools-1.3.1.tar.bz && \
    cd samtools-1.3.1 && \
    make && \
    make install && \
    make clean

# bwa
RUN cd /opt && \
    wget 'https://github.com/lh3/bwa/archive/0.7.12.tar.gz' -O bwa-0.7.12.tar.gz && \
    tar xf bwa-0.7.12.tar.gz && \
    rm bwa-0.7.12.tar.gz && \
    cd bwa-0.7.12 && \
    make
ENV PATH /opt/bwa-0.7.12:$PATH

# Copy jars to container
COPY bin/* /opt/

# CLEAN UP
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/
