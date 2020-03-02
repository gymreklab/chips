FROM ubuntu:16.04

# Update necessary packages
RUN apt-get update && apt-get install -qqy \
    autoconf \
    autoconf-archive \
    build-essential \
    default-jre \
    git \
    libtool libtool-bin \
    libbz2-dev \
    liblzma-dev \
    make \
    pkg-config \
    wget \
    zlib1g-dev

# Install htslib
RUN mkdir /dependencies
WORKDIR /dependencies
RUN wget -O htslib-1.8.tar.bz2 https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2
RUN tar -xjvf htslib-1.8.tar.bz2
WORKDIR htslib-1.8
RUN ./configure --prefix=/usr/local --disable-lzma --disable-bz2 && make && make install

# Add samtools
WORKDIR /dependencies
RUN wget -O samtools-1.9.tar.bz2 https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
RUN tar -xjf samtools-1.9.tar.bz2
WORKDIR samtools-1.9
RUN ./configure --without-curses && make && make install

# Install Chips
WORKDIR /dependencies
COPY chips-1.10.1.tar.gz /dependencies/chips-1.10.1.tar.gz
RUN tar -xzvf chips-1.10.1.tar.gz
WORKDIR /dependencies/chips-1.10.1
RUN ./configure && make && make install
ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/local/lib"

# Install Picard
RUN wget -O /dependencies/picard.jar https://github.com/broadinstitute/picard/releases/download/2.18.27/picard.jar
ENV PICARD=/dependencies/picard.jar

