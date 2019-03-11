FROM ubuntu:16.04

# Update necessary packages
RUN apt-get update && apt-get install -qqy \
    autoconf \
    autoconf-archive \
    build-essential \
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

# Install ChIPmunk
WORKDIR /dependencies
RUN git clone https://github.com/gymreklab/ChIPmunk
WORKDIR ChIPmunk
RUN ./reconf && ./configure --prefix=/usr/local && make && make install
RUN ldconfig