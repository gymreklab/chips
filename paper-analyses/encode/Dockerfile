FROM gymreklab/chipmunk-1.6

RUN apt-get update

RUN apt-get install -y \
    default-jre \
    python3 \
    python3-pip \
    unzip

RUN pip3 install awscli

# Add samtools
WORKDIR /dependencies
RUN wget -O samtools-1.9.tar.bz2 https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
RUN tar -xjf samtools-1.9.tar.bz2
WORKDIR samtools-1.9
RUN ./configure --without-curses && make && make install

# Add the fetch_and_run.sh script
ADD fetch_and_run.sh /usr/local/bin/fetch_and_run.sh

# Get set up to run
WORKDIR /tmp
ENTRYPOINT ["/usr/local/bin/fetch_and_run.sh"]