FROM gymreklab/chipmunk-1.9

RUN apt-get install -y \
    curl \
    python \
    python-dev \
    python-pip \
    unzip

RUN pip install awscli

# Add power analysis script
ADD power-analysis.py /usr/local/bin/power-analysis.py
RUN chmod +x /usr/local/bin/power-analysis.py

# Install BWA
WORKDIR /dependencies
RUN git clone https://github.com/lh3/bwa.git
WORKDIR bwa
RUN make
RUN cp bwa /usr/local/bin/bwa

# Install MACS
RUN pip install numpy
RUN pip install MACS2

# Add the fetch_and_run.sh script
ADD fetch_and_run.sh /usr/local/bin/fetch_and_run.sh

# Get set up to run
WORKDIR /tmp
ENTRYPOINT ["/usr/local/bin/fetch_and_run.sh"]
