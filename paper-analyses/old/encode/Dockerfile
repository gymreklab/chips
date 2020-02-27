FROM gymreklab/chipmunk-1.9

RUN apt-get update

RUN apt-get install -y \
    curl \
    python3 \
    python3-pip \
    unzip

RUN pip3 install awscli

# Add the fetch_and_run.sh script
ADD fetch_and_run.sh /usr/local/bin/fetch_and_run.sh

# Get set up to run
WORKDIR /tmp
ENTRYPOINT ["/usr/local/bin/fetch_and_run.sh"]