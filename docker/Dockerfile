FROM ubuntu:jammy

RUN apt-get update && \
    apt-get install -y \
    build-essential \
    wget \
    git \
    python3 \
    python3-pip

# Install all software under /opt/software:
RUN mkdir -p /opt/software 
ADD dragon.py /opt/software/ 
FROM ubuntu:jammy

RUN apt-get update && \
    apt-get install -y \
    build-essential \
    wget \
    git \
    python3 \
    python3-pip

# Install all software under /opt/software:
RUN mkdir -p /opt/software && \
    mkdir -p /opt/software/resources 

ADD requirements.txt /opt/software/
RUN pip3 install -r /opt/software/requirements.txt

ADD preprocess.py /opt/software/
ADD reformatted_probe_mapping.hm450.tsv /opt/software/resources/
