FROM ubuntu:16.04
RUN apt-get update \
  && apt-get install -y hmmer \
  && apt-get install -y mafft \
  && apt-get autoremove \
  && apt-get clean

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev \
  && cd /usr/local/bin \
  && pip3 install --upgrade pip \
  && apt-get autoremove \
  && apt-get clean

RUN pip3 install 'snakemake==5.5.4'
RUN pip3 install 'pandas==0.24.2'
RUN pip3 install 'matplotlib==3.0.3'
RUN pip3 install 'scikit-learn==0.21.3'

ADD scripts /usr/local/bin
RUN chmod +x /usr/local/bin/*
