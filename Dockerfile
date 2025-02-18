FROM ubuntu:22.04

LABEL maintainer="cschu1981@gmail.com"
LABEL version="2.18.3"
LABEL description="gffquant - functional profiling of metagenomic/transcriptomic wgs samples"


ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get install -y wget python3-pip git gawk bwa minimap2 samtools
# RUN apt-get upgrade -y

# RUN apt install -y wget python3-pip git dirmngr gnupg ca-certificates build-essential libssl-dev libcurl4-gnutls-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev bwa minimap2
# RUN apt clean

ADD LICENSE DESCRIPTION.md README.md environment.yml setup.py /opt/software/gffquant/
ADD gffquant /opt/software/gffquant/gffquant/

RUN cd /opt/software/gffquant && pip install .


# ARG CACHEBUST=1
# RUN mkdir -p /opt/software && \
# 	cd /opt/software && \
# 	git clone https://github.com/cschu/gff_quantifier && \
# 	cd gff_quantifier && \
# 	pip install .
  
CMD ["gffquant"]
