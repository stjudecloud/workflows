FROM ubuntu:22.04

RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install wget zip gcc zlib1g-dev libbz2-dev liblzma-dev jq python3 python3-pip parallel -y \
    && rm -r /var/lib/apt/lists/* \
    && python3 -m pip --no-cache-dir install csvkit pysam
