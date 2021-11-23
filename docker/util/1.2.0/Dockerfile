FROM ubuntu:18.04

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install wget zip gcc zlib1g-dev jq python3 python3-pip -y && \
    rm -r /var/lib/apt/lists/* && \
    python3 -m pip --no-cache-dir install csvkit
