FROM ubuntu:18.04

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install wget zip gcc zlib1g-dev python3 -y && \
    rm -r /var/lib/apt/lists/*