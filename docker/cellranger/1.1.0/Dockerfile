# Supply a valid download link and md5sum from "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest"

FROM ubuntu:20.04

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install curl -y && \
    rm -r /var/lib/apt/lists/*

WORKDIR /opt

ARG cellranger_url
ARG cellranger_md5

RUN curl -o cellranger.tar.gz \
    ${cellranger_url} \
    && echo "${cellranger_md5}  cellranger.tar.gz" > cellranger.tar.gz.md5 \
    && md5sum -c cellranger.tar.gz.md5 \
    && tar -xzvf cellranger.tar.gz \
    && mv cellranger-* cellranger

ENV PATH "/opt/cellranger:$PATH"
