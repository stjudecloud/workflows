# Supply a valid download link and md5sum from "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest"

FROM ubuntu:20.04

ARG CELLRANGER_URL
ARG CELLRANGER_MD5

RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install curl -y \
    && rm -r /var/lib/apt/lists/*

WORKDIR /opt

RUN curl -o cellranger.tar.gz \
    ${CELLRANGER_URL} \
    && echo "${CELLRANGER_MD5}  cellranger.tar.gz" > cellranger.tar.gz.md5 \
    && md5sum -c cellranger.tar.gz.md5 \
    && tar -xzvf cellranger.tar.gz \
    && mv cellranger-* cellranger \
    && rm cellranger.tar.gz*

ENV PATH "/opt/cellranger:$PATH"

ENTRYPOINT [ "cellranger" ]
