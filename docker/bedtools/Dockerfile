FROM ubuntu:22.04 AS builder

RUN apt-get update && apt-get install -y build-essential zlib1g-dev libbz2-dev liblzma-dev wget python3

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz \
    && tar -xzf bedtools-2.31.1.tar.gz \
    && cd bedtools2 \
    && make static \
    && make all \
    && mv bin/bedtools.static bin/bedtools \
    && cp bin/* /usr/local/bin/ \
    && cd .. \
    && rm -rf bedtools-2.31.1.tar.gz bedtools2

FROM python:3.13.1-alpine

COPY --from=builder /usr/local/bin/ /usr/local/bin/
COPY --from=scripts --chmod=777 hic/filter_hic.py /usr/local/bin/filter_hic.py

RUN apk add --no-cache bc bash