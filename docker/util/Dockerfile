FROM ubuntu:24.04

RUN apt-get update \
    && apt-get upgrade -y \
    &&  apt-get install wget zip gcc zlib1g-dev libbz2-dev \
      liblzma-dev jq python3 python3-pip parallel \
      python3-csvkit python3-pysam -y \
    && rm -r /var/lib/apt/lists/*

COPY --from=scripts --chmod=777 util/calc_global_phred_scores.py /scripts/util/calc_global_phred_scores.py
