FROM stjudecloud/conda-base:1.1.1 AS builder

RUN conda create -n bwa \
    bwa==0.7.17 \
    samtools==1.11 \
    -y \
    && conda clean --all

ENV PATH /opt/conda/envs/bwa/bin:$PATH

ENTRYPOINT [ "bwa" ]
