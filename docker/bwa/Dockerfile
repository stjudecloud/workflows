FROM stjudecloud/conda-base:1.0.0 AS builder

RUN conda create -n bwa \
        bwa==0.7.17 \
        -y \
    && conda clean --all

ENV PATH /opt/conda/envs/bwa/bin:$PATH

ENTRYPOINT [ "bwa" ]
