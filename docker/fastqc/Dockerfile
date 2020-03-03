FROM stjudecloud/conda-base:1.0.0-alpha AS builder

RUN conda create -n fastqc \
        fastqc==0.11.8 \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/fastqc/bin:$PATH

ENTRYPOINT [ "fastqc" ]