FROM continuumio/miniconda3:4.10.3 AS builder

RUN echo 'channels:\n  - conda-forge\n  - bioconda\n' > /root/.condarc \
    && conda config --set channel_priority strict \
    && conda create -n fastq_screen \
    fastq-screen==0.13.0=pl526_1 \
    tree==2.0.0 \
    -y \
    && conda clean --all -y

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/fastq_screen/bin /opt/conda/envs/fastq_screen/bin
COPY --from=builder /opt/conda/envs/fastq_screen/lib/ /opt/conda/envs/fastq_screen/lib/
COPY --from=builder /opt/conda/envs/fastq_screen/share/fastq-screen-0.13.0-1/ /opt/conda/envs/fastq_screen/share/fastq-screen-0.13.0-1/
ENV PATH /opt/conda/envs/fastq_screen/bin:$PATH

COPY ./fastq_screen.conf /home/

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install wget -y && \
    rm -r /var/lib/apt/lists/*

ENTRYPOINT [ "fastq_screen" ]
