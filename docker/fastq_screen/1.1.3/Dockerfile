FROM stjudecloud/conda-base:1.1.1 AS builder

RUN conda create -n fastq_screen \
    fastq-screen==0.13.0=pl526_1 \
    -y \
    && conda clean --all

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
