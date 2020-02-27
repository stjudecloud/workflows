FROM conda-base AS builder

RUN conda create -n fastq_screen \
        fastq-screen==0.13.0 \
        -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/fastq_screen/bin /opt/conda/envs/fastq_screen/bin
COPY --from=builder /opt/conda/envs/fastq_screen/lib/ /opt/conda/envs/fastq_screen/lib/
COPY --from=builder /opt/conda/envs/fastq_screen/share/fastq-screen-0.13.0-1/fastq_screen /opt/conda/envs/fastq_screen/bin/fastq_screen
ENV PATH /opt/conda/envs/fastq_screen/bin:$PATH

ENTRYPOINT [ "fastq_screen" ]