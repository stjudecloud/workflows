FROM continuumio/miniconda3:4.10.3

RUN echo 'channels:\n  - conda-forge\n  - bioconda\n' > /root/.condarc \
    && conda config --set channel_priority strict \
    && conda create -n bwa \
    bwa==0.7.17 \
    samtools==1.11 \
    -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/bwa/bin/ /opt/conda/envs/bwa/bin/
COPY --from=builder /opt/conda/envs/bwa/lib/ /opt/conda/envs/bwa/lib/
ENV PATH /opt/conda/envs/bwa/bin:$PATH

ENTRYPOINT [ "bwa" ]
