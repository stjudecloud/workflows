FROM stjudecloud/conda-base:1.1.1 AS builder

RUN conda create -n bwa \
    bwa==0.7.17 \
    samtools==1.11 \
    -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/bwa/bin/ /opt/conda/envs/bwa/bin/
COPY --from=builder /opt/conda/envs/bwa/lib/ /opt/conda/envs/bwa/lib/
ENV PATH /opt/conda/envs/bwa/bin:$PATH

ENTRYPOINT [ "bwa" ]
