FROM ghcr.io/stjudecloud/conda-base:branch-gh-packages-1.1.1 AS builder

RUN conda create -n samtools \
    samtools==1.9 \
    -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/samtools/bin/ /opt/conda/envs/samtools/bin/
COPY --from=builder /opt/conda/envs/samtools/lib/ /opt/conda/envs/samtools/lib/
ENV PATH /opt/conda/envs/samtools/bin:$PATH

ENTRYPOINT [ "samtools" ]
