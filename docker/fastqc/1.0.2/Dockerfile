FROM ghcr.io/stjudecloud/conda-base:branch-gh-packages-1.1.1 AS builder

RUN conda create -n fastqc \
    fastqc==0.11.8 \
    -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/fastqc/bin/ /opt/conda/envs/fastqc/bin/
COPY --from=builder /opt/conda/envs/fastqc/opt/ /opt/conda/envs/fastqc/opt/
COPY --from=builder /opt/conda/envs/fastqc/lib/ /opt/conda/envs/fastqc/lib/
COPY --from=builder /opt/conda/envs/fastqc/etc/fonts/ /opt/conda/envs/fastqc/etc/fonts/
COPY --from=builder /opt/conda/envs/fastqc/fonts/ /opt/conda/envs/fastqc/fonts/
ENV PATH /opt/conda/envs/fastqc/bin:$PATH

ENTRYPOINT [ "fastqc" ]
