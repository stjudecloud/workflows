FROM ghcr.io/stjudecloud/conda-base:1.1.1 AS builder

RUN conda create -n gtfparse \
    gtfparse==1.2.1 \
    -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/gtfparse/bin/ /opt/conda/envs/gtfparse/bin/
COPY --from=builder /opt/conda/envs/gtfparse/lib/ /opt/conda/envs/gtfparse/lib/
ENV PATH /opt/conda/envs/gtfparse/bin:$PATH

ENTRYPOINT [ "python" ]
