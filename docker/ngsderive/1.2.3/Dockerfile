FROM ghcr.io/stjudecloud/conda-base:1.1.1 AS builder

RUN conda create -n ngsderive \
    ngsderive==2.2.0 \
    tabix==1.11 \
    -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/ngsderive/bin/ /opt/conda/envs/ngsderive/bin/
COPY --from=builder /opt/conda/envs/ngsderive/lib/ /opt/conda/envs/ngsderive/lib/
ENV PATH /opt/conda/envs/ngsderive/bin:$PATH

ENTRYPOINT [ "ngsderive" ]
