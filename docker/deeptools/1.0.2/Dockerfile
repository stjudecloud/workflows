FROM stjudecloud/conda-base:1.1.1 AS builder

RUN conda create -n deeptools \
    deeptools==3.3.1 \
    -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/deeptools/bin /opt/conda/envs/deeptools/bin
COPY --from=builder /opt/conda/envs/deeptools/lib/ /opt/conda/envs/deeptools/lib/
ENV PATH /opt/conda/envs/deeptools/bin:$PATH
