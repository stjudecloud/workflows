FROM ghcr.io/stjudecloud/conda-base:branch-gh-packages-1.1.1 AS builder

RUN conda create -n picard \
    picard==2.20.2=0 \
    -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/picard/bin/ /opt/conda/envs/picard/bin/
COPY --from=builder /opt/conda/envs/picard/lib/ /opt/conda/envs/picard/lib/
COPY --from=builder /opt/conda/envs/picard/share/picard-2.20.2-0/ /opt/conda/envs/picard/share/picard-2.20.2-0/
ENV PATH /opt/conda/envs/picard/bin:$PATH

ENTRYPOINT [ "picard" ]
