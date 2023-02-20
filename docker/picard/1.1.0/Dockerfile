FROM ghcr.io/stjudecloud/conda-base:1.1.1 AS builder

RUN conda create -n picard \
    picard==2.27.5=hdfd78af_0 \
    -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/picard/bin/ /opt/conda/envs/picard/bin/
COPY --from=builder /opt/conda/envs/picard/lib/ /opt/conda/envs/picard/lib/
COPY --from=builder /opt/conda/envs/picard/share/picard-2.27.5-0/ /opt/conda/envs/picard/share/picard-2.27.5-0/
ENV PATH /opt/conda/envs/picard/bin:$PATH

ENTRYPOINT [ "picard" ]
