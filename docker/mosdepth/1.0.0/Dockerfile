FROM ghcr.io/stjudecloud/conda-base:1.1.1 AS builder

RUN conda create -n mosdepth \
    mosdepth==0.3.3 \
    -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/mosdepth/bin/ /opt/conda/envs/mosdepth/bin/
ENV PATH /opt/conda/envs/mosdepth/bin:$PATH

ENTRYPOINT [ "mosdepth" ]
