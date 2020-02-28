FROM stjudecloud/conda-base:1.0.0-alpha AS builder

RUN conda create -n star \
        star==2.7.1a \
        -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/star/bin /opt/conda/envs/star/bin
ENV PATH /opt/conda/envs/star/bin:$PATH
ENTRYPOINT [ "STAR" ]