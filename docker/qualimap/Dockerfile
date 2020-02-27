FROM conda-base AS builder

RUN conda create -n qualimap \
        qualimap==2.2.2c \
        -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/qualimap/bin/ /opt/conda/envs/qualimap/bin/
COPY --from=builder /opt/conda/envs/qualimap/lib/ /opt/conda/envs/qualimap/lib/
COPY --from=builder /opt/conda/envs/qualimap/share/qualimap-2.2.2c-1/ /opt/conda/envs/qualimap/share/qualimap-2.2.2c-1/
ENV PATH /opt/conda/envs/qualimap/bin:$PATH

ENTRYPOINT [ "qualimap" ]