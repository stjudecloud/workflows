FROM stjudecloud/conda-base:1.1.0 AS builder

RUN conda create -n qualimap \
    qualimap==2.2.2c=1 \
    -y \
    && conda clean --all

FROM debian:10-slim
RUN apt-get update && \
    apt-get install -y libfontconfig && \
    rm -r /var/lib/apt/lists/*
COPY --from=builder /opt/conda/envs/qualimap/bin/ /opt/conda/envs/qualimap/bin/
COPY --from=builder /opt/conda/envs/qualimap/lib/ /opt/conda/envs/qualimap/lib/
COPY --from=builder /opt/conda/envs/qualimap/share/qualimap-2.2.2c-1/ /opt/conda/envs/qualimap/share/qualimap-2.2.2c-1/
COPY --from=builder /opt/conda/envs/qualimap/etc/fonts/ /opt/conda/envs/qualimap/etc/fonts/

ENV PATH /opt/conda/envs/qualimap/bin:$PATH

ENTRYPOINT [ "qualimap" ]