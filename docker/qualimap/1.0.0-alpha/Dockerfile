FROM stjudecloud/conda-base:1.0.0-alpha AS builder

RUN conda create -n qualimap \
        qualimap==2.2.2c \
        -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/qualimap/bin/ /opt/conda/envs/qualimap/bin/
COPY --from=builder /opt/conda/envs/qualimap/lib/ /opt/conda/envs/qualimap/lib/
COPY --from=builder /opt/conda/envs/qualimap/share/qualimap-2.2.2c-1/ /opt/conda/envs/qualimap/share/qualimap-2.2.2c-1/
COPY --from=builder /opt/conda/envs/qualimap/etc/fonts/ /opt/conda/envs/qualimap/etc/fonts/
COPY --from=builder /opt/conda/envs/qualimap/fonts/ /opt/conda/envs/qualimap/fonts/
ENV PATH /opt/conda/envs/qualimap/bin:$PATH

ENTRYPOINT [ "qualimap" ]