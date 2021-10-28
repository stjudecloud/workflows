FROM ghcr.io/stjudecloud/conda-base:branch-gh-packages-1.1.1 AS builder

RUN conda create -n multiqc \
    multiqc==1.10.1 \
    -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/multiqc/bin/ /opt/conda/envs/multiqc/bin/
COPY --from=builder /opt/conda/envs/multiqc/lib/ /opt/conda/envs/multiqc/lib/
COPY .multiqc_config.yaml /home/
ENV PATH /opt/conda/envs/multiqc/bin:$PATH

ENTRYPOINT [ "multiqc" ]
