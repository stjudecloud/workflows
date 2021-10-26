FROM ghcr.io/stjudecloud/conda-base:branch-gh-packages-1.1.1 AS builder

RUN conda create -n htseq \
    htseq==0.11.2 \
    -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/htseq/bin/ /opt/conda/envs/htseq/bin/
COPY --from=builder /opt/conda/envs/htseq/lib/ /opt/conda/envs/htseq/lib/
ENV PATH /opt/conda/envs/htseq/bin:$PATH

ENTRYPOINT [ "htseq-count" ]
