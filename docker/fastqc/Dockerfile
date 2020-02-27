FROM conda-base AS builder

RUN conda create -n fastqc \
        fastqc==0.11.8 \
        -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/fastqc/bin/ /opt/conda/envs/fastqc/bin/
COPY --from=builder /opt/conda/envs/fastqc/opt/ /opt/conda/envs/fastqc/opt/
COPY --from=builder /opt/conda/envs/fastqc/lib/ /opt/conda/envs/fastqc/lib/
ENV PATH /opt/conda/envs/fastqc/bin:$PATH

ENTRYPOINT [ "fastqc" ]