FROM stjudecloud/conda-base:1.0.0-alpha AS builder

RUN conda create -n rseqc \
        rseqc==3.0.0 \
        -y \
    && conda clean --all

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/rseqc/bin/ /opt/conda/envs/rseqc/bin/
ENV PATH /opt/conda/envs/rseqc/bin:$PATH