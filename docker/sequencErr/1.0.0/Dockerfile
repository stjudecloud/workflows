FROM stjudecloud/conda-base:1.0.0 AS builder

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install build-essential pkg-config cmake zlib1g-dev -y && \
    rm -r /var/lib/apt/lists/*

RUN conda create -n sequencErr \
    meson==0.56.2 \
    ninja==1.10.2 \
    -y \
    && conda clean --all

ARG GITHUB_TOKEN
RUN git clone --depth 1 --branch v2.0.7 https://${GITHUB_TOKEN}@github.com/stjude/sequencErr.git

ENV PATH="/opt/conda/envs/sequencErr/bin:${PATH}"

RUN cd sequencErr && \
    meson --buildtype release build && \
    ninja -C build install

FROM debian:10-slim

COPY --from=builder /usr/local/bin/sequencerr /usr/local/bin/sequencerr

ENTRYPOINT [ "sequencerr" ]
