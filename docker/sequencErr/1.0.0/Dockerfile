FROM registry.access.redhat.com/ubi7/python-36 AS builder

ARG GITHUB_TOKEN
RUN git clone --depth 1 --branch v2.0.7 https://${GITHUB_TOKEN}@github.com/stjude/sequencErr.git

RUN pip3 install meson==0.56.2 ninja==1.10

RUN cd sequencErr && \
    meson -Dprefix=$HOME --buildtype release build && \
    ninja -C build install && \
    ./build/sequencerr_tests

FROM registry.access.redhat.com/ubi7/ubi-minimal

COPY --from=builder /opt/app-root/src/bin/sequencerr /usr/local/bin/sequencerr

ENTRYPOINT [ "sequencerr" ]
