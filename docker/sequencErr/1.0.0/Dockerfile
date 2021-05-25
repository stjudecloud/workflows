# syntax=docker/dockerfile:1.0.0-experimental

FROM centos:7 AS builder

RUN yum -y install gcc-c++ git openssh-clients python3 zlib-devel \
    && yum clean all

RUN pip3 install meson==0.57.2 ninja==1.10.0.post2

RUN --mount=type=ssh \
    cd /opt \
    && mkdir -p -m 0600 ~/.ssh \
    && ssh-keyscan github.com >> ~/.ssh/known_hosts \
    && git clone --depth 1 --branch v2.0.7 git@github.com:stjude/sequencErr.git \
    && cd sequencErr \
    && meson --buildtype release build \
    && ninja -C build install

FROM centos:7

COPY --from=builder /opt/sequencErr/build/sequencerr /opt/sequencerr/bin/

ENV PATH=/opt/sequencerr/bin:$PATH

ENTRYPOINT ["/opt/sequencerr/bin/sequencerr"]
