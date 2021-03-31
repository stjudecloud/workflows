FROM rust:1.50.0 as fqlib-builder

RUN cargo install \
      --git https://github.com/stjude/fqlib.git \
      --tag v0.6.0 \
      --root /opt/fqlib/

FROM ubuntu:18.04 as builder

COPY --from=fqlib-builder /opt/fqlib/bin/fq /usr/local/bin/

ENTRYPOINT [ "fq" ]