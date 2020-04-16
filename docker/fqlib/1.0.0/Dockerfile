FROM rust:1.40.0 as fqlib-builder

RUN cargo install \
      --git https://github.com/stjude/fqlib.git \
      --root /opt/fqlib/

FROM ubuntu:18.04 as builder

COPY --from=fqlib-builder /opt/fqlib/bin/fq /usr/local/bin/

ENTRYPOINT [ "fq" ]