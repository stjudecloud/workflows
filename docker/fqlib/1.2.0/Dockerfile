FROM rust:1.59.0 as fqlib-builder

RUN cargo install \
      --git https://github.com/stjude-rust-labs/fq.git \
      --tag v0.9.1 \
      --root /opt/fqlib/

FROM ubuntu:20.04 as builder

COPY --from=fqlib-builder /opt/fqlib/bin/fq /usr/local/bin/

ENTRYPOINT [ "fq" ]
