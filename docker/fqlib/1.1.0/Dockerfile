FROM rust:1.56.0 as fqlib-builder

RUN cargo install \
      --git https://github.com/stjude-rust-labs/fq.git \
      --tag v0.8.0 \
      --root /opt/fqlib/

FROM ubuntu:18.04 as builder

COPY --from=fqlib-builder /opt/fqlib/bin/fq /usr/local/bin/

ENTRYPOINT [ "fq" ]