FROM quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2 AS bedtools
FROM python:3.9.19

COPY --from=bedtools /usr/local/bin/ /usr/local/bin/

ENTRYPOINT [ "bash" ]