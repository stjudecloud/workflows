FROM quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2 AS bedtools
FROM python:3.13.0

COPY --from=bedtools /usr/local/bin/ /usr/local/bin/
COPY --from=scripts --chmod=777 hic/filter_hic.py /usr/local/bin/filter_hic.py

ENTRYPOINT [ "bash" ]