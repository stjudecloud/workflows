FROM ghcr.io/stjudecloud/bedtools:2.31.1-0 AS bedtools
FROM python:3.13.1-alpine

RUN apk add --no-cache bash

COPY --from=bedtools /usr/local/bin/ /usr/local/bin/

COPY --from=scripts --chmod=777 hic/qc_hic.py /scripts/hic/qc_hic.py
COPY --from=scripts --chmod=777 hic/filter_hic.py /scripts/hic/filter_hic.py