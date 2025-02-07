FROM ghcr.io/stjudecloud/bedtools:branch-hic_workflow-2.31.1-0 AS bedtools
FROM python:3.13.1-alpine

RUN apk add --no-cache bash

COPY --from=bedtools /usr/local/bin/ /usr/local/bin/

COPY --from=scripts --chmod=777 hic/qc_hic.py /usr/local/bin/qc_hic.py
COPY --from=scripts --chmod=777 hic/filter_hic.py /usr/local/bin/filter_hic.py