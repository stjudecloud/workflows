FROM python:3.13.1-alpine

RUN apk add --no-cache bash

COPY --from=scripts --chmod=777 hic/qc_hic.py /usr/local/bin/qc_hic.py