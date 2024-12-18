FROM quay.io/biocontainers/pandas:2.2.1

COPY --from=scripts --chmod=777 methylation/combine.py /usr/local/bin/combine.py
COPY --from=scripts --chmod=777 methylation/filter.py /usr/local/bin/filter.py