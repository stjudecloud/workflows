FROM quay.io/biocontainers/pandas:2.2.1

COPY --from=scripts --chmod=777 methylation/combine.py /scripts/methylation/combine.py
COPY --from=scripts --chmod=777 methylation/filter.py /scripts/methylation/filter.py