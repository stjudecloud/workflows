FROM quay.io/biocontainers/samtools:1.17--h00cdaf9_0 AS samtools
FROM quay.io/biocontainers/minimap2:2.30--h577a1d6_0

COPY --from=samtools /usr/local/bin/ /usr/local/bin/
COPY --from=samtools /usr/local/lib/ /usr/local/lib/
COPY --from=samtools /usr/local/libexec/ /usr/local/libexec/

ENTRYPOINT [ "minimap2" ]