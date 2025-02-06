FROM quay.io/biocontainers/samtools:1.17--h00cdaf9_0 AS samtools
FROM quay.io/biocontainers/bwa:0.7.17--h7132678_9

COPY --from=samtools /usr/local/bin/ /usr/local/bin/
COPY --from=samtools /usr/local/lib/ /usr/local/lib/
COPY --from=samtools /usr/local/libexec/ /usr/local/libexec/

ENTRYPOINT [ "bwa" ]
