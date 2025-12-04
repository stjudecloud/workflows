FROM quay.io/biocontainers/samtools:1.17--h00cdaf9_0 AS samtools
FROM quay.io/biocontainers/hisat2:2.2.1--h503566f_8

COPY --from=samtools /usr/local/bin/ /usr/local/bin/
COPY --from=samtools /usr/local/lib/ /usr/local/lib/
COPY --from=samtools /usr/local/libexec/ /usr/local/libexec/

ENTRYPOINT [ "hisat2" ]