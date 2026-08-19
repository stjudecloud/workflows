# Salmon test fixtures

`transcripts.fasta` — synthetic transcripts built by concatenating real R1 + reverse-complemented R2 sequences from the shared `fastqs/test_R1.fq.gz`/`fastqs/test_R2.fq.gz` fixtures, ensuring genuine alignment for testing purposes.

`salmon_index.tar.gz` — built by running the `build_salmon_index` task against `transcripts.fasta` above.
