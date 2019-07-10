FROM stjudecloud/bio-base:bleeding-edge

ENV GRCH38_FA_URL "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" 
ENV GENCODE_v30_GTF_URL "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz"

# Create genomics reference directory
RUN mkdir -p /opt/genomics/reference
WORKDIR /opt/genomics/reference/

# Download GRCh38_no_alt FASTA and verify contents
RUN wget ${GRCH38_FA_URL} -O GRCh38_no_alt.fa.gz && \
    gunzip GRCh38_no_alt.fa.gz && \
    echo "a6da8681616c05eb542f1d91606a7b2f  GRCh38_no_alt.fa" > GRCh38_no_alt.fa.md5 && \
    md5sum -c GRCh38_no_alt.fa.md5

# Download GENCODE v30 GTF and verify contents
RUN wget ${GENCODE_v30_GTF_URL} -O gencode.v30.annotation.gtf.gz && \
    gunzip gencode.v30.annotation.gtf.gz && \
    echo "63770a3d2c6adb4d9d1bbc9ba3bd4adf  gencode.v30.annotation.gtf" > gencode.v30.annotation.gtf.md5 && \
    md5sum -c gencode.v30.annotation.gtf.md5

# Generate STAR genome dir "STARDB"
RUN mkdir STARDB/ && \
    STAR --runMode genomeGenerate \
         --genomeDir STARDB/ \
         --runThreadN `nproc` \
         --genomeFastaFiles GRCh38_no_alt.fa \
         --sjdbGTFfile gencode.v30.annotation.gtf \
         --limitGenomeGenerateRAM 9000000000 \
         --sjdbOverhang 125 && \
    rm -rf GRCh38_no_alt* gencode.v30*