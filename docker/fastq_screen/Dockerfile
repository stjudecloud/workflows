FROM continuumio/miniconda3:4.7.12-alpine

RUN /opt/conda/bin/conda create -n fastq_screen \
                 -c bioconda \
                 -c conda-forge \
                 fastq-screen==0.13.0 \
                 -y \
    && /opt/conda/bin/conda clean --all
ENV PATH /opt/conda/envs/fastq_screen/bin:$PATH

ENTRYPOINT [ "/usr/bin/env", "bash", "-c", "fastq_screen"]