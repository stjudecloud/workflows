FROM continuumio/miniconda3:latest

RUN /opt/conda/bin/conda create -n fastq_screen \
                 -c bioconda \
                 -c conda-forge \
                 fastq-screen==0.13.0 \
                 -y \
    && /opt/conda/bin/conda clean --all
ENV PATH /opt/conda/envs/fastq_screen/bin:$PATH

CMD [ "-h" ]
ENTRYPOINT [ "fastq_screen" ]