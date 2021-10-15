FROM continuumio/miniconda3:4.10.3

RUN apt-get update && \
    apt-get upgrade -y && \
    rm -r /var/lib/apt/lists/*

RUN echo 'channels:\n  - conda-forge\n  - bioconda\n' > /root/.condarc \
    && conda config --set channel_priority strict