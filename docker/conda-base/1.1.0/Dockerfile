FROM continuumio/miniconda3:4.7.12

RUN apt-get --allow-releaseinfo-change update && \
    apt-get upgrade -y && \
    rm -r /var/lib/apt/lists/*

RUN echo 'channels:\n  - conda-forge\n  - bioconda\n' > /root/.condarc \
    && conda config --set channel_priority strict