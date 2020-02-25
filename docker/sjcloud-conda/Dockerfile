FROM continuumio/miniconda3:4.7.12

RUN echo 'channels:\n  - conda-forge\n  - bioconda\n' > /root/.condarc \
    && conda config --set channel_priority strict