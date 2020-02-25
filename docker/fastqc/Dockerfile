FROM continuumio/miniconda3:4.7.12

RUN conda create -n fastqc \
                -c bioconda \
                -c conda-forge \
                fastqc==0.11.8 \
                -y \
    && conda clean --all
ENV PATH /opt/conda/envs/fastqc/bin:$PATH

ENTRYPOINT [ "fastqc" ]