FROM docker.pkg.github.com/stjudecloud/workflows/conda-base:0.1.0

RUN conda create -n fastqc \
        -c bioconda \
        -c conda-forge \
        fastqc==0.11.8 \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/fastqc/bin:$PATH

ENTRYPOINT [ "fastqc" ]