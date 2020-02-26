FROM docker.pkg.github.com/stjudecloud/workflows/conda-base:0.1.0

RUN conda create -n fastq_screen \
        fastq-screen==0.13.0 \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/fastq_screen/bin:$PATH

ENTRYPOINT [ "fastq_screen" ]