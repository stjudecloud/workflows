FROM docker.pkg.github.com/stjudecloud/workflows/conda-base:0.1.0

RUN conda create -n samtools \
        samtools==1.9 \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/samtools/bin:$PATH

ENTRYPOINT [ "samtools" ]