FROM docker.pkg.github.com/stjudecloud/workflows/conda-base:0.1.0

RUN conda create -n picard \
        picard==2.20.2 \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/picard/bin:$PATH

ENTRYPOINT [ "picard" ]