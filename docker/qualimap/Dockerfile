FROM docker.pkg.github.com/stjudecloud/workflows/conda-base:0.1.0

RUN conda create -n qualimap \
        qualimap==2.2.2c \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/qualimap/bin:$PATH

ENTRYPOINT [ "qualimap" ]