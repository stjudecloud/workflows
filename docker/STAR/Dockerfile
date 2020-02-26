FROM docker.pkg.github.com/stjudecloud/workflows/conda-base:0.1.0

RUN conda create -n star \
        star==2.7.1a \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/star/bin:$PATH

ENTRYPOINT [ "STAR" ]