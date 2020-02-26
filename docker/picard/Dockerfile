FROM conda-base

RUN conda create -n picard \
        picard==2.20.2 \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/picard/bin:$PATH

ENTRYPOINT [ "picard" ]