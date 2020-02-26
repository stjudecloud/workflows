FROM conda-base

RUN conda create -n qualimap \
        qualimap==2.2.2c \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/qualimap/bin:$PATH

ENTRYPOINT [ "qualimap" ]