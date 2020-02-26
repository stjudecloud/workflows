FROM conda-base

RUN conda create -n star \
        star==2.7.1a -y \
    && conda clean --all
ENV PATH /opt/conda/envs/star/bin:$PATH

ENTRYPOINT [ "STAR" ]