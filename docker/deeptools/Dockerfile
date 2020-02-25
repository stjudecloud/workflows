FROM sjcloud-conda

RUN conda create -n deeptools \
                 deeptools==3.3.1 \
                 -y \
    && conda clean --all
ENV PATH /opt/conda/envs/deeptools/bin:$PATH