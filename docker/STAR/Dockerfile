FROM continuumio/miniconda3:4.7.12

RUN /opt/conda/bin/conda create -n star \
                 -c bioconda \
                 star==2.7.1a -y \
    && /opt/conda/bin/conda clean --all
ENV PATH /opt/conda/envs/star/bin:$PATH

ENTRYPOINT [ "STAR" ]