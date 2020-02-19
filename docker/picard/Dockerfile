FROM continuumio/miniconda3:4.7.12

RUN /opt/conda/bin/conda create -n picard \
                 -c bioconda \
                 -c conda-forge \
                 picard==2.20.2 \
                 bash \
                 -y \
    && /opt/conda/bin/conda clean --all
ENV PATH /opt/conda/envs/picard/bin:$PATH

ENTRYPOINT [ "/usr/bin/env", "bash", "-c", "bash picard" ]