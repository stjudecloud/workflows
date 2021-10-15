FROM stjudecloud/conda-base:1.1.1 AS builder

RUN conda create -n star \
    star==2.7.1a \
    -y \
    && conda clean --all

ENV PATH /opt/conda/envs/star/bin:$PATH

COPY ./sort_star_input.py /home

ENTRYPOINT [ "STAR" ]
