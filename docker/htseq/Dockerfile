FROM docker.pkg.github.com/stjudecloud/workflows/conda-base:0.1.0

RUN conda create -n htseq \
        htseq==0.11.2 \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/htseq/bin:$PATH

ENTRYPOINT [ "htseq-count" ]