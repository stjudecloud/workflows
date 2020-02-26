FROM sjcloud-conda

RUN conda create -n htseq \
        htseq==0.11.2 \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/htseq/bin:$PATH

ENTRYPOINT [ "htseq-count" ]