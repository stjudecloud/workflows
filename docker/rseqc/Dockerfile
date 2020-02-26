FROM conda-base

RUN conda create -n rseqc \
        rseqc==3.0.0 \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/rseqc/bin:$PATH