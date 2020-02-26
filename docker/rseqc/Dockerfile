FROM docker.pkg.github.com/stjudecloud/workflows/conda-base:0.1.0

RUN conda create -n rseqc \
        rseqc==3.0.0 \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/rseqc/bin:$PATH