FROM docker.pkg.github.com/stjudecloud/workflows/conda-base:0.1.0

RUN conda create -n multiqc \
        multiqc==1.7 \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/multiqc/bin:$PATH

ENTRYPOINT [ "multiqc" ]