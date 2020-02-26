FROM conda-base

RUN conda create -n multiqc \
        multiqc==1.7 \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/multiqc/bin:$PATH

ENTRYPOINT [ "multiqc" ]