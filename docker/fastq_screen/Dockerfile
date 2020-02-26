FROM sjcloud-conda

RUN conda create -n fastq_screen \
        fastq-screen==0.13.0 \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/fastq_screen/bin:$PATH

ENTRYPOINT [ "fastq_screen" ]