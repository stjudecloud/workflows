FROM docker.pkg.github.com/stjudecloud/workflows/conda-base:0.1.0

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install gcc zlib1g-dev libbz2-dev libcurl4-openssl-dev build-essential libssl-dev -y && \
    rm -r /var/lib/apt/lists/*

RUN conda create -n ngsderive \
        -c anaconda \
        python=3.7 \ 
        cython \
        -y \
    && conda clean --all
ENV PATH /opt/conda/envs/ngsderive/bin:$PATH

RUN pip install --no-cache-dir ngsderive==1.0.1

ENTRYPOINT [ "ngsderive" ]