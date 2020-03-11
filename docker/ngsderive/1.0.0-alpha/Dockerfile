FROM stjudecloud/conda-base:1.0.0-alpha

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install gcc zlib1g-dev libbz2-dev libcurl4-openssl-dev build-essential libssl-dev -y && \
    rm -r /var/lib/apt/lists/*

RUN pip install --no-cache-dir ngsderive==1.0.2

ENTRYPOINT [ "ngsderive" ]