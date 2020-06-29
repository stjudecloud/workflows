FROM stjudecloud/conda-base:1.0.0

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install gcc zlib1g-dev libbz2-dev libcurl4-openssl-dev build-essential libssl-dev tabix -y && \
    rm -r /var/lib/apt/lists/*

RUN pip install --no-cache-dir ngsderive==1.1.0

ENTRYPOINT [ "ngsderive" ]