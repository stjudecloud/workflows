FROM ubuntu:20.04

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install curl -y && \
    rm -r /var/lib/apt/lists/*

WORKDIR /opt

RUN curl -o cellranger-7.0.0.tar.gz \
    "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.0.0.tar.gz?Expires=1659509820&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NTk1MDk4MjB9fX1dfQ__&Signature=oFGrWvbHjvKV3WJX2KCypqmAPVmYtD6204Onee1Qlxwyqvh-2h9Mjny-M4sdNsZMSnvBZV-PlINrXi00FAg~padTWrZe7f-gLcxMphyGDLjmUBm~ae2qkfFTnYlSwU76MD3zkR3hy-B6QawuGYQXHYVkcJe1632GG~Nf-u1QgIBvYhQuE6jEE8KcHCu0t80BT18yfu9Yo1XV6VAij13umUjD12vwJes5giJ0XeuBNRTIMCvAzucsQXnwCXSD7wnjyEMVlX701JeT87tTZGY4a8mFoCyfM9AqArBhipvz4lOe~LpVASBgteaiPNvaIRpxpPtAlwrvg5B7mxSEoIBMpQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" \
    && echo '30855cb96a097c9cab6b02bdb520423f  cellranger-7.0.0.tar.gz' > cellranger-7.0.0.tar.gz.md5 \
    && md5sum -c cellranger-7.0.0.tar.gz.md5 \
    && tar -xzvf cellranger-7.0.0.tar.gz 

ENV PATH "/opt/cellranger-7.0.0:$PATH"