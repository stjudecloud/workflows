FROM ubuntu:20.04

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install curl -y && \
    rm -r /var/lib/apt/lists/*

WORKDIR /opt

RUN curl -o cellranger-7.0.1.tar.gz \
    "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.0.1.tar.gz?Expires=1667269124&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjAuMS50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NjcyNjkxMjR9fX1dfQ__&Signature=jTFY9Y0Phw8Sj621d9hgGxVSHt0u0nbiI3D5-~zMkCDqdX4LlW62ObLwTjRFWB5X6E~wXsksJivmS~UphIwyIdTQ62cKqWQ-90779rz1~aVKDuDpe7BRYm9~AqVumyu7d1V7ceOC8OR6zEk-T~7NiN8kFavL3Y9tbTGxV2Sn84aBiiLVGbY295myszqgTRbgfF6AzRYt3c8xybJGk96f02-yM-dp4FwdEvptvBALAVUS~lsdn4Zzx1SY27fsOJvBXM0ZJvSoXF382xN07u6FY4~6euqzq1pJCs~I8t~cpA6jkYxa-l1bkYWplnJ7iVkc-CegMx71M0yCHwiZcPmnsQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" \
    && echo '98ffff631adb7298d38194073626817f  cellranger-7.0.1.tar.gz' > cellranger-7.0.1.tar.gz.md5 \
    && md5sum -c cellranger-7.0.1.tar.gz.md5 \
    && tar -xzvf cellranger-7.0.1.tar.gz

ENV PATH "/opt/cellranger-7.0.1:$PATH"
