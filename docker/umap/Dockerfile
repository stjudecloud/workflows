FROM python:3.13.1-alpine

# Downstream requires llvm15, but that is out-of-date
# llvm15 installs with the -15 suffix, so we need to symlink
RUN apk add build-base llvm15-dev \
    && ln -s /usr/bin/llvm-config-15 /usr/bin/llvm-config

RUN pip install umap-learn==0.5.7 pandas

COPY --from=scripts --chmod=777 methylation/generate_umap.py /usr/local/bin/generate_umap.py