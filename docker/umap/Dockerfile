FROM python:3.14.3-alpine

# Downstream requires llvm20, but that is out-of-date
# llvm20 installs with the -20 suffix, so we need to symlink
RUN apk add --no-cache build-base llvm20-dev bash cmake llvm20 llvm20-libs llvm20-static llvm20-gtest \
    && ln -s /usr/bin/llvm20-config /usr/bin/llvm-config

RUN pip install umap-learn==0.5.11 pandas

COPY --from=scripts --chmod=777 methylation/generate_umap.py /scripts/methylation/generate_umap.py