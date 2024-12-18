FROM python:3.12.7

RUN pip install umap-learn==0.5.7 pandas

COPY --from=scripts --chmod=777 methylation/generate_umap.py /usr/local/bin/generate_umap.py