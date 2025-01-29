FROM python:3.13.0

RUN pip install pandas numpy matplotlib seaborn plotly

COPY --from=scripts --chmod=777 methylation/plot_umap.py /usr/local/bin/plot_umap.py