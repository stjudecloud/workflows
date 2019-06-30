FROM ubuntu:18.04

ENV MINICONDA_URL "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install wget -y

RUN wget ${MINICONDA_URL} -O miniconda.sh
RUN bash miniconda.sh -b -p /opt/conda/
RUN rm -r miniconda.sh

RUN conda update -n base -c defaults conda -y
RUN conda install \
    -c conda-forge \
    -c bioconda \
    coreutils==8.31 \
    picard==2.20.2 \
    samtools==1.9 \
    bwa==0.7.17 \
    star==2.7.1a \
    fastqc==0.11.8 \
    qualimap==2.2.2c \
    multiqc==1.7 \
    rseqc==3.0.0 \
    -y