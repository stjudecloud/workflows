FROM quay.io/biocontainers/bioconductor-minfi:1.48.0--r43hdfd78af_0

RUN R --no-save <<SCRIPT
    install.packages("argparse", repos='http://cran.us.r-project.org')
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos='http://cran.us.r-project.org')
    BiocManager::install("IlluminaHumanMethylationEPICmanifest")
    BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
SCRIPT

COPY --from=scripts --chmod=777 methylation/methylation-preprocess.R /scripts/methylation/methylation-preprocess.R