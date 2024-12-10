FROM quay.io/biocontainers/bioconductor-minfi:1.48.0--r43hdfd78af_0

RUN R --no-save <<SCRIPT
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos='http://cran.us.r-project.org')
    BiocManager::install("IlluminaHumanMethylationEPICmanifest")
    BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
SCRIPT
