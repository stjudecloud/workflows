FROM ghcr.io/r-hub/r-minimal/r-minimal:4.4.2

RUN apk upgrade --no-cache

RUN R -e 'install.packages("estimate", repos="http://r-forge.r-project.org", dependencies=TRUE)'
