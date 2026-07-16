# Reproducible environment for the PAH lung transcriptomic meta-analysis pipeline.
# Pins R 4.6.1 (Bioconductor 3.22/3.23) and restores the exact package library from renv.lock.
FROM rocker/r-ver:4.6.1

# system libraries required by the Bioconductor / single-cell / plotting stack
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev libcurl4-openssl-dev libssl-dev libhdf5-dev \
    libpng-dev libjpeg-dev libtiff5-dev libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev zlib1g-dev \
    libglpk-dev libgsl-dev cmake && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /pipeline
COPY renv.lock renv.lock

# restore the exact pinned package set (CRAN + Bioconductor) recorded at analysis time
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org'); \
          renv::restore(lockfile='renv.lock', prompt=FALSE)"

# the scripts are copied into the image; the GEO input data is mounted at run time
COPY . /pipeline
# example:
#   docker run -v /path/to/project:/data -w /data <image> Rscript submission/scripts/03a_meta_analysis.R
CMD ["R"]
