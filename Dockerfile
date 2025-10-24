FROM rocker/tidyverse:4.3.1

LABEL org.opencontainers.image.title="scQC-flow QC reporting"
LABEL org.opencontainers.image.description="Reproducible environment for summarising pre/post QC Seurat objects."
LABEL org.opencontainers.image.source="https://github.com/Parisarahbar97/scQC-flow"

ENV DEBIAN_FRONTEND=noninteractive
ENV RENV_DEFAULT_REPOS=https://cloud.r-project.org
ENV MAKEFLAGS=-j4

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libfftw3-dev \
    libgsl-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libfreetype6-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "options(Ncpus = parallel::detectCores()); pkgs <- c('Seurat','SeuratDisk','optparse','patchwork','scales','stringr','knitr','readr','dplyr','tidyr','ggplot2','purrr','tibble'); install.packages(pkgs, repos='https://cloud.r-project.org', dependencies=TRUE)"

WORKDIR /workspace

COPY scripts/qc_summary_report.R /usr/local/bin/qc_summary_report.R
RUN chmod +x /usr/local/bin/qc_summary_report.R

ENTRYPOINT ["Rscript", "/usr/local/bin/qc_summary_report.R"]
CMD ["--help"]
