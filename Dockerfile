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
    libcairo2-dev \
    libxt-dev \
    libv8-dev \
    libudunits2-dev \
    libgeos-dev \
    libproj-dev \
    libgdal-dev \
    libmagick++-dev \
    libglpk-dev \
    libmysqlclient-dev \
    libpq-dev \
    libbz2-dev \
    liblzma-dev \
    cmake \
    gfortran \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install(c('rtracklayer'), ask=FALSE, update=FALSE)"
RUN R -e "options(Ncpus = parallel::detectCores()); install.packages('Matrix', repos='https://cloud.r-project.org')"
RUN R -e "options(Ncpus = parallel::detectCores()); pkgs <- c('BPCells','presto','SeuratObject','SeuratDisk','Seurat','optparse','patchwork','scales','stringr','knitr','readr','dplyr','tidyr','ggplot2','purrr','tibble'); install.packages(pkgs, repos='https://cloud.r-project.org', dependencies=TRUE); missing <- setdiff(pkgs, rownames(installed.packages())); if (length(missing) > 0) stop('Failed to install packages: ', paste(missing, collapse=', '))"

WORKDIR /workspace

COPY scripts/qc_summary_report.R /usr/local/bin/qc_summary_report.R
RUN chmod +x /usr/local/bin/qc_summary_report.R

ENTRYPOINT ["Rscript", "/usr/local/bin/qc_summary_report.R"]
CMD ["--help"]
