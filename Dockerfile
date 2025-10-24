FROM rocker/tidyverse:4.3.1

LABEL org.opencontainers.image.title="scQC-flow QC reporting"
LABEL org.opencontainers.image.description="Reproducible environment for summarising pre/post QC Seurat objects."
LABEL org.opencontainers.image.source="https://github.com/Parisarahbar97/scQC-flow"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

ENV RENV_DEFAULT_REPOS=https://cloud.r-project.org
ENV CRAN_MIRROR=https://packagemanager.posit.co/cran/__linux__/jammy/latest

RUN R -e "options(repos='${CRAN_MIRROR}'); install.packages(c('optparse','Seurat','SeuratDisk','patchwork','scales','stringr','knitr'), dependencies=TRUE)"

WORKDIR /workspace

COPY scripts/qc_summary_report.R /usr/local/bin/qc_summary_report.R
RUN chmod +x /usr/local/bin/qc_summary_report.R

ENTRYPOINT ["Rscript", "/usr/local/bin/qc_summary_report.R"]
CMD ["--help"]
