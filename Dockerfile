FROM satijalab/seurat:latest

LABEL org.opencontainers.image.title="scQC-flow QC reporting"
LABEL org.opencontainers.image.description="Reproducible environment for summarising pre/post QC Seurat objects."
LABEL org.opencontainers.image.source="https://github.com/Parisarahbar97/scQC-flow"

ENV RENV_DEFAULT_REPOS=https://cloud.r-project.org
ENV MAKEFLAGS=-j4

RUN R -e "options(Ncpus = parallel::detectCores()); pkgs <- c('Seurat','SeuratObject','SeuratDisk','optparse','patchwork','scales','stringr','knitr','readr','dplyr','tidyr','ggplot2','purrr','tibble'); install.packages(pkgs, repos='https://cloud.r-project.org', dependencies=TRUE); missing <- setdiff(pkgs, rownames(installed.packages())); if (length(missing) > 0) stop('Failed to install packages: ', paste(missing, collapse=', '))"

WORKDIR /workspace

COPY scripts/qc_summary_report.R /usr/local/bin/qc_summary_report.R
RUN chmod +x /usr/local/bin/qc_summary_report.R

ENTRYPOINT ["Rscript", "/usr/local/bin/qc_summary_report.R"]
CMD ["--help"]
