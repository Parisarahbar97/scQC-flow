FROM satijalab/seurat:4.3.0

LABEL org.opencontainers.image.title="scQC-flow QC reporting"
LABEL org.opencontainers.image.description="Reproducible environment for summarising pre/post QC Seurat objects."
LABEL org.opencontainers.image.source="https://github.com/Parisarahbar97/scQC-flow"

ENV RENV_DEFAULT_REPOS=https://cloud.r-project.org

RUN R -e "pkgs <- c('optparse','patchwork','scales','stringr','knitr'); to_install <- setdiff(pkgs, rownames(installed.packages())); if (length(to_install) > 0) install.packages(to_install, repos='https://cloud.r-project.org')"

WORKDIR /workspace

COPY scripts/qc_summary_report.R /usr/local/bin/qc_summary_report.R
RUN chmod +x /usr/local/bin/qc_summary_report.R

ENTRYPOINT [\"Rscript\", \"/usr/local/bin/qc_summary_report.R\"]
CMD [\"--help\"]
