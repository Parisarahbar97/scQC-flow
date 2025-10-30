# scVI / scANVI helper scripts

These scripts reproduce the full analysis workflow on the latest scQC-flow outputs.
Nothing runs automatically—execute each step manually on the cluster once the
requirements (`scanpy`, `scvi-tools`, `anndata`, `matplotlib`, etc.) are available.

## 1. Build per-sample post-QC AnnData files

```bash
python scripts/scvi_pipeline/01_build_per_sample_h5ad.py \
  --healthy-root /home/pr422/RDS/live/Users/Parisa/EPILEP/healthy/qc/output_latest \
  --diseased-root /home/pr422/RDS/live/Users/Parisa/EPILEP/diseased/qc/output_latest \
  --out-root /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/scvi_integration \
  --clear-inputs
```

Outputs: `inputs/h5ad_files/{healthy,diseased}/*_postqc.h5ad` plus
`inputs/build_summary.tsv`.

## 2. Sanity checks

```bash
python scripts/scvi_pipeline/02_sanity_checks.py \
  --inputs-root /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/scvi_integration/inputs
```

Outputs: `inputs/sanity_report.tsv` and `inputs/sanity_report.txt` with detailed
QC flags (integer-like counts, key metadata columns, etc.).

## 3. Create 2k-HVG dataset (tutorial style)

```bash
python scripts/scvi_pipeline/03_make_hvg2k.py \
  --inputs-root /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/scvi_integration/inputs \
  --out-root /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/scvi_integration
```

Outputs: `outputs_2/adata_hvg2k_TM.h5ad` and `outputs_2/hvg2k_TM_genes.tsv`.

## 4. Train scVI (2k HVGs)

```bash
python scripts/scvi_pipeline/04_run_scvi.py \
  --hvg-path /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/scvi_integration/outputs_2/adata_hvg2k_TM.h5ad \
  --out-root /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/scvi_integration \
  --model-dir /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/scvi_integration/models_2/scvi_TM_hvg2k \
  --max-epochs 150
```

Outputs: `outputs_2/adata_scvi_TM2k_integrated.h5ad`, UMAP plots under
`outputs_2/plots/`, and the saved model in `models_2/scvi_TM_hvg2k/`.

## 5. Run scANVI label transfer

```bash
python scripts/scvi_pipeline/05_run_scanvi.py \
  --reference /mnt/data/pr422/projects/epilep/outputs/scanvi/ref_LIBD_hippocampus_snRNA.h5ad \
  --query /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/scvi_integration/outputs_2/adata_hvg2k_TM.h5ad \
  --out-root /rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/scvi_integration \
  --labels-key cell.type2
```

Outputs: `outputs_2/scanvi/query_scanvi_annotated.h5ad`, joint object,
confidence summary, and UMAP/diagnostic plots in
`outputs_2/scanvi/plots/`.

Feel free to tweak epochs, latent dimensions, or UMAP/Leiden settings via the
optional CLI flags provided in each script.

