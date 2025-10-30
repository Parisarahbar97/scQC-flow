#!/usr/bin/env python3
import argparse
from pathlib import Path
import scanpy as sc
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description='Make 2k HVG dataset for scVI (tutorial-style)')
    ap.add_argument('--inputs-root', type=Path, required=True, help='scvi_integration/inputs root path')
    ap.add_argument('--out-root', type=Path, required=True, help='Base scvi_integration path (will write outputs_2/)')
    ap.add_argument('--n-top', type=int, default=2000)
    args = ap.parse_args()

    files = sorted((args.inputs_root / 'h5ad_files').glob('**/*_postqc.h5ad'))
    if not files:
        raise SystemExit('No *_postqc.h5ad files found under inputs/h5ad_files')

    # Load and concat
    adatas = [sc.read_h5ad(f) for f in files]
    # Ensure 'sample' exists
    for a in adatas:
        if 'sample' not in a.obs:
            a.obs['sample'] = Path(a.uns.get('sample_name', 'unknown')).name if 'sample_name' in a.uns else 'unknown'
    adata = adatas[0].concatenate(*adatas[1:], batch_key='concat_batch', batch_categories=[p.stem for p in files], index_unique=None)

    # Use provided sample column as batch for HVG selection
    if 'sample' not in adata.obs:
        adata.obs['sample'] = adata.obs['concat_batch']

    # Keep counts in layers['counts'] as the raw counts source
    if 'counts' not in adata.layers:
        adata.layers['counts'] = adata.X.copy()

    # HVGs tutorial-style (Seurat v3 flavor, per-sample)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=args.n_top,
        flavor='seurat_v3',
        batch_key='sample',
        inplace=True
    )

    hvg_mask = adata.var['highly_variable'].values
    adata_hvg = adata[:, hvg_mask].copy()

    out2 = args.out_root / 'outputs_2'
    out2.mkdir(parents=True, exist_ok=True)

    adata_hvg.write_h5ad(out2 / 'adata_hvg2k_TM.h5ad')
    pd.Series(adata_hvg.var_names, name='gene').to_csv(out2 / 'hvg2k_TM_genes.tsv', index=False)
    print(f"Saved {adata_hvg.n_obs} cells x {adata_hvg.n_vars} genes to {out2 / 'adata_hvg2k_TM.h5ad'}")


if __name__ == '__main__':
    raise SystemExit(main())

