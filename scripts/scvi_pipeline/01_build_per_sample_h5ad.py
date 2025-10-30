#!/usr/bin/env python3
import argparse
from pathlib import Path
import sys
import os
import json
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def is_integer_like_matrix(mtx, check_n=100000):
    """Heuristic: sample up to check_n entries and verify they are close to integers."""
    if sparse.issparse(mtx):
        data = mtx.data
        if data.size == 0:
            return True
        if data.size > check_n:
            idx = np.random.choice(data.size, size=check_n, replace=False)
            data = data[idx]
    else:
        flat = mtx.ravel()
        if flat.size == 0:
            return True
        if flat.size > check_n:
            idx = np.random.choice(flat.size, size=check_n, replace=False)
            data = flat[idx]
        else:
            data = flat
    return np.allclose(data, np.round(data))


def load_counts_from_sample_dir(sample_dir: Path, sample: str):
    """Load counts AnnData for a sample from preferred sources.

    Preference order:
    1) <sample>_cellbender_output_seurat.h5 (10x-like H5)
    2) <sample>_cellbender_output/cellbender_out_filtered.h5
    3) <sample>_cellbender_output/cellbender_out.h5
    """
    # 1) seurat-compatible H5
    h5 = sample_dir / f"{sample}_cellbender_output_seurat.h5"
    if h5.exists():
        return sc.read_10x_h5(h5)

    # 2) cellbender filtered raw
    cb_dir = sample_dir / f"{sample}_cellbender_output"
    h5b = cb_dir / "cellbender_out_filtered.h5"
    if h5b.exists():
        return sc.read_10x_h5(h5b)

    # 3) cellbender unfiltered output
    h5c = cb_dir / "cellbender_out.h5"
    if h5c.exists():
        return sc.read_10x_h5(h5c)

    raise FileNotFoundError(
        f"No counts H5 found for {sample} in {sample_dir}. Expected one of: "
        f"{h5.name}, {h5b}, {h5c}"
    )


def attach_qc_and_filter(adata, sample_dir: Path, sample: str,
                         min_features: int, min_counts: int,
                         max_mito: float, min_nuclear: float):
    # Ensure var and obs names unique
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # Compute counts-based QC
    counts = adata.X
    if sparse.issparse(counts):
        total_counts = np.asarray(counts.sum(axis=1)).ravel()
        n_genes = np.asarray((counts > 0).sum(axis=1)).ravel()
    else:
        total_counts = counts.sum(axis=1)
        n_genes = (counts > 0).sum(axis=1)
    adata.obs['total_counts'] = total_counts
    adata.obs['n_genes_by_counts'] = n_genes

    # Mito percent
    mito_mask = adata.var_names.str.upper().str.startswith('MT-')
    if mito_mask.any():
        if sparse.issparse(counts):
            mito_counts = np.asarray(counts[:, mito_mask].sum(axis=1)).ravel()
        else:
            mito_counts = counts[:, mito_mask].sum(axis=1)
        adata.obs['percent_mt'] = 100.0 * mito_counts / np.maximum(total_counts, 1)
    else:
        adata.obs['percent_mt'] = np.nan

    # Attach DropletQC nuclear_fraction if present
    dq_csv = sample_dir / f"{sample}_dropletqc_metrics.csv"
    if dq_csv.exists():
        dq = pd.read_csv(dq_csv)
        # Harmonize barcode column
        bc_col = None
        for c in ['cell_barcode', 'barcode', 'Barcode', 'cell', 'cell_id']:
            if c in dq.columns:
                bc_col = c
                break
        if bc_col is None:
            print(f"Warning: No barcode column found in {dq_csv}")
        else:
            dq['barcode_base'] = dq[bc_col].astype(str).str.replace(r"-\d+$", "", regex=True)
            obs_base = pd.Index(adata.obs_names.str.replace(r"-\d+$", "", regex=True), name='barcode_base')
            dq = dq.set_index('barcode_base')
            # nuclear_fraction key
            if 'nuclear_fraction' in dq.columns:
                adata.obs['nuclear_fraction'] = pd.Series(index=adata.obs_names, dtype=float)
                adata.obs.loc[:, 'nuclear_fraction'] = dq.reindex(obs_base)['nuclear_fraction'].values
            else:
                print(f"Warning: 'nuclear_fraction' not found in {dq_csv}")
    else:
        print(f"Note: DropletQC metrics not found for {sample}")

    # Attach scDblFinder metadata if present
    sd_csv = sample_dir / f"{sample}_scdbl_metrics.csv"
    if sd_csv.exists():
        sd = pd.read_csv(sd_csv)
        bc_col = None
        for c in ['cell_barcode', 'barcode', 'Barcode', 'cell', 'cell_id']:
            if c in sd.columns:
                bc_col = c
                break
        if bc_col is None:
            print(f"Warning: No barcode column found in {sd_csv}")
        else:
            sd['barcode_base'] = sd[bc_col].astype(str).str.replace(r"-\d+$", "", regex=True)
            obs_base = pd.Index(adata.obs_names.str.replace(r"-\d+$", "", regex=True), name='barcode_base')
            sd = sd.set_index('barcode_base')
            if 'doublet_score' in sd.columns:
                adata.obs['doublet_score'] = sd.reindex(obs_base)['doublet_score'].values
            if 'doublet_class' in sd.columns:
                adata.obs['doublet_class'] = sd.reindex(obs_base)['doublet_class'].astype('category').values
    else:
        print(f"Note: scDblFinder metrics not found for {sample}")

    # Filtering
    keep = np.ones(adata.n_obs, dtype=bool)
    keep &= adata.obs['n_genes_by_counts'].values >= min_features
    keep &= adata.obs['total_counts'].values >= min_counts
    if np.isfinite(adata.obs['percent_mt']).any():
        keep &= (adata.obs['percent_mt'].fillna(0).values <= max_mito)
    if 'nuclear_fraction' in adata.obs.columns:
        keep &= (pd.to_numeric(adata.obs['nuclear_fraction'], errors='coerce').fillna(0).values >= min_nuclear)
    if 'doublet_class' in adata.obs.columns:
        dc = pd.Series(adata.obs['doublet_class']).astype(str).str.lower().values
        keep &= (dc != 'doublet')

    filtered = adata[keep].copy()

    return adata, filtered


def main():
    ap = argparse.ArgumentParser(description="Build per-sample postQC h5ad files from scQC-flow outputs")
    ap.add_argument('--healthy-root', type=Path, required=True, help='Path to healthy output_latest directory')
    ap.add_argument('--diseased-root', type=Path, required=True, help='Path to diseased output_latest directory')
    ap.add_argument('--out-root', type=Path, required=True, help='Base output dir for scvi_integration (will write inputs/h5ad_files/...)')
    ap.add_argument('--clear-inputs', action='store_true', help='Clear out-root/inputs before writing')
    ap.add_argument('--min-features', type=int, default=200)
    ap.add_argument('--min-counts', type=int, default=500)
    ap.add_argument('--max-mito', type=float, default=10.0)
    ap.add_argument('--min-nuclear', type=float, default=0.4)
    args = ap.parse_args()

    inputs_dir = args.out_root / 'inputs'
    h5ad_healthy_dir = inputs_dir / 'h5ad_files' / 'healthy'
    h5ad_diseased_dir = inputs_dir / 'h5ad_files' / 'diseased'

    if args.clear_inputs and inputs_dir.exists():
        # Be cautious: only remove files under inputs
        for p in inputs_dir.glob('**/*'):
            try:
                if p.is_file():
                    p.unlink()
            except Exception:
                pass
    h5ad_healthy_dir.mkdir(parents=True, exist_ok=True)
    h5ad_diseased_dir.mkdir(parents=True, exist_ok=True)

    summary = []

    def process_condition(root: Path, condition: str, out_dir: Path):
        sample_dirs = sorted([p for p in root.iterdir() if p.is_dir()])
        for sd in sample_dirs:
            sample = sd.name
            print(f"\n=== Processing {condition} sample: {sample} ===")
            try:
                ad = load_counts_from_sample_dir(sd, sample)
            except Exception as e:
                print(f"ERROR loading counts for {sample}: {e}")
                summary.append({
                    'sample': sample, 'condition': condition, 'status': 'load_error', 'error': str(e)
                })
                continue

            # store raw counts as layer
            if ad.X is None:
                print(f"ERROR: {sample} has empty matrix")
                summary.append({'sample': sample, 'condition': condition, 'status': 'empty_matrix'})
                continue
            ad.layers['counts'] = ad.X.copy()
            ad.obs['sample'] = sample
            ad.obs['condition'] = condition

            ad_pre, ad_post = attach_qc_and_filter(
                ad, sd, sample,
                min_features=args.min_features, min_counts=args.min_counts,
                max_mito=args.max_mito, min_nuclear=args.min_nuclear
            )

            # Sanity: integer-like counts
            int_like = is_integer_like_matrix(ad.layers['counts'])

            # Save postQC
            out_path = out_dir / f"{sample}_postqc.h5ad"
            ad_post.write_h5ad(out_path)

            # Collect summary
            s = {
                'sample': sample,
                'condition': condition,
                'n_obs_pre': int(ad_pre.n_obs),
                'n_var': int(ad_pre.n_vars),
                'n_obs_post': int(ad_post.n_obs),
                'integer_like_counts': bool(int_like),
                'has_percent_mt': bool('percent_mt' in ad_post.obs.columns),
                'has_nuclear_fraction': bool('nuclear_fraction' in ad_post.obs.columns),
                'has_doublet_score': bool('doublet_score' in ad_post.obs.columns),
                'has_doublet_class': bool('doublet_class' in ad_post.obs.columns),
                'output': str(out_path)
            }
            summary.append(s)

    process_condition(args.healthy_root, 'healthy', h5ad_healthy_dir)
    process_condition(args.diseased_root, 'diseased', h5ad_diseased_dir)

    # Write summary
    inputs_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summary).to_csv(inputs_dir / 'build_summary.tsv', sep='\t', index=False)
    with open(inputs_dir / 'README.txt', 'w') as fh:
        fh.write(
            "Per-sample postQC h5ad files built from scQC-flow outputs.\n"
            "Counts sourced from CellBender-derived 10x H5 when available.\n"
            f"Filters: n_genes>={args.min_features}, total_counts>={args.min_counts}, "
            f"percent_mt<={args.max_mito}, nuclear_fraction>={args.min_nuclear}, exclude doublets.\n"
        )
    print(f"\nWrote summary to {inputs_dir / 'build_summary.tsv'}")


if __name__ == '__main__':
    sys.exit(main())

