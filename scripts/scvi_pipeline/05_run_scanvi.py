#!/usr/bin/env python3
"""Run scANVI label transfer using a LIBD hippocampus reference and the query dataset.

This mirrors the previous manual workflow:
    * load reference h5ad (with cell type labels)
    * load query 2k-HVG dataset (counts in layers["counts"], batch=sample)
    * intersect genes, concatenate, run SCVI pretrain + SCANVI fine-tuning
    * write annotated query-only object and diagnostic plots

Run manually, e.g.:

    python scripts/scvi_pipeline/05_run_scanvi.py \
        --reference /mnt/.../ref_LIBD_hippocampus_snRNA.h5ad \
        --query /rds/.../scvi_integration/outputs_2/adata_hvg2k_TM.h5ad \
        --out-root /rds/.../scvi_integration \
        --labels-key cell.type2

"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi


def ensure_counts_layer(adata, label):
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
        print(f"[{label}] counts layer missing -> copied from X")


def make_umap(adata, color, out_path, **kwargs):
    sc.pl.umap(adata, color=color, show=False, **kwargs)
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()


def main():
    ap = argparse.ArgumentParser(description="scANVI label transfer")
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--query", type=Path, required=True)
    ap.add_argument("--out-root", type=Path, required=True)
    ap.add_argument("--out-subdir", type=str, default="scanvi")
    ap.add_argument("--labels-key", type=str, default="cell.type2")
    ap.add_argument("--unlabeled-category", type=str, default="Unknown")
    ap.add_argument("--n-layers", type=int, default=2)
    ap.add_argument("--n-latent", type=int, default=30)
    ap.add_argument("--pretrain-epochs", type=int, default=30)
    ap.add_argument("--scanvi-epochs", type=int, default=30)
    ap.add_argument("--n-samples-per-label", type=int, default=50)
    ap.add_argument("--umap-min-dist", type=float, default=0.3)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    scvi.settings.seed = args.seed
    scvi.settings.dl_num_workers = 0

    ref = sc.read_h5ad(args.reference)
    query = sc.read_h5ad(args.query)

    ensure_counts_layer(ref, "reference")
    ensure_counts_layer(query, "query")

    if args.labels_key not in ref.obs.columns:
        raise ValueError(f"Reference does not have labels key '{args.labels_key}'")

    if "sample" not in query.obs.columns:
        raise ValueError("Query AnnData must contain 'sample' in .obs")

    ref.obs["dataset"] = "reference"
    query.obs["dataset"] = "query"

    # Harmonise counts dtype to float32 CSR for scvi-tools friendliness
    for ad in (ref, query):
        if ad.layers["counts"] is not None and ad.layers["counts"].dtype != np.float32:
            ad.layers["counts"] = ad.layers["counts"].astype(np.float32)

    # Gene intersection
    common = sorted(set(ref.var_names).intersection(set(query.var_names)))
    if len(common) < 500:
        raise ValueError(f"Only {len(common)} shared genes detected. Check gene naming / HVG selection.")
    ref = ref[:, common].copy()
    query = query[:, common].copy()

    ref.obs["celltype_scanvi"] = ref.obs[args.labels_key].astype("category")
    query.obs["celltype_scanvi"] = args.unlabeled_category

    # Build batch column: keep query sample granularity, reference -> 'reference'
    ref.obs["batch"] = "reference"
    query.obs["batch"] = query.obs["sample"].astype(str)

    adata = ref.concatenate(query, batch_key="concat_source", batch_categories=["reference", "query"], index_unique=None)

    ensure_counts_layer(adata, "joint")
    if "batch" not in adata.obs.columns:
        raise RuntimeError("Batch column missing after concatenation")

    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")

    model = scvi.model.SCVI(adata, n_layers=args.n_layers, n_latent=args.n_latent)
    model.train(max_epochs=args.pretrain_epochs)

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        model,
        labels_key="celltype_scanvi",
        unlabeled_category=args.unlabeled_category,
    )
    scanvi_model.train(max_epochs=args.scanvi_epochs, n_samples_per_label=args.n_samples_per_label)

    adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()

    preds, conf = scanvi_model.predict(soft=True, return_proba=True)
    adata.obs["scanvi_pred"] = preds
    adata.obs["scanvi_confidence"] = conf.max(axis=1)

    sc.pp.neighbors(adata, use_rep="X_scANVI")
    sc.tl.umap(adata, min_dist=args.umap_min_dist)

    out_dir = Path(args.out_root) / "outputs_2" / "scanvi" / args.out_subdir
    plots_dir = out_dir / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Joint plots
    make_umap(adata, color="dataset", out_path=plots_dir / "umap_scanvi_joint_dataset.png")
    make_umap(adata, color="scanvi_pred", out_path=plots_dir / "umap_scanvi_joint_pred.png")
    make_umap(adata, color="scanvi_confidence", out_path=plots_dir / "umap_scanvi_joint_conf.png")

    # Query-only view
    query_mask = adata.obs["dataset"] == "query"
    adata_query = adata[query_mask].copy()
    sc.pp.neighbors(adata_query, use_rep="X_scANVI")
    sc.tl.umap(adata_query, min_dist=args.umap_min_dist)

    make_umap(adata_query, color="scanvi_pred", out_path=plots_dir / "umap_scanvi_query_pred.png")
    make_umap(adata_query, color="scanvi_confidence", out_path=plots_dir / "umap_scanvi_query_conf.png")
    if "sample" in adata_query.obs:
        make_umap(adata_query, color="sample", out_path=plots_dir / "umap_scanvi_query_sample.png")
    if "condition" in adata_query.obs:
        make_umap(adata_query, color="condition", out_path=plots_dir / "umap_scanvi_query_condition.png")

    # Confidence histogram
    plt.figure(figsize=(4, 3))
    adata_query.obs["scanvi_confidence"].hist(bins=50)
    plt.xlabel("scANVI confidence")
    plt.ylabel("Cell count")
    plt.tight_layout()
    plt.savefig(plots_dir / "conf_hist_query.png", dpi=200)
    plt.close()

    # Confidence summary CSV (query only)
    conf_summary = (
        adata_query.obs.groupby("scanvi_pred")["scanvi_confidence"]
        .agg(["count", "mean", "median", "min", "max"])
        .reset_index()
    )
    conf_summary.to_csv(out_dir / "scanvi_confidence_by_label.csv", index=False)

    adata.write_h5ad(out_dir / "scanvi_joint.h5ad")
    adata_query.write_h5ad(out_dir / "query_scanvi_annotated.h5ad")

    print(f"Joint scANVI object saved to {out_dir / 'scanvi_joint.h5ad'}")
    print(f"Query annotations saved to {out_dir / 'query_scanvi_annotated.h5ad'}")
    print(f"Plots written under {plots_dir}")


if __name__ == "__main__":
    raise SystemExit(main())
