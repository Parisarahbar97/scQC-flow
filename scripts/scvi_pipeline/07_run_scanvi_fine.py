#!/usr/bin/env python3
"""Run scANVI label transfer using fine.cell.class labels."""

from pathlib import Path
import subprocess
import sys


def main() -> None:
    repo_scripts = Path(__file__).resolve().parent
    scvi_root = Path("/rds/general/user/pr422/projects/puklandmarkproject/live/Users/Parisa/EPILEP/scvi_integration")
    reference = scvi_root / "outputs_2/scanvi/scanvi_ref/ref_LIBD_hippocampus_snRNA.h5ad"
    query = scvi_root / "outputs_2/adata_hvg2k_TM.h5ad"

    if not reference.exists():
        raise FileNotFoundError(f"Reference file not found: {reference}")
    if not query.exists():
        raise FileNotFoundError(f"Query file not found: {query}")

    cmd = [
        sys.executable,
        str(repo_scripts / "05_run_scanvi.py"),
        "--reference",
        str(reference),
        "--query",
        str(query),
        "--out-root",
        str(scvi_root),
        "--out-subdir",
        "scanvi_31oct_fine",
        "--labels-key",
        "fine.cell.class",
        "--pretrain-epochs",
        "30",
        "--scanvi-epochs",
        "30",
        "--n-samples-per-label",
        "50",
        "--n-layers",
        "2",
        "--n-latent",
        "30",
        "--umap-min-dist",
        "0.3",
    ]

    subprocess.check_call(cmd)


if __name__ == "__main__":
    main()

