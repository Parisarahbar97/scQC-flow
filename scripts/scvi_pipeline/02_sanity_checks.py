#!/usr/bin/env python3
import argparse
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def is_integer_like_matrix(mtx, check_n=200000):
    if mtx is None:
        return False
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


def check_file(p: Path):
    try:
        a = ad.read_h5ad(p)
    except Exception as e:
        return {'file': str(p), 'status': 'read_error', 'error': str(e)}

    checks = {
        'file': str(p),
        'status': 'ok',
        'n_obs': int(a.n_obs),
        'n_var': int(a.n_vars),
        'obs_names_unique': bool(a.obs_names.is_unique),
        'var_names_unique': bool(a.var_names.is_unique),
        'has_counts_layer': bool('counts' in a.layers),
        'counts_integer_like': bool(is_integer_like_matrix(a.layers['counts']) if 'counts' in a.layers else False),
        'has_percent_mt': bool('percent_mt' in a.obs.columns),
        'has_nuclear_fraction': bool('nuclear_fraction' in a.obs.columns),
        'has_doublet_score': bool('doublet_score' in a.obs.columns),
        'has_doublet_class': bool('doublet_class' in a.obs.columns),
        'has_sample': bool('sample' in a.obs.columns),
        'has_condition': bool('condition' in a.obs.columns),
    }
    return checks


def main():
    ap = argparse.ArgumentParser(description='Sanity-check all per-sample h5ad files')
    ap.add_argument('--inputs-root', type=Path, required=True, help='scvi_integration/inputs root path')
    args = ap.parse_args()

    h5ads = sorted(list((args.inputs_root / 'h5ad_files').glob('**/*_postqc.h5ad')))
    results = [check_file(p) for p in h5ads]

    if not results:
        msg = (
            "No *_postqc.h5ad files found under inputs/h5ad_files. "
            "Run 01_build_per_sample_h5ad.py first."
        )
        (args.inputs_root / 'sanity_report.txt').write_text(msg + '\n')
        print(msg)
        return

    df = pd.DataFrame(results)
    df.to_csv(args.inputs_root / 'sanity_report.tsv', sep='\t', index=False)

    # Human-readable summary
    ok = (df['status'] == 'ok').sum()
    total = len(df)
    lines = [
        f"Sanity checks for {total} files: {ok} OK, {total-ok} issues",
        "Criteria: unique obs/var names, counts layer present & integer-like, key QC columns present.",
        "",
    ]
    # Flag any failures
    issues = df[df['status'] != 'ok']
    if not issues.empty:
        lines.append('Read errors:')
        for _, r in issues.iterrows():
            lines.append(f"- {r['file']}: {r.get('error','unknown error')}")
        lines.append('')

    # Summaries
    missing = []
    for col in ['obs_names_unique', 'var_names_unique', 'has_counts_layer', 'counts_integer_like',
                'has_percent_mt', 'has_nuclear_fraction', 'has_doublet_score', 'has_doublet_class',
                'has_sample', 'has_condition']:
        bad = df[(df['status'] == 'ok') & (~df[col])]
        if not bad.empty:
            missing.append((col, bad['file'].tolist()))

    if missing:
        lines.append('Issues in OK files:')
        for col, files in missing:
            lines.append(f"- {col}: {len(files)} failing")
            for f in files[:10]:
                lines.append(f"  * {f}")
        lines.append('')

    (args.inputs_root / 'sanity_report.txt').write_text('\n'.join(lines))
    print('\n'.join(lines))


if __name__ == '__main__':
    raise SystemExit(main())
