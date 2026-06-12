"""
Tier 3 evidence — bulk RNA-seq expression.

Reports the two genes' expression across bulk samples and their correlation (a single
transcript should track together). Expects a matrix with a gene-id column and one column
per sample; bulk uses LOC ids, so no id-map is needed for the americana run. Skips
gracefully when no bulk matrix is configured.
"""
from __future__ import annotations
import os
import sys

import numpy as np
import pandas as pd

from ..config import Config


def _load_id_map(path):
    m = {}
    if path and os.path.exists(path):
        for line in open(path):
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                m[parts[0]] = parts[1]
    return m


def run(hit, cfg: Config, log=sys.stdout) -> dict:
    if not cfg.bulk_matrix or not os.path.exists(cfg.bulk_matrix):
        return {"status": "skipped", "reason": "no bulk matrix provided"}

    sep = "\t" if cfg.bulk_matrix.endswith((".tsv", ".txt", ".gz")) else ","
    mat = pd.read_csv(cfg.bulk_matrix, sep=sep, index_col=0)
    id_map = _load_id_map(cfg.bulk_id_map)

    def find(loc):
        for cand in (id_map.get(loc, loc), loc):
            if cand in mat.index:
                return cand
        return None
    g1, g2 = find(hit.gene_1), find(hit.gene_2)
    if g1 is None or g2 is None:
        missing = [g for g, v in ((hit.gene_1, g1), (hit.gene_2, g2)) if v is None]
        return {"status": "skipped", "reason": f"gene(s) not in bulk matrix: {missing}"}

    # restrict to the configured sample set (default: head only, to avoid the
    # head/thorax tissue confound that inflates the pooled correlation)
    samples = mat.columns.tolist()
    sample_set = cfg.bulk_sample_filter or "all"
    if cfg.bulk_sample_filter:
        samples = [c for c in samples if cfg.bulk_sample_filter.lower() in str(c).lower()]
    if not samples:
        return {"status": "skipped",
                "reason": f"no bulk samples match filter {cfg.bulk_sample_filter!r}"}
    v1 = mat.loc[g1, samples].astype(float)
    v2 = mat.loc[g2, samples].astype(float)
    # Pearson on log1p to damp count-scale skew
    l1, l2 = np.log1p(v1.values), np.log1p(v2.values)
    corr = float(np.corrcoef(l1, l2)[0, 1]) if l1.std() and l2.std() else float("nan")

    table = pd.DataFrame({"sample": samples, "g1": v1.values, "g2": v2.values})
    outdir = os.path.join(cfg.outdir, hit.name)
    os.makedirs(outdir, exist_ok=True)
    table.to_csv(os.path.join(outdir, f"{hit.name}_bulk_expression.csv"), index=False)

    print(f"[bulk] {g1}/{g2} [{sample_set}]: mean {v1.mean():.1f}/{v2.mean():.1f}, "
          f"log1p Pearson r={corr:.3f} over {len(samples)} {sample_set} samples",
          file=log, flush=True)
    return {"status": "ok",
            "summary": {"bulk_g1": g1, "bulk_g2": g2, "sample_set": sample_set,
                        "n_samples": len(samples),
                        "mean_g1": round(float(v1.mean()), 2), "mean_g2": round(float(v2.mean()), 2),
                        "log1p_pearson_r": round(corr, 3)},
            "tables": {"bulk_expression": table}}
