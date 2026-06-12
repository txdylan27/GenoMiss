"""
Tier 3 evidence — single-cell co-occurrence.

If the locus is one gene, its two halves should co-express in the same cells; two genes
may differ. Loads the atlas and asks, per cell type, what fraction of cells express each
gene and what fraction express both. The atlas uses egapxtmp ids (+ finalized symbols)
while the hit uses LOC ids, so genes are mapped LOC->atlas-id by genomic overlap in the
new GFF (both GFFs share the assembly) unless an explicit atlas_id_map is supplied.
"""
from __future__ import annotations
import os
import sys

import numpy as np
import pandas as pd

from ..config import Config
from .. import gff_model

_CELLTYPE_CANDIDATES = ["cell_type", "celltype", "annotation", "seurat_clusters", "leiden"]


def _load_id_map(path):
    m = {}
    if path and os.path.exists(path):
        for line in open(path):
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                m[parts[0]] = parts[1]
    return m


def _resolve_var(loc_id, model, cfg, var_set, id_map):
    """Map a hit gene to an atlas var name: explicit map -> direct LOC -> egapxtmp gene
    overlapping the LOC locus in the new GFF -> overlapping gene's display id."""
    if loc_id in id_map and id_map[loc_id] in var_set:
        return id_map[loc_id]
    if loc_id in var_set:
        return loc_id
    if model is not None and cfg.new_gff and os.path.exists(cfg.new_gff):
        overl = gff_model.models_in_window(cfg.new_gff, model.chrom, model.start, model.end,
                                           source="new")
        # pick the overlapping gene with the largest reciprocal overlap whose id is in the
        # atlas. The new GFF uses egapxtmp_NNNNNN (underscore) while the sa3 atlas uses
        # egapxtmp-NNNNNN (dash), so try the dash-normalized form too.
        best = None
        for m in overl:
            ov = min(model.end, m.end) - max(model.start, m.start)
            if ov <= 0:
                continue
            base = m.gene_id.replace("gene-", "")
            for cand in (m.gene_id, base, base.replace("_", "-")):
                if cand in var_set and (best is None or ov > best[0]):
                    best = (ov, cand)
        if best:
            return best[1]
    return None


def run(hit, cfg: Config, log=sys.stdout) -> dict:
    if not cfg.atlas_h5ad or not os.path.exists(cfg.atlas_h5ad):
        return {"status": "skipped", "reason": "no atlas h5ad provided"}
    import anndata as ad

    adata = ad.read_h5ad(cfg.atlas_h5ad)
    var_set = set(map(str, adata.var_names))
    id_map = _load_id_map(cfg.atlas_id_map)
    v1 = _resolve_var(hit.gene_1, hit.gene1_model, cfg, var_set, id_map)
    v2 = _resolve_var(hit.gene_2, hit.gene2_model, cfg, var_set, id_map)
    if v1 is None or v2 is None:
        missing = [g for g, v in ((hit.gene_1, v1), (hit.gene_2, v2)) if v is None]
        return {"status": "skipped",
                "reason": f"gene(s) not found in atlas: {missing} (supply atlas_id_map)"}
    if v1 == v2:
        # both LOC genes overlap a single gene in the new annotation -> already merged there;
        # co-occurrence is degenerate (same gene). Report the merge, not a fake jaccard.
        print(f"[coexpression] {hit.gene_1}+{hit.gene_2} both map to {v1} "
              f"-> merged in new annotation; co-occurrence N/A", file=log, flush=True)
        return {"status": "ok", "tables": {},
                "summary": {"atlas_var_g1": v1, "atlas_var_g2": v2,
                            "merged_in_new_annotation": True,
                            "pct_cells_coexpr": None, "jaccard_coexpr": None}}

    def col(v):
        x = adata[:, v].X
        x = x.toarray().ravel() if hasattr(x, "toarray") else np.asarray(x).ravel()
        return x
    e1, e2 = col(v1), col(v2)
    b1, b2 = e1 > 0, e2 > 0
    n = adata.n_obs

    ctcol = next((c for c in _CELLTYPE_CANDIDATES if c in adata.obs.columns), None)
    per_type = pd.DataFrame()
    if ctcol:
        df = pd.DataFrame({"ct": adata.obs[ctcol].astype(str).values,
                           "b1": b1, "b2": b2, "both": b1 & b2})
        g = df.groupby("ct")
        per_type = pd.DataFrame({
            "n_cells": g.size(),
            "pct_expr_g1": (g["b1"].mean() * 100).round(1),
            "pct_expr_g2": (g["b2"].mean() * 100).round(1),
            "pct_coexpr": (g["both"].mean() * 100).round(1),
        }).reset_index().sort_values("n_cells", ascending=False)
        outdir = os.path.join(cfg.outdir, hit.name)
        os.makedirs(outdir, exist_ok=True)
        per_type.to_csv(os.path.join(outdir, f"{hit.name}_coexpression_by_celltype.csv"), index=False)

    n_either = int((b1 | b2).sum())
    summary = {
        "atlas_var_g1": v1, "atlas_var_g2": v2,
        "pct_cells_expr_g1": round(100 * b1.mean(), 2),
        "pct_cells_expr_g2": round(100 * b2.mean(), 2),
        "pct_cells_coexpr": round(100 * (b1 & b2).mean(), 2),
        "jaccard_coexpr": round(int((b1 & b2).sum()) / n_either, 3) if n_either else 0.0,
        "celltype_col": ctcol,
    }
    print(f"[coexpression] {v1}/{v2}: coexpr {summary['pct_cells_coexpr']}% of cells, "
          f"jaccard {summary['jaccard_coexpr']}", file=log, flush=True)
    return {"status": "ok", "summary": summary, "tables": {"coexpression_by_celltype": per_type}}
