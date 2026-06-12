"""
Tier 1 evidence — cross-species fused-protein hits.

Does a single homolog in other insects span BOTH halves of the fused query? Runs
DIAMOND twice: capped (max_target_seqs=200, reproduces GenoMiss's organism_count)
and uncapped (true taxonomic breadth). Reports per-half alignment coverage
(g1_frac/g2_frac) so a one-domain hit that slipped past the lenient +/-10 aa boundary
rule is visible, and the best hit to the focal genome (e.g. Drosophila).
"""
from __future__ import annotations
import os
import sys

import pandas as pd

from ..config import Config
from .. import diamond_utils, ncbi_symbol


def _per_half_coverage(df: pd.DataFrame, g1: int, g2len: int) -> pd.DataFrame:
    """g1_frac = fraction of gene_1 covered up to the boundary; g2_frac = fraction of
    gene_2 covered past it. Low g2_frac => partial / one-domain hit."""
    out = df.copy()
    qs = out["start_of_alignment_in_query"]
    qe = out["end_of_alignment_in_query"]
    out["g1_frac"] = ((g1 - qs).clip(lower=0) / g1).round(3)
    out["g2_frac"] = ((qe - g1).clip(lower=0) / g2len).round(3)
    return out


def run(hit, cfg: Config, log=sys.stdout) -> dict:
    outdir = os.path.join(cfg.outdir, hit.name)
    os.makedirs(outdir, exist_ok=True)
    query_faa = os.path.join(outdir, f"{hit.name}_fused_query.faa")
    diamond_utils.build_fused_query(hit.p1_seq, hit.p2_seq, hit.fused_id, query_faa)

    result = {"status": "ok", "tables": {}, "summary": {}}
    breadth = {}
    for tag, mts in (("capped", cfg.max_target_capped), ("uncapped", cfg.max_target_uncapped)):
        tsv = os.path.join(outdir, f"{hit.name}_diamond_{tag}.tsv")
        diamond_utils.run_diamond(cfg, query_faa, tsv, mts, log=log)
        df = diamond_utils.read_diamond_tsv(tsv)
        filt = diamond_utils.genomiss_filter(df, hit.fused_id, hit.gene_1_len_aa, cfg)
        filt = _per_half_coverage(filt, hit.gene_1_len_aa, hit.gene_2_len_aa)
        filt = filt.sort_values("bit_score", ascending=False)
        orgs = sorted(o for o in filt["organism"].dropna().unique())
        filt.to_csv(os.path.join(outdir, f"{hit.name}_filtered_hits_{tag}.csv"), index=False)
        with open(os.path.join(outdir, f"{hit.name}_organisms_{tag}.txt"), "w") as fo:
            fo.write("\n".join(orgs) + ("\n" if orgs else ""))
        breadth[tag] = dict(raw_hits=int((df["fused_protein"] == hit.fused_id).sum()),
                            filtered_hits=len(filt), n_organisms=len(orgs),
                            organisms=orgs, hits=filt)
        print(f"[cross_species] {tag}: filtered={len(filt)} organisms={len(orgs)}",
              file=log, flush=True)

    cap = breadth["capped"]
    result["summary"] = {
        "organism_count_capped": cap["n_organisms"],
        "organism_count_genomiss": hit.organism_count,
        "reproduces_genomiss": cap["n_organisms"] == hit.organism_count,
        "organism_count_uncapped": breadth["uncapped"]["n_organisms"],
    }
    # both-half coverage from the best uncapped hit
    un = breadth["uncapped"]["hits"]
    if not un.empty:
        top = un.iloc[0]
        result["summary"].update(
            best_hit=str(top["subject_id"]), best_hit_organism=str(top.get("organism", "")),
            best_pident=float(top["percentage_of_identical_matches"]),
            best_qcov=float(top["query_coverage"]),
            best_g1_frac=float(top["g1_frac"]), best_g2_frac=float(top["g2_frac"]))

    # focal genome (D. melanogaster) best hit -> resolve gene symbol via NCBI (fallback:
    # cleaned protein title)
    foc = un[un["subject_title"].str.contains(cfg.focal_organism, na=False)] if not un.empty else un
    if foc is not None and not foc.empty:
        fr = foc.sort_values("bit_score", ascending=False).iloc[0]
        acc = str(fr["subject_id"])
        cache_path = os.path.join(cfg.outdir, "_symbol_cache.json")
        sym = (ncbi_symbol.gene_symbol(acc, cfg.entrez_email, cache_path)
               or ncbi_symbol.clean_title(fr["subject_title"]))
        result["summary"]["focal_organism"] = cfg.focal_organism
        result["summary"]["focal_symbol"] = sym
        result["summary"]["focal_best_hit"] = acc
        result["summary"]["focal_pident"] = float(fr["percentage_of_identical_matches"])
        result["summary"]["focal_g1_frac"] = float(fr["g1_frac"])
        result["summary"]["focal_g2_frac"] = float(fr["g2_frac"])
    else:
        result["summary"]["focal_best_hit"] = None
        result["summary"]["focal_symbol"] = None

    result["breadth"] = breadth
    result["tables"]["cross_species_top_hits"] = un.head(25) if not un.empty else un
    return result
