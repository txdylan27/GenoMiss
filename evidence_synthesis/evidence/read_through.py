"""
Tier 2 evidence — RNA-seq read-through junction.

A splice junction whose intron starts in one of the pair's genes and ends in the other
is direct evidence of a single transcript spanning the locus. Junctions are taken from
STAR's SJ.out.tab (motif, novelty, unique-read support, overhang) and cross-checked
against the gene1<->gene2 junctions found in the BAM scan. Canonical motif, strand
match, and novelty are reported as annotations (not used to drop anything).
"""
from __future__ import annotations
import os
import sys

import pandas as pd

from ..config import Config

# STAR SJ.out.tab fields
_SJ_COLS = ["chrom", "intron_start", "intron_end", "strand_code", "motif_code",
            "annotated", "unique_reads", "multi_reads", "max_overhang"]
_MOTIF = {0: "non-canonical", 1: "GT/AG", 2: "CT/AC", 3: "GC/AG", 4: "CT/GC",
          5: "AT/AC", 6: "GT/AT"}
_STRAND = {0: ".", 1: "+", 2: "-"}


def _in(pos, span, tol):
    return (span[0] - tol) <= pos <= (span[1] + tol)


def run(hit, cfg: Config, read_level_result: dict = None, log=sys.stdout) -> dict:
    if hit.gene1_model is None or hit.gene2_model is None:
        return {"status": "skipped", "reason": "gene models unavailable (GFF lookup failed)"}
    g1 = (hit.gene1_model.start, hit.gene1_model.end)
    g2 = (hit.gene2_model.start, hit.gene2_model.end)
    tol = cfg.junction_tol_bp

    sj_rows = pd.DataFrame(columns=_SJ_COLS)
    if cfg.sj_tab and os.path.exists(cfg.sj_tab):
        sj = pd.read_csv(cfg.sj_tab, sep="\t", names=_SJ_COLS)
        sj = sj[sj["chrom"] == hit.chrom]
        keep = sj.apply(
            lambda r: (_in(r.intron_start, g1, tol) and _in(r.intron_end, g2, tol)) or
                      (_in(r.intron_start, g2, tol) and _in(r.intron_end, g1, tol)), axis=1)
        sj_rows = sj[keep].copy()
        if not sj_rows.empty:
            sj_rows["motif"] = sj_rows["motif_code"].map(_MOTIF)
            sj_rows["junction_strand"] = sj_rows["strand_code"].map(_STRAND)
            sj_rows["canonical"] = sj_rows["motif_code"].isin([1, 2])
            sj_rows["strand_matches_gene"] = sj_rows["junction_strand"] == hit.strand
            sj_rows["novel"] = sj_rows["annotated"] == 0
    elif cfg.sj_tab:
        print(f"[read_through] SJ.out.tab not found: {cfg.sj_tab}", file=log, flush=True)

    # BAM-derived junctions (donor+1 == SJ intron_start; acceptor == SJ intron_end)
    bam_junctions = pd.DataFrame()
    if read_level_result and read_level_result.get("status") == "ok":
        import collections
        tally = collections.Counter()
        umi = collections.defaultdict(set)
        for r in read_level_result["reads"]:
            key = (r.donor + 1, r.acceptor)
            tally[key] += 1
            if r.barcode and r.umi:
                umi[key].add((r.barcode, r.umi))
        bam_junctions = pd.DataFrame(
            [dict(intron_start=k[0], intron_end=k[1], bam_reads=n, bam_distinct_umis=len(umi[k]))
             for k, n in tally.items()])

    outdir = os.path.join(cfg.outdir, hit.name)
    os.makedirs(outdir, exist_ok=True)
    if not sj_rows.empty:
        sj_rows.to_csv(os.path.join(outdir, f"{hit.name}_sj_junctions.csv"), index=False)

    n_canonical = int(sj_rows["canonical"].sum()) if "canonical" in sj_rows else 0
    summary = {
        "n_sj_junctions": len(sj_rows),
        "n_canonical": n_canonical,
        "max_unique_reads_sj": int(sj_rows["unique_reads"].max()) if not sj_rows.empty else 0,
        "n_bam_junctions": len(bam_junctions),
    }
    if not sj_rows.empty:
        best = sj_rows.sort_values("unique_reads", ascending=False).iloc[0]
        summary.update(top_junction=f"{int(best.intron_start)}-{int(best.intron_end)}",
                       top_motif=best.motif, top_canonical=bool(best.canonical),
                       top_strand_match=bool(best.strand_matches_gene),
                       top_novel=bool(best.novel), top_unique_reads=int(best.unique_reads))
    print(f"[read_through] SJ junctions={len(sj_rows)} (canonical={n_canonical}) "
          f"BAM junctions={len(bam_junctions)}", file=log, flush=True)
    status = "ok" if (not sj_rows.empty or not bam_junctions.empty) else "skipped"
    out = {"status": status, "summary": summary,
           "tables": {"sj_junctions": sj_rows, "bam_junctions": bam_junctions}}
    if status == "skipped":
        out["reason"] = "no gene1<->gene2 junction found in SJ.out.tab or BAM"
    return out
