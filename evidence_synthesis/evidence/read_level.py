"""
Tier 2 evidence — read-level uniqueness of the bridging reads.

Confirms the gene1<->gene2 bridging reads are real, independent molecules: each passes
the hard gates (primary, MAPQ>=min, NH==1, min overhang) and we deduplicate by distinct
(barcode, UMI). Reports raw bridging-read count and distinct-UMI count, and retains the
per-read alignment blocks for the zoomed junction figure. (Single-cell MALBAC data is
chimera-prone, so independent UMIs matter.)
"""
from __future__ import annotations
import os
import sys

import pandas as pd

from ..config import Config
from .. import bam_utils


def run(hit, cfg: Config, log=sys.stdout) -> dict:
    if not cfg.bam or not os.path.exists(cfg.bam):
        return {"status": "skipped", "reason": "no BAM provided"}
    if hit.gene1_model is None or hit.gene2_model is None:
        return {"status": "skipped", "reason": "gene models unavailable (GFF lookup failed)"}

    reads, stats = bam_utils.scan_bridging_reads(hit, cfg, log=log)
    umis = {(r.barcode, r.umi) for r in reads if r.barcode and r.umi}
    barcodes = {r.barcode for r in reads if r.barcode}

    # concordance uses the library-corrected TRANSCRIPT strand, not the read flag strand
    n_sense = sum(r.tx_strand == hit.strand for r in reads)
    rows = [dict(qname=r.qname, barcode=r.barcode, umi=r.umi, mapq=r.mapq, nh=r.nh,
                 read_flag_strand=("-" if r.is_reverse else "+"), tx_strand=r.tx_strand,
                 gene_strand=hit.strand, strand_concordant=(r.tx_strand == hit.strand),
                 direction=r.direction, donor=r.donor, acceptor=r.acceptor,
                 left_overhang=r.left_overhang, right_overhang=r.right_overhang)
            for r in reads]
    table = pd.DataFrame(rows)
    outdir = os.path.join(cfg.outdir, hit.name)
    os.makedirs(outdir, exist_ok=True)
    table.to_csv(os.path.join(outdir, f"{hit.name}_bridging_reads.tsv"), sep="\t", index=False)

    print(f"[read_level] bridging reads={len(reads)} distinct_UMIs={len(umis)} "
          f"distinct_barcodes={len(barcodes)}", file=log, flush=True)
    return {
        "status": "ok",
        "reads": reads,                # consumed by the figure + read_through cross-check
        "tables": {"bridging_reads": table},
        "summary": {
            "n_bridging_reads": len(reads),
            "n_strand_concordant": n_sense,
            "n_strand_discordant": len(reads) - n_sense,
            "n_distinct_umis": len(umis),
            "n_distinct_barcodes": len(barcodes),
            "scan_inspected": stats["inspected"],
            "scan_rejected_mapq": stats["rej_mapq"],
            "scan_rejected_nh": stats["rej_nh"],
            "scan_rejected_overhang": stats["rej_overhang"],
        },
    }
