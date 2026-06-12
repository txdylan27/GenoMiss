"""
Tier 2 evidence — RNA-seq read-through junctions (BAM-derived).

The bridging reads found in the BAM scan collapse onto distinct splice coordinates
(donor->acceptor). This module reports the junction-level view: how many distinct
gene1<->gene2 junctions there are and how the reads/UMIs concentrate on the top one —
many reads on one junction is a clean read-through; reads scattered across several is
artifact-like. (STAR's SJ.out.tab is intentionally not used: it is a filtered/collapsed
junction list and drops exactly these long-intron, non-canonical, low-count junctions
even though the read alignments remain in the BAM.)
"""
from __future__ import annotations
import collections
import os
import sys

import pandas as pd

from ..config import Config


def run(hit, cfg: Config, read_level_result: dict = None, log=sys.stdout) -> dict:
    if not read_level_result or read_level_result.get("status") != "ok":
        return {"status": "skipped", "reason": "no bridging-read scan available"}
    reads = read_level_result.get("reads", [])
    if not reads:
        return {"status": "skipped", "reason": "no gene1<->gene2 bridging reads in the BAM"}

    # collapse bridging reads onto distinct (donor, acceptor) junctions; tally reads/UMIs
    tally = collections.Counter()
    umis = collections.defaultdict(set)
    for r in reads:
        key = (r.donor + 1, r.acceptor)   # 1-based intron start, 0-based first exonic base after
        tally[key] += 1
        if r.barcode and r.umi:
            umis[key].add((r.barcode, r.umi))
    rows = [dict(intron_start=k[0], intron_end=k[1], reads=n, distinct_umis=len(umis[k]))
            for k, n in tally.items()]
    table = pd.DataFrame(rows).sort_values("reads", ascending=False).reset_index(drop=True)

    outdir = os.path.join(cfg.outdir, hit.name)
    os.makedirs(outdir, exist_ok=True)
    table.to_csv(os.path.join(outdir, f"{hit.name}_bridging_junctions.csv"), index=False)

    top = table.iloc[0]
    total = len(reads)
    print(f"[read_through] distinct junctions={len(table)} | top "
          f"{int(top.intron_start)}-{int(top.intron_end)} reads={int(top.reads)} "
          f"umis={int(top.distinct_umis)} ({int(top.reads)}/{total} reads)", file=log, flush=True)
    return {
        "status": "ok",
        "tables": {"bridging_junctions": table},
        "summary": {
            "n_junctions": len(table),
            "top_junction": f"{int(top.intron_start)}-{int(top.intron_end)}",
            "top_reads": int(top.reads),
            "top_umis": int(top.distinct_umis),
            "reads_at_top": int(top.reads),
            "total_bridging_reads": total,
        },
    }
