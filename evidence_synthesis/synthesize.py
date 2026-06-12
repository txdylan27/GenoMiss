"""
Orchestrator: synthesize all evidence lines for one (or many) GenoMiss hits.

Runs the evidence modules in tier order with graceful degradation (a module returns
status='skipped' with a reason if its inputs are missing), renders the figures, and
hands everything to report.py for the Excel + SUMMARY.md bundle. Everything is logged to
a per-hit log file beside the outputs. No automated one-gene-vs-two verdict is made.
"""
from __future__ import annotations
import os
import sys

from .config import Config, AMERICANA
from . import hit_resolver, figures, report
from .evidence import (cross_species, taxonomy, read_through, read_level,
                       coexpression, bulk_expression)

# evidence lines, in tier order
LINES = ["cross_species", "taxonomy", "read_through", "read_level",
         "coexpression", "bulk_expression"]


class _Tee:
    """Write log lines to both a file and the console, flushing each time."""
    def __init__(self, *streams):
        self.streams = streams

    def write(self, s):
        for st in self.streams:
            st.write(s)
            st.flush()

    def flush(self):
        for st in self.streams:
            st.flush()


def synthesize_hit(cfg: Config, hits_csv: str, fused_id: str, name=None) -> dict:
    hit = hit_resolver.resolve(cfg, hits_csv, fused_id, name=name)
    outdir = os.path.join(cfg.outdir, hit.name)
    os.makedirs(outdir, exist_ok=True)
    logf = open(os.path.join(outdir, f"{hit.name}.log"), "w")
    log = _Tee(sys.stdout, logf)
    print(f"=== {hit.name} | {fused_id} | {hit.gene_1}+{hit.gene_2} | "
          f"score {hit.composite_score} ===", file=log, flush=True)

    results = {}
    results["cross_species"] = cross_species.run(hit, cfg, log=log)
    results["taxonomy"] = taxonomy.run(hit, cfg, results["cross_species"], log=log)
    results["read_level"] = read_level.run(hit, cfg, log=log)
    results["read_through"] = read_through.run(hit, cfg, results["read_level"], log=log)
    results["coexpression"] = coexpression.run(hit, cfg, log=log)
    results["bulk_expression"] = bulk_expression.run(hit, cfg, log=log)

    # figures
    figs = {}
    figs["overview"] = figures.overview_figure(
        hit, cfg, os.path.join(outdir, f"{hit.name}_overview.png"))
    reads = results["read_level"].get("reads") if results["read_level"].get("status") == "ok" else None
    if reads:
        figs["zoom"] = figures.zoom_figure(
            hit, cfg, reads, os.path.join(outdir, f"{hit.name}_junction_zoom.png"))
    else:
        print("[figures] no bridging reads -> skipping zoom figure", file=log, flush=True)

    report.write_summary_md(hit, results, figs,
                            os.path.join(outdir, f"{hit.name}_SUMMARY.md"))
    logf.close()
    return {"hit": hit, "results": results, "figs": figs}


def synthesize(cfg: Config, hits_csv: str, fused_ids, names=None, xlsx_name="evidence_report.xlsx"):
    """Synthesize one or more hits and build a combined workbook (index + one sheet/hit)."""
    names = names or [None] * len(fused_ids)
    bundles = []
    for fid, nm in zip(fused_ids, names):
        bundles.append(synthesize_hit(cfg, hits_csv, fid, name=nm))
    out_xlsx = os.path.join(cfg.outdir, xlsx_name)
    report.build_workbook(bundles, out_xlsx, cfg)
    print(f"\nWrote workbook -> {out_xlsx}")
    return bundles
