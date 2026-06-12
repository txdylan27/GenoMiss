"""
Report assembly: per-hit Excel sheet (dashboard + embedded figures), an index sheet
for batch runs, and a SUMMARY.md. Presents every evidence line and its key numbers;
makes no automated one-gene-vs-two call.
"""
from __future__ import annotations
import os
import shutil
import time

from openpyxl import Workbook
from openpyxl.drawing.image import Image as XLImage
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

_HDR = Font(bold=True, size=12)
_SUB = Font(bold=True, size=10, color="FFFFFF")
_SUBFILL = PatternFill("solid", fgColor="4472C4")
_OK = PatternFill("solid", fgColor="E2EFDA")
_SKIP = PatternFill("solid", fgColor="F2F2F2")

_LINE_TITLES = {
    "cross_species": "Cross-species fused-protein hits",
    "taxonomy": "Taxonomic breadth",
    "read_through": "RNA-seq read-through junction",
    "read_level": "Read-level uniqueness (bridging reads)",
    "coexpression": "Single-cell co-occurrence",
    "bulk_expression": "Bulk expression",
}


def _line_metrics(name, res):
    """Return a list of (key, value) headline metrics for one evidence line."""
    if res.get("status") != "ok":
        return [("status", f"skipped — {res.get('reason', '')}")]
    s = res.get("summary", {})
    if name == "cross_species":
        return [("organisms (capped)", f"{s.get('organism_count_capped')} "
                 f"(GenoMiss {s.get('organism_count_genomiss')}; "
                 f"{'reproduced' if s.get('reproduces_genomiss') else 'DIFFERS'})"),
                ("organisms (uncapped)", s.get("organism_count_uncapped")),
                ("best hit", f"{s.get('best_hit')} [{s.get('best_hit_organism')}] "
                 f"pid={s.get('best_pident')} qcov={s.get('best_qcov')}"),
                ("best hit both-half cov", f"g1={s.get('best_g1_frac')} g2={s.get('best_g2_frac')}"),
                ("focal best hit", f"{s.get('focal_best_hit')} "
                 f"(g1={s.get('focal_g1_frac')} g2={s.get('focal_g2_frac')})"
                 if s.get("focal_best_hit") else "none")]
    if name == "taxonomy":
        rk = s.get("rank", "order")
        return [(f"{rk}s with hits", s.get(f"{rk}s_with_hits")),
                ("DB universe (n species)", s.get("db_universe_n")),
                ("overall hit rate %", s.get("overall_hitrate_pct"))]
    if name == "read_through":
        return [("SJ junctions (gene1<->gene2)", s.get("n_sj_junctions")),
                ("canonical motif", s.get("n_canonical")),
                ("top junction", f"{s.get('top_junction','-')} motif={s.get('top_motif','-')} "
                 f"strand_match={s.get('top_strand_match','-')} novel={s.get('top_novel','-')} "
                 f"unique_reads={s.get('top_unique_reads','-')}"),
                ("BAM-confirmed junctions", s.get("n_bam_junctions"))]
    if name == "read_level":
        return [("bridging reads", s.get("n_bridging_reads")),
                ("distinct UMIs", s.get("n_distinct_umis")),
                ("distinct barcodes", s.get("n_distinct_barcodes")),
                ("reads inspected", s.get("scan_inspected"))]
    if name == "coexpression":
        return [("atlas genes", f"{s.get('atlas_var_g1')} / {s.get('atlas_var_g2')}"),
                ("% cells co-express", s.get("pct_cells_coexpr")),
                ("Jaccard co-expression", s.get("jaccard_coexpr"))]
    if name == "bulk_expression":
        ss = s.get("sample_set", "all")
        return [("samples", f"{s.get('n_samples')} ({ss})"),
                ("mean g1 / g2", f"{s.get('mean_g1')} / {s.get('mean_g2')}"),
                (f"log1p Pearson r ({ss})", s.get("log1p_pearson_r"))]
    return [(k, v) for k, v in s.items()]


def add_hit_sheet(wb, bundle):
    hit, results, figs = bundle["hit"], bundle["results"], bundle["figs"]
    ws = wb.create_sheet(title=hit.name[:31])
    ws["A1"] = hit.fused_product or hit.name
    ws["A1"].font = _HDR
    meta = [("fused_protein", hit.fused_id), ("genes", f"{hit.gene_1} + {hit.gene_2}"),
            ("products", f"{hit.product_1} + {hit.product_2}"),
            ("composite_score", hit.composite_score),
            ("locus", f"{hit.chrom}:{hit.locus_window[1]:,}-{hit.locus_window[2]:,} ({hit.strand})")]
    r = 2
    for k, v in meta:
        ws.cell(r, 1, k).font = Font(bold=True)
        ws.cell(r, 2, v)
        r += 1

    r += 1
    for name in ["cross_species", "taxonomy", "read_through", "read_level",
                 "coexpression", "bulk_expression"]:
        res = results.get(name, {"status": "skipped", "reason": "not run"})
        c = ws.cell(r, 1, _LINE_TITLES[name])
        c.font = _SUB
        c.fill = _SUBFILL
        ws.cell(r, 2, "").fill = _SUBFILL
        ok = res.get("status") == "ok"
        r += 1
        for k, v in _line_metrics(name, res):
            ws.cell(r, 1, k).alignment = Alignment(indent=1)
            ws.cell(r, 2, v)
            ws.cell(r, 1).fill = _OK if ok else _SKIP
            ws.cell(r, 2).fill = _OK if ok else _SKIP
            r += 1
        r += 1

    ws.column_dimensions["A"].width = 30
    ws.column_dimensions["B"].width = 70

    # embed figures to the right of the dashboard
    img_col = "D"
    img_row = 2
    for key in ("overview", "zoom"):
        p = figs.get(key)
        if p and os.path.exists(p):
            img = XLImage(p)
            img.anchor = f"{img_col}{img_row}"
            ws.add_image(img)
            img_row += 34
    return ws


def build_workbook(bundles, out_xlsx):
    if os.path.exists(out_xlsx):  # backup before overwrite (user verifies formatting)
        shutil.copy2(out_xlsx, out_xlsx + f".bak_{time.strftime('%Y%m%d_%H%M%S')}")
    wb = Workbook()
    idx = wb.active
    idx.title = "Index"
    headers = ["hit", "genes", "score", "organisms(capped)", "reproduces_GenoMiss",
               "orders_hit", "bridging_reads", "distinct_UMIs", "%coexpr"]
    for j, h in enumerate(headers, 1):
        c = idx.cell(1, j, h)
        c.font = _SUB
        c.fill = _SUBFILL
    for i, b in enumerate(bundles, 2):
        hit, res = b["hit"], b["results"]
        cs = res["cross_species"].get("summary", {})
        tx = res["taxonomy"].get("summary", {})
        rl = res["read_level"].get("summary", {})
        co = res["coexpression"].get("summary", {})
        row = [hit.name, f"{hit.gene_1}+{hit.gene_2}", hit.composite_score,
               cs.get("organism_count_capped"), cs.get("reproduces_genomiss"),
               tx.get(f"{tx.get('rank','order')}s_with_hits"),
               rl.get("n_bridging_reads"), rl.get("n_distinct_umis"),
               co.get("pct_cells_coexpr")]
        for j, v in enumerate(row, 1):
            idx.cell(i, j, v)
        add_hit_sheet(wb, b)
    for j in range(1, len(headers) + 1):
        idx.column_dimensions[get_column_letter(j)].width = 18
    wb.save(out_xlsx)
    return out_xlsx


def write_summary_md(hit, results, figs, path):
    L = []
    L.append(f"# Evidence synthesis — {hit.fused_product or hit.name}\n")
    L.append(f"**Fused protein:** `{hit.fused_id}`  ")
    L.append(f"**Genes:** {hit.gene_1} ({hit.product_1}) + {hit.gene_2} ({hit.product_2})  ")
    L.append(f"**Locus:** {hit.chrom}:{hit.locus_window[1]:,}-{hit.locus_window[2]:,} "
             f"({hit.strand})  ")
    L.append(f"**GenoMiss composite score:** {hit.composite_score}\n")
    L.append("Evidence below is presented without an automated one-gene-vs-two verdict.\n")
    for name in ["cross_species", "taxonomy", "read_through", "read_level",
                 "coexpression", "bulk_expression"]:
        res = results.get(name, {"status": "skipped", "reason": "not run"})
        L.append(f"## {_LINE_TITLES[name]}")
        if res.get("status") != "ok":
            L.append(f"_skipped — {res.get('reason','')}_\n")
            continue
        for k, v in _line_metrics(name, res):
            L.append(f"- **{k}:** {v}")
        L.append("")
    if figs.get("overview"):
        L.append(f"## Figures")
        L.append(f"- Locus overview (dual annotation): `{os.path.basename(figs['overview'])}`")
        if figs.get("zoom"):
            L.append(f"- Read-through junction zoom: `{os.path.basename(figs['zoom'])}`")
        L.append("")
    with open(path, "w") as fo:
        fo.write("\n".join(L))
    return path
