"""
Figures for an evidence report.

1) Overview: dual-annotation, intron-compressed genomic locus map. Both GFFs share the
   assembly, so the old (NCBI/LOC) and new (egapxtmp) annotations are overlaid on one
   shared pseudo-genomic axis with introns squished to a fixed width (exons stay visible
   across the ~hundreds-of-kb locus). Exons are chevrons pointing in transcription
   direction, labelled Exon 1..n per gene; a gene arrow runs under each.
2) Zoom: the read-through junction at true genomic scale (long intron compressed), with
   individual bridging reads drawn as aligned blocks joined by a thin splice connector.
"""
from __future__ import annotations
import os
from typing import List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrow

from .config import Config
from . import gff_model

_COLOR_G1, _COLOR_G2, _COLOR_OTHER = "#F2A93B", "#5B9BD5", "#9E9E9E"


def _overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def _gene_color(model, g1span, g2span):
    """gene_1 locus -> orange, gene_2 locus -> blue (consistent across both annotation
    tracks via genomic overlap); any other feature (e.g. a tRNA between them) -> gray."""
    o1 = _overlap(model.start, model.end, *g1span) if g1span else 0
    o2 = _overlap(model.start, model.end, *g2span) if g2span else 0
    if o1 == 0 and o2 == 0:
        return _COLOR_OTHER
    return _COLOR_G1 if o1 >= o2 else _COLOR_G2


class CoordCompressor:
    """Piecewise-linear genomic->plot-x map: anchor (exon/read) intervals keep true bp
    width; gaps between them are squished to a fixed plot width."""
    def __init__(self, anchors: List[Tuple[int, int]], window: Tuple[int, int],
                 gap_w: float = 6.0, bp_per_unit: float = None):
        w0, w1 = window
        merged = []
        for s, e in sorted(anchors):
            s, e = max(s, w0), min(e, w1)
            if e <= s:
                continue
            if merged and s <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        if not merged:
            merged = [(w0, w1)]
        total_bp = sum(e - s for s, e in merged) or 1
        self.scale = 1.0 / (bp_per_unit or (total_bp / 100.0))  # ~100 plot units of exon
        self.segments = []  # (g0, g1, x0, x1, kind)
        x = 0.0
        cur = w0
        for s, e in merged:
            if s > cur:                                   # leading gap
                self.segments.append((cur, s, x, x + gap_w, "gap")); x += gap_w
            xe = x + (e - s) * self.scale
            self.segments.append((s, e, x, xe, "exon")); x = xe
            cur = e
        if cur < w1:
            self.segments.append((cur, w1, x, x + gap_w, "gap")); x += gap_w
        self.xmax = x

    def x(self, pos: int) -> float:
        for g0, g1, x0, x1, _ in self.segments:
            if g0 <= pos <= g1:
                frac = (pos - g0) / (g1 - g0) if g1 > g0 else 0.0
                return x0 + frac * (x1 - x0)
        return self.xmax if pos > self.segments[-1][1] else 0.0

    def xi(self, s: int, e: int) -> Tuple[float, float]:
        a, b = self.x(s), self.x(e)
        return (a, b) if a <= b else (b, a)


def _chevron(x0, x1, yc, h, strand, tipfrac=0.4):
    tip = min(tipfrac * (x1 - x0), 0.6 * h)
    if strand == "-":
        return [(x0, yc), (x0 + tip, yc - h / 2), (x1, yc - h / 2),
                (x1, yc + h / 2), (x0 + tip, yc + h / 2)]
    return [(x0, yc - h / 2), (x1 - tip, yc - h / 2), (x1, yc),
            (x1 - tip, yc + h / 2), (x0, yc + h / 2)]


def _draw_track(ax, models, comp, ytop, label, g1span=None, g2span=None,
                exon_h=0.5, lab_exons=True):
    """Draw one annotation track of gene models; returns the y used for the gene arrows."""
    for gi, m in enumerate(models):
        color = _gene_color(m, g1span, g2span)
        exs = sorted(m.exons)
        labels = m.exon_labels()
        # intron backbone
        gx0, gx1 = comp.x(m.start), comp.x(m.end)
        ax.plot([gx0, gx1], [ytop, ytop], color=color, lw=0.8, zorder=1)
        for ei, (s, e) in enumerate(exs):
            x0, x1 = comp.xi(s, e)
            if x1 - x0 < 0.6:            # ensure tiny exons stay visible
                x1 = x0 + 0.6
            ax.add_patch(Polygon(_chevron(x0, x1, ytop, exon_h, m.strand),
                                 closed=True, facecolor=color, edgecolor="black",
                                 lw=0.4, zorder=3))
            if lab_exons:
                ax.text((x0 + x1) / 2, ytop + exon_h * 0.75, labels[ei], ha="center",
                        va="bottom", fontsize=6, rotation=90, color="#333333")
        # gene arrow + name under the exons
        ya = ytop - exon_h * 1.6
        ax.add_patch(FancyArrow(gx0 if m.strand != "-" else gx1, ya,
                                (gx1 - gx0) * (1 if m.strand != "-" else -1), 0,
                                width=exon_h * 0.18, head_width=exon_h * 0.5,
                                head_length=min(2.0, 0.15 * (gx1 - gx0) + 0.5),
                                length_includes_head=True, facecolor=color,
                                edgecolor="black", lw=0.3, zorder=2))
        ax.text((gx0 + gx1) / 2, ya - exon_h * 0.8, m.gene_id, ha="center", va="top",
                fontsize=7, fontweight="bold")
    ax.text(-1.5, ytop, label, ha="right", va="center", fontsize=8, fontweight="bold")


def overview_figure(hit, cfg: Config, out_png: str) -> Optional[str]:
    if hit.chrom is None:
        return None
    chrom, w0, w1 = hit.locus_window
    pad = int(0.02 * (w1 - w0))
    w0, w1 = w0 - pad, w1 + pad
    # keep only same-strand genes as the pair of interest (drop opposite-sense clutter)
    old = [m for m in gff_model.models_in_window(cfg.run_gff, chrom, w0, w1, "NCBI")
           if m.strand == hit.strand]
    new = [m for m in gff_model.models_in_window(cfg.new_gff, chrom, w0, w1, "new")
           if m.strand == hit.strand]
    anchors = [ex for m in (old + new) for ex in m.exons]
    comp = CoordCompressor(anchors, (w0, w1), gap_w=6.0)
    g1span = (hit.gene1_model.start, hit.gene1_model.end)
    g2span = (hit.gene2_model.start, hit.gene2_model.end)

    fig, ax = plt.subplots(figsize=(max(9, comp.xmax / 9), 4.2))
    _draw_track(ax, old, comp, ytop=3.0, label="NCBI\n(LOC)", g1span=g1span, g2span=g2span)
    _draw_track(ax, new, comp, ytop=0.6, label="sheina2025\n(egapxtmp)",
                g1span=g1span, g2span=g2span)
    ax.set_xlim(-3, comp.xmax + 1)
    ax.set_ylim(-1.2, 4.4)
    ax.axis("off")
    ax.set_title(f"{hit.fused_product}\n{chrom}:{w0:,}-{w1:,}  (introns compressed)",
                 fontsize=9)
    fig.text(0.5, 0.01, "Exons to scale; introns/intergenic squished to fixed width",
             ha="center", fontsize=6, color="#666666")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_png


def _coord_ticks(ax, comp, y):
    """Genomic-coordinate labels at the start of each anchored (exon/read) segment, so a
    compressed axis stays interpretable."""
    last_x = -1e9
    for g0, g1, x0, x1, kind in comp.segments:
        if kind == "exon" and (x0 - last_x) > 4:
            ax.plot([x0, x0], [y, y + 0.12], color="#888888", lw=0.5)
            ax.text(x0, y - 0.05, f"{g0:,}", ha="center", va="top", fontsize=5,
                    rotation=90, color="#777777")
            last_x = x0


def zoom_figure(hit, cfg: Config, reads, out_png: str,
                junction: Optional[Tuple[int, int]] = None) -> Optional[str]:
    """Read-through reads on a single compressed axis anchored on the read donor/acceptor
    sites (the long gaps between sites are squished). Reads are colored by strand
    concordance with the gene orientation (informational only — the annotation is the
    reference being questioned, not a filter, so intronic ends are not penalized).
    `reads` is a list of bam_utils.BridgingRead."""
    if hit.chrom is None or not reads:
        return None
    chrom = hit.chrom
    # reads arrive sense-only (read_level excludes opposite-strand bridging reads); guard
    # defensively in case the figure is called directly.
    reads = [r for r in reads if r.tx_strand == hit.strand]
    if not reads:
        return None
    read_blocks = [b for r in reads for b in r.blocks]
    pad = 800
    # anchor on the read sites; pull in the nearest exon on each flank for context
    anchors = [(s - pad, e + pad) for s, e in read_blocks]
    w0 = min(a[0] for a in anchors) - 6000
    w1 = max(a[1] for a in anchors) + 6000
    # keep only same-strand genes (drop opposite-sense neighbours, e.g. tRNAs)
    old = [m for m in gff_model.models_in_window(cfg.run_gff, chrom, w0, w1, "NCBI")
           if m.strand == hit.strand]
    new = [m for m in gff_model.models_in_window(cfg.new_gff, chrom, w0, w1, "new")
           if m.strand == hit.strand]
    all_exons = [ex for m in (old + new) for ex in m.exons]
    for (s, e) in read_blocks:                       # nearest exon within 12 kb of a read end
        near = [ex for ex in all_exons if min(ex[1], e + 12000) - max(ex[0], s - 12000) > 0]
        anchors.extend(near)
    window = (min(a[0] for a in anchors), max(a[1] for a in anchors))
    comp = CoordCompressor(anchors, window, gap_w=9.0)

    n = len(reads)
    g1span = (hit.gene1_model.start, hit.gene1_model.end)
    g2span = (hit.gene2_model.start, hit.gene2_model.end)
    fig, ax = plt.subplots(figsize=(max(10, comp.xmax / 7), 3.8 + 0.18 * n))
    ytop = 1.8 + 0.2 * n
    _draw_track(ax, old, comp, ytop=ytop + 1.7, label="NCBI", g1span=g1span, g2span=g2span,
                exon_h=0.45)
    _draw_track(ax, new, comp, ytop=ytop, label="new", g1span=g1span, g2span=g2span,
                exon_h=0.45)
    _coord_ticks(ax, comp, ytop - 1.0)

    # individual reads: aligned blocks joined by splice arcs that return to the baseline
    y = ytop - 1.6
    for r in sorted(reads, key=lambda r: r.donor):
        blocks = sorted(r.blocks)
        for (s, e) in blocks:
            x0, x1 = comp.xi(s, e)
            ax.plot([x0, max(x1, x0 + 0.4)], [y, y], color="#1A237E", lw=2.8,
                    solid_capstyle="butt", zorder=4)          # aligned blocks
        for i in range(len(blocks) - 1):                      # splice arc per gap
            xe, xs = comp.x(blocks[i][1]), comp.x(blocks[i + 1][0])
            ax.plot([xe, (xe + xs) / 2, xs], [y, y + 0.16, y], color="#9E9E9E",
                    lw=0.6, zorder=3)
        y -= 0.2

    umis = len({(r.barcode, r.umi) for r in reads if r.barcode and r.umi})
    ax.set_title(f"Bridging reads across the gene1–gene2 region — {len(reads)} reads, "
                 f"{umis} distinct UMIs  ({chrom}; introns compressed)", fontsize=9)
    ax.set_xlim(-3, comp.xmax + 1)
    ax.set_ylim(y - 0.5, ytop + 3.1)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_png
