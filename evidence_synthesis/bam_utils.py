"""
BAM scan for read-through (bridging) reads — shared by the read_through and
read_level evidence modules so the 19 GB BAM is fetched only once per hit.

A bridging read has an `N` (splice) CIGAR op whose donor falls in one of the pair's
genes and whose acceptor falls in the other (within tolerance). Hard gates: primary
alignment only, MAPQ >= cfg.mapq_min, NH==1 (if required), min overhang on each side.
Barcode/UMI parsed from the read name; the per-read aligned blocks are retained for
the zoomed junction figure.
"""
from __future__ import annotations
import re
import sys
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import pysam

# pysam CIGAR op codes
_BAM_CMATCH, _BAM_CREF_SKIP = 0, 3
_CONSUME_REF = {0, 2, 3, 7, 8}          # M D N = X consume reference
_ALN_BLOCK = {0, 7, 8}                   # M = X are aligned (matched) blocks


@dataclass
class BridgingRead:
    qname: str
    barcode: Optional[str]
    umi: Optional[str]
    mapq: int
    nh: Optional[int]
    is_reverse: bool
    donor: int                 # 0-based first intronic base of the gene1<->gene2 N-gap
    acceptor: int              # 0-based first exonic base after the intron
    direction: str             # "gene1->gene2" or "gene2->gene1"
    left_overhang: int
    right_overhang: int
    blocks: List[Tuple[int, int]] = field(default_factory=list)  # aligned ref blocks (0-based half-open)


def _in_span(pos: int, span: Tuple[int, int], tol: int) -> bool:
    return (span[0] - tol) <= pos <= (span[1] + tol)


def scan_bridging_reads(hit, cfg, log=sys.stdout):
    """Return (reads, stats). `reads` is a list[BridgingRead] passing all hard gates;
    `stats` reports how many reads were inspected/rejected and why."""
    bc_umi = re.compile(cfg.readname_bc_umi_regex)
    g1 = (hit.gene1_model.start, hit.gene1_model.end)   # 1-based GFF spans
    g2 = (hit.gene2_model.start, hit.gene2_model.end)
    chrom = hit.chrom
    w0, w1 = min(g1[0], g2[0]) - 1, max(g1[1], g2[1])   # fetch window (0-based start)

    stats = dict(inspected=0, has_njunction=0, kept=0,
                 rej_secondary=0, rej_mapq=0, rej_nh=0, rej_overhang=0)
    reads: List[BridgingRead] = []
    bam = pysam.AlignmentFile(cfg.bam, "rb")            # .csi auto-detected
    for r in bam.fetch(chrom, w0, w1):
        stats["inspected"] += 1
        if r.is_secondary or r.is_supplementary or r.is_unmapped:
            stats["rej_secondary"] += 1
            continue
        cigar = r.cigartuples
        if not cigar or not any(op == _BAM_CREF_SKIP for op, _ in cigar):
            continue
        # walk CIGAR, collect aligned blocks and N-gaps with adjacent overhangs
        ref = r.reference_start
        blocks: List[Tuple[int, int]] = []
        njunctions = []  # (donor, acceptor, left_overhang, right_overhang)
        prev_aln_len = 0
        for i, (op, length) in enumerate(cigar):
            if op in _ALN_BLOCK:
                blocks.append((ref, ref + length))
                prev_aln_len = length
            if op == _BAM_CREF_SKIP:
                donor, acceptor = ref, ref + length
                # right overhang = length of the next aligned block
                right = 0
                for op2, len2 in cigar[i + 1:]:
                    if op2 in _ALN_BLOCK:
                        right = len2
                        break
                njunctions.append((donor, acceptor, prev_aln_len, right))
            if op in _CONSUME_REF:
                ref += length

        # keep only N-gaps bridging gene1<->gene2
        match = None
        for donor, acceptor, lo, ro in njunctions:
            if _in_span(donor + 1, g1, cfg.junction_tol_bp) and _in_span(acceptor, g2, cfg.junction_tol_bp):
                match = (donor, acceptor, lo, ro, "gene1->gene2")
                break
            if _in_span(donor + 1, g2, cfg.junction_tol_bp) and _in_span(acceptor, g1, cfg.junction_tol_bp):
                match = (donor, acceptor, lo, ro, "gene2->gene1")
                break
        if match is None:
            continue
        stats["has_njunction"] += 1

        # hard gates
        if r.mapping_quality < cfg.mapq_min:
            stats["rej_mapq"] += 1
            continue
        nh = r.get_tag("NH") if r.has_tag("NH") else None
        if cfg.require_nh1 and nh is not None and nh != 1:
            stats["rej_nh"] += 1
            continue
        donor, acceptor, lo, ro, direction = match
        if lo < cfg.min_overhang_bp or ro < cfg.min_overhang_bp:
            stats["rej_overhang"] += 1
            continue

        m = bc_umi.search(r.query_name)
        barcode, umi = (m.group(1), m.group(2)) if m else (None, None)
        reads.append(BridgingRead(qname=r.query_name, barcode=barcode, umi=umi,
                                  mapq=r.mapping_quality, nh=nh, is_reverse=r.is_reverse,
                                  donor=donor, acceptor=acceptor, direction=direction,
                                  left_overhang=lo, right_overhang=ro, blocks=blocks))
        stats["kept"] += 1
    bam.close()
    print(f"[bam] inspected={stats['inspected']} bridging_Ngap={stats['has_njunction']} "
          f"kept={stats['kept']} (rej mapq={stats['rej_mapq']} nh={stats['rej_nh']} "
          f"overhang={stats['rej_overhang']})", file=log, flush=True)
    return reads, stats
