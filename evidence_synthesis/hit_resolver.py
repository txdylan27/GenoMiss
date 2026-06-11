"""
Resolve a GenoMiss hit into a fully-specified Hit object.

Input: a fused_hits.csv + a `fused_protein` id (XP_a_XP_b). Pulls P1/P2, gene_1/gene_2
LOC ids, gene_1_len_aa/gene_2_len_aa, products, organism_count, score, intron length
straight from the row (no hand-edited PAIR block); protein sequences from the proteome;
genomic coords/strand/exons of the two genes from the run GFF.
"""
from __future__ import annotations
import re
from dataclasses import dataclass, field
from typing import List, Optional

import pandas as pd

from .config import Config
from . import diamond_utils, gff_model


@dataclass
class Hit:
    fused_id: str
    name: str                      # output slug
    p1: str
    p2: str
    gene_1: str                    # LOC id
    gene_2: str
    gene_1_len_aa: int
    gene_2_len_aa: int
    product_1: str
    product_2: str
    fused_product: str
    organism_count: int            # GenoMiss's count (expected to reproduce)
    composite_score: float
    theorized_intron_length_bp: Optional[int]
    p1_seq: str = ""
    p2_seq: str = ""
    # genomic (from run GFF); gene_1 = first by genomic coord is NOT assumed -> filled per LOC
    chrom: Optional[str] = None
    strand: Optional[str] = None
    gene1_model: Optional[gff_model.GeneModel] = None
    gene2_model: Optional[gff_model.GeneModel] = None

    @property
    def fused_seq(self) -> str:
        return self.p1_seq + self.p2_seq

    @property
    def locus_window(self):
        """(chrom, start, end) spanning both genes from the run annotation."""
        starts = [m.start for m in (self.gene1_model, self.gene2_model) if m]
        ends = [m.end for m in (self.gene1_model, self.gene2_model) if m]
        return self.chrom, (min(starts) if starts else None), (max(ends) if ends else None)


def _slug(s: str) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "_", s).strip("_")


def resolve(cfg: Config, hits_csv: str, fused_id: str,
            name: Optional[str] = None) -> Hit:
    df = pd.read_csv(hits_csv)
    rows = df[df["fused_protein"] == fused_id]
    if rows.empty:
        raise ValueError(f"fused_protein {fused_id!r} not found in {hits_csv}")
    r = rows.iloc[0]

    # pretty default name from the two products, else from gene ids
    prod = str(r.get("fused_product", "")) or fused_id
    default_name = _slug(prod)[:60] or _slug(fused_id)

    hit = Hit(
        fused_id=fused_id,
        name=name or default_name,
        p1=str(r["product_1"]), p2=str(r["product_2"]),
        gene_1=str(r["gene_1"]), gene_2=str(r["gene_2"]),
        gene_1_len_aa=int(r["gene_1_len_aa"]), gene_2_len_aa=int(r["gene_2_len_aa"]),
        product_1=str(r["product_1"]), product_2=str(r["product_2"]),
        fused_product=str(r.get("fused_product", "")),
        organism_count=int(r["organism_count"]),
        composite_score=float(r["composite_score"]),
        theorized_intron_length_bp=(int(r["theorized_intron_length_bp"])
                                    if pd.notna(r.get("theorized_intron_length_bp")) else None),
    )

    # protein sequences (validate gene_1_len matches the proteome, as recover_organisms did)
    seqs = diamond_utils.read_fasta(cfg.proteome)
    hit.p1_seq, hit.p2_seq = seqs[hit.p1], seqs[hit.p2]
    if len(hit.p1_seq) != hit.gene_1_len_aa:
        raise ValueError(f"{hit.p1} length {len(hit.p1_seq)} != gene_1_len_aa {hit.gene_1_len_aa}")

    # genomic coords + exon models from the run GFF
    found = gff_model.find_genes(cfg.run_gff, {hit.gene_1, hit.gene_2})
    if hit.gene_1 in found and hit.gene_2 in found:
        chrom = found[hit.gene_1]["chrom"]
        hit.chrom = chrom
        hit.strand = found[hit.gene_1]["strand"]
        lo = min(found[hit.gene_1]["start"], found[hit.gene_2]["start"])
        hi = max(found[hit.gene_1]["end"], found[hit.gene_2]["end"])
        models = gff_model.models_in_window(cfg.run_gff, chrom, lo, hi, source="NCBI")
        by_name = {m.gene_id: m for m in models}
        hit.gene1_model = by_name.get(found[hit.gene_1]["name"]) or by_name.get(hit.gene_1)
        hit.gene2_model = by_name.get(found[hit.gene_2]["name"]) or by_name.get(hit.gene_2)
    return hit
