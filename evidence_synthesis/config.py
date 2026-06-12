"""
Configuration for the GenoMiss evidence-synthesis tool.

v1 targets Schistocerca americana: paths/params below are the verified americana
defaults. Every field can be overridden from the CLI (see cli.py). For now we hardcode
what the Dec-3 GenoMiss run used; americana-specific assumptions are isolated here so a
later pass can generalize them.
"""
from __future__ import annotations
import os
from dataclasses import dataclass, field, replace
from typing import Optional

# Repo root = two levels up from this file (.../Gene-Misannotation-Tool/)
TOOL = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@dataclass
class Config:
    # --- core (Tier 1: inherits GenoMiss's own -p/-a/-db/-d inputs) ---
    tool_root: str = TOOL
    diamond_bin: str = os.path.expanduser("~/programs/anaconda3/envs/genomiss/bin/diamond")
    diamond_db: str = f"{TOOL}/updated_insecta_refseq_wtaxid.dmnd"
    proteome: str = f"{TOOL}/input_files/americana.faa"
    # GFF GenoMiss actually ran on (LOC ids) -> drives gene-ID linkage + "old" annotation track
    run_gff: str = f"{TOOL}/input_files/schistocercaAmericana.gff"
    # newer, more accurate annotation (egapxtmp ids) -> "new" annotation track; SAME assembly
    new_gff: str = f"{TOOL}/input_files/sheina2025_schistocercaAmericana.gff"
    taxon_exclude: str = "7009"            # S. americana (self); from run_gff ##species
    focal_organism: str = "Drosophila melanogaster"
    entrez_email: str = "david.bellini@bcm.edu"   # required by NCBI E-utilities

    # DIAMOND breadth: capped reproduces GenoMiss exactly; uncapped = true taxonomic breadth
    max_target_capped: int = 200           # GenoMiss.py:871
    max_target_uncapped: int = 25000
    diamond_threads: int = 8
    diamond_evalue: str = "1e-5"
    boundary_slack_aa: int = 10            # GenoMiss ±10 aa boundary rule
    min_qcov: float = 50.0
    min_pident: float = 50.0
    taxonomic_rank: str = "order"         # ete3 rank for the breadth table

    # --- Tier 2 (sequencing): optional, enables read-through + read-level + zoom figure ---
    bam: Optional[str] = (
        "/home/davidbellini/files/workspace/americana_atlas3/mapping_pipeline/"
        "work/f5/26df0fd8196f970ce1f36809d47686/Aligned.sortedByCoord.out.bam"
    )
    sj_tab: Optional[str] = (
        "/home/davidbellini/files/workspace/americana_atlas3/mapping_pipeline/"
        "work/f5/26df0fd8196f970ce1f36809d47686/SJ.out.tab"
    )
    # bridging-read gates (hard filters)
    mapq_min: int = 250                    # STAR unique == 255; >=250 per user
    require_nh1: bool = True
    min_overhang_bp: int = 8               # min aligned bp flanking the splice on each side
    junction_tol_bp: int = 10              # tolerance matching an N-gap to the gene1->gene2 intron
    # featureCounts -s strandedness of the library (mapping_pipeline nextflow.config: 2).
    # Single-end + reverse-stranded => transcript strand is OPPOSITE the read flag strand
    # (this is the source of IGV's apparent strand flip). 1=forward, 2=reverse, 0=unstranded.
    library_strandedness: int = 2
    # barcode/UMI parsed from read name: <illumina_id>_<barcode>_<UMI>
    readname_bc_umi_regex: str = r"_([ACGTN]+)_([ACGTN]+)$"

    # --- Tier 3 (expression): optional, needs atlas/bulk + id-map ---
    atlas_h5ad: Optional[str] = (
        "/home/davidbellini/files/workspace/americana_atlas3/R_analysis/sa3_labeled_apr14.h5ad"
    )
    # gene-ID map for the atlas namespace (egapxtmp) <-> LOC; None -> try LOC ids directly
    atlas_id_map: Optional[str] = None
    # Maeva bulk counts (LOC ids; 20 samples AME G/S x head/thorax) — no id-map needed
    bulk_matrix: Optional[str] = ("/home/davidbellini/OneDrive/gabbiani/PROJECTS/"
                                  "genomiss_validations/supporting/americana_bulk_counts_maeva.csv")
    bulk_id_map: Optional[str] = None
    # substring selecting which bulk sample columns to correlate over (None = all).
    # "head" -> head samples only, avoiding the head/thorax tissue confound in the pooled r.
    bulk_sample_filter: Optional[str] = "head"

    # --- output ---
    outdir: str = f"{TOOL}/evidence_synthesis/output"
    # cached DB organism universe (denominator); derived from diamond_db if missing
    db_universe: str = f"{TOOL}/evidence_synthesis/output/db_all_organisms.txt"

    def with_overrides(self, **kw) -> "Config":
        """Return a copy with non-None overrides applied."""
        clean = {k: v for k, v in kw.items() if v is not None}
        return replace(self, **clean)


# DIAMOND outfmt-6 field order (matches GenoMiss.py:214-218) and the column names we parse them into.
DIAMOND_OUTFMT = [
    "qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "sstart", "send",
    "pident", "nident", "mismatch", "evalue", "bitscore", "length",
    "qcovhsp", "scovhsp", "qtitle", "stitle",
]
DIAMOND_COLUMNS = [
    "fused_protein", "query_length", "subject_id", "subject_length",
    "start_of_alignment_in_query", "end_of_alignment_in_query",
    "start_of_alignment_in_subject", "end_of_alignment_in_subject",
    "percentage_of_identical_matches", "number_of_identical_matches",
    "number_of_mismatches", "expected_value", "bit_score", "alignment_length",
    "query_coverage", "subject_coverage", "query_title", "subject_title",
]

AMERICANA = Config()
