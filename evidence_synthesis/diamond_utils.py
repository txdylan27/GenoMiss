"""
Shared DIAMOND helpers — consolidates the logic copy-pasted across
recover_organisms.py / run_uncapped.py / analyze_pair.py into one place.

Faithfully reproduces the GenoMiss fused-hit invocation and filter so that
organism_count is recovered exactly (GenoMiss.py:208-263, max_target_seqs at :871).
"""
from __future__ import annotations
import os
import subprocess
import sys
import pandas as pd

from .config import Config, DIAMOND_OUTFMT, DIAMOND_COLUMNS


def read_fasta(path: str) -> dict:
    """Minimal FASTA reader keyed by the first whitespace token of the header."""
    seqs, cur, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def build_fused_query(p1_seq: str, p2_seq: str, fused_id: str, out_faa: str) -> int:
    """Write protein_1 + protein_2 concatenated (no separator -> matches GenoMiss fusion).
    Returns the gene_1 length (aa boundary used by the GenoMiss filter)."""
    os.makedirs(os.path.dirname(out_faa), exist_ok=True)
    with open(out_faa, "w") as out:
        out.write(f">{fused_id}\n{p1_seq}{p2_seq}\n")
    return len(p1_seq)


def run_diamond(cfg: Config, query_faa: str, out_tsv: str, max_target_seqs: int,
                log=sys.stdout) -> str:
    """Run diamond blastp with GenoMiss's exact flags. Mirrors GenoMiss's graceful
    fallback when the DB lacks taxonomy info for --taxon-exclude (GenoMiss.py:188-191)."""
    base = [cfg.diamond_bin, "blastp", "--db", cfg.diamond_db, "--query", query_faa,
            "--out", out_tsv, "--outfmt", "6", *DIAMOND_OUTFMT,
            "--header", "--evalue", cfg.diamond_evalue,
            "--threads", str(cfg.diamond_threads),
            "--max-target-seqs", str(max_target_seqs)]
    cmd = base + (["--taxon-exclude", cfg.taxon_exclude] if cfg.taxon_exclude else [])
    print(f"[diamond] max_target_seqs={max_target_seqs} taxon_exclude={cfg.taxon_exclude}",
          file=log, flush=True)
    r = subprocess.run(cmd, text=True, capture_output=True)
    if r.returncode != 0:
        if cfg.taxon_exclude and "taxonomy" in (r.stderr or "").lower():
            print("[diamond] DB lacks taxonomy tree for --taxon-exclude; retrying without it",
                  file=log, flush=True)
            r = subprocess.run(base, text=True, capture_output=True)
        if r.returncode != 0:
            sys.stderr.write(r.stderr or "")
            r.check_returncode()
    if r.stderr:
        print(r.stderr, file=log, flush=True)
    return out_tsv


def read_diamond_tsv(tsv: str) -> pd.DataFrame:
    """Parse a --header DIAMOND tsv into the GenoMiss column names (skiprows=3)."""
    return pd.read_csv(tsv, sep="\t", skiprows=3, names=DIAMOND_COLUMNS)


def extract_organisms(df: pd.DataFrame) -> pd.Series:
    """Organism name from the trailing [Organism] bracket — identical to GenoMiss.py:634."""
    return df["subject_title"].str.extract(r"\[([^\]]+)\]$")[0]


def genomiss_filter(df: pd.DataFrame, fused_id: str, gene_1_len: int, cfg: Config) -> pd.DataFrame:
    """The GenoMiss fused-hit filter (GenoMiss.py:257-263):
       alignment crosses the gene1/gene2 boundary by >= slack aa on each side,
       qcov & pident above cutoffs, and subject not titled 'uncharacterized'."""
    sub = df[df["fused_protein"] == fused_id].copy()
    g1, slack = gene_1_len, cfg.boundary_slack_aa
    filt = sub[
        (sub["start_of_alignment_in_query"] < (g1 - slack)) &
        (sub["end_of_alignment_in_query"] > (g1 + slack)) &
        (sub["query_coverage"] >= cfg.min_qcov) &
        (sub["percentage_of_identical_matches"] >= cfg.min_pident) &
        (~sub["subject_title"].str.contains("uncharacterized", case=False, na=False))
    ].copy()
    filt["organism"] = extract_organisms(filt)
    return filt
