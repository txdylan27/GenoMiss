"""
Retrospective entrypoint: synthesize evidence for one or more GenoMiss hits from an
existing fused_hits.csv.

  python -m evidence_synthesis.cli \
      --hits-csv annotation_comparison/inputs/fused_hits.csv \
      --fused-id XP_046994343.1_XP_046994344.1 [--fused-id ...] \
      [--name dipep1_dipep1] [--outdir ...] [--no-bam] [--no-atlas]

Any Config field can be overridden with --set field=value.
"""
from __future__ import annotations
import argparse
import sys

from .config import AMERICANA, Config
from . import synthesize


def _build_cfg(args) -> Config:
    overrides = {}
    for kv in args.set or []:
        k, v = kv.split("=", 1)
        overrides[k] = v
    if args.outdir:
        overrides["outdir"] = args.outdir
    if args.no_bam:
        overrides["bam"] = None
    if args.no_atlas:
        overrides["atlas_h5ad"] = None
    if args.diamond_db:
        overrides["diamond_db"] = args.diamond_db
    if args.proteome:
        overrides["proteome"] = args.proteome
    if args.run_gff:
        overrides["run_gff"] = args.run_gff
    if args.new_gff:
        overrides["new_gff"] = args.new_gff
    return AMERICANA.with_overrides(**overrides)


def main(argv=None):
    ap = argparse.ArgumentParser(description="GenoMiss evidence synthesis (retrospective).")
    ap.add_argument("--hits-csv", required=True, help="GenoMiss fused_hits.csv")
    ap.add_argument("--fused-id", action="append", required=True,
                    help="fused_protein id (XP_a_XP_b); repeatable")
    ap.add_argument("--name", action="append", default=None,
                    help="output slug per --fused-id (optional, same order)")
    ap.add_argument("--outdir")
    ap.add_argument("--xlsx", default="evidence_report.xlsx")
    ap.add_argument("--no-bam", action="store_true", help="skip Tier-2 sequencing evidence")
    ap.add_argument("--no-atlas", action="store_true", help="skip single-cell co-occurrence")
    ap.add_argument("--diamond-db")
    ap.add_argument("--proteome")
    ap.add_argument("--run-gff")
    ap.add_argument("--new-gff")
    ap.add_argument("--set", action="append", help="override any Config field: field=value")
    args = ap.parse_args(argv)

    cfg = _build_cfg(args)
    names = args.name if args.name else None
    if names and len(names) != len(args.fused_id):
        ap.error("--name count must match --fused-id count")
    synthesize.synthesize(cfg, args.hits_csv, args.fused_id, names=names, xlsx_name=args.xlsx)


if __name__ == "__main__":
    main(sys.argv[1:])
