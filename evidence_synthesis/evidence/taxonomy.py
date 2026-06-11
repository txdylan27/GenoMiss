"""
Tier 1 evidence — taxonomic breadth.

Builds the per-rank hit-rate table: for each taxonomic rank (default order), how many
DB species exist (denominator universe) vs how many the fused query hits. The universe
is derived ONCE from the user's own DIAMOND DB and cached (so the denominator always
matches the DB in use), rather than shipping a fixed species list.
"""
from __future__ import annotations
import os
import re
import subprocess
import sys

from ..config import Config
from .. import taxonomy_utils


def ensure_db_universe(cfg: Config, log=sys.stdout) -> str:
    """Return path to the cached DB organism universe; derive from the DB if absent."""
    if os.path.exists(cfg.db_universe) and os.path.getsize(cfg.db_universe) > 0:
        return cfg.db_universe
    # seed from the legacy file derived from this same DB, if available
    legacy = os.path.join(cfg.tool_root,
                          "netrin_dipep_organism_recovery/output/db_all_organisms.txt")
    os.makedirs(os.path.dirname(cfg.db_universe), exist_ok=True)
    if os.path.exists(legacy) and os.path.getsize(legacy) > 0:
        print(f"[taxonomy] seeding DB universe from {legacy}", file=log, flush=True)
        with open(legacy) as a, open(cfg.db_universe, "w") as b:
            b.write(a.read())
        return cfg.db_universe
    # derive from the DB: dump headers, extract trailing [Organism]
    print("[taxonomy] deriving DB organism universe from DIAMOND DB (one-time)",
          file=log, flush=True)
    p = subprocess.Popen([cfg.diamond_bin, "getseq", "--db", cfg.diamond_db],
                         stdout=subprocess.PIPE, text=True)
    orgs = set()
    for line in p.stdout:
        if line.startswith(">"):
            m = re.search(r"\[([^\]]+)\]\s*$", line.rstrip())
            if m:
                orgs.add(m.group(1))
    p.wait()
    with open(cfg.db_universe, "w") as fo:
        fo.write("\n".join(sorted(orgs)) + "\n")
    return cfg.db_universe


def run(hit, cfg: Config, cross_species_result: dict, log=sys.stdout) -> dict:
    if cross_species_result.get("status") != "ok":
        return {"status": "skipped", "reason": "cross_species evidence unavailable"}
    universe_path = ensure_db_universe(cfg, log=log)
    with open(universe_path) as fh:
        db_orgs = [l.strip() for l in fh if l.strip()]
    hit_orgs = cross_species_result["breadth"]["uncapped"]["organisms"]

    table = taxonomy_utils.rank_hitrate_table(db_orgs, hit_orgs, rank=cfg.taxonomic_rank,
                                              hit_label="hits")
    outdir = os.path.join(cfg.outdir, hit.name)
    os.makedirs(outdir, exist_ok=True)
    table.to_csv(os.path.join(outdir, f"{hit.name}_{cfg.taxonomic_rank}_hitrate.csv"), index=False)

    total = table[table[cfg.taxonomic_rank] == "TOTAL"].iloc[0]
    n_ranks_hit = int(((table[cfg.taxonomic_rank] != "TOTAL") & (table["hits"] > 0)).sum())
    print(f"[taxonomy] {cfg.taxonomic_rank}s hit: {n_ranks_hit} | overall "
          f"{int(total['hits'])}/{int(total['DB_n'])} ({total['hits_rate_pct']}%)",
          file=log, flush=True)
    return {
        "status": "ok",
        "tables": {f"{cfg.taxonomic_rank}_hitrate": table},
        "summary": {
            "rank": cfg.taxonomic_rank,
            f"{cfg.taxonomic_rank}s_with_hits": n_ranks_hit,
            "db_universe_n": int(total["DB_n"]),
            "overall_hitrate_pct": float(total["hits_rate_pct"]),
        },
    }
