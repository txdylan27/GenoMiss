"""
Shared taxonomy helpers via ete3 NCBITaxa — consolidates the get_order()/order()
function that was duplicated three times (classify_orders.py, run_uncapped.py,
analyze_pair.py). Rank is configurable (default 'order') for later generality.
"""
from __future__ import annotations
import collections
import functools
from typing import Iterable

import pandas as pd
from ete3 import NCBITaxa

_NCBI = None


def _ncbi() -> NCBITaxa:
    global _NCBI
    if _NCBI is None:
        _NCBI = NCBITaxa()
    return _NCBI


@functools.lru_cache(maxsize=None)
def assign_rank(name: str, rank: str = "order") -> str:
    """Map an organism name to its taxon at `rank`. Falls back full name ->
    genus+species -> genus before giving up (handles strain/partial names)."""
    ncbi = _ncbi()
    for cand in (name, " ".join(name.split()[:2]), name.split()[0]):
        d = ncbi.get_name_translator([cand])
        if d:
            taxid = d[cand][0]
            lineage = ncbi.get_lineage(taxid)
            ranks = ncbi.get_rank(lineage)
            names = ncbi.get_taxid_translator(lineage)
            for tid in lineage:
                if ranks[tid] == rank:
                    return names[tid]
            return f"No{rank.capitalize()}Rank"
    return "Unmatched"


def rank_hitrate_table(db_organisms: Iterable[str], hit_organisms: Iterable[str],
                       rank: str = "order", hit_label: str = "hits") -> pd.DataFrame:
    """Per-rank table: DB_n (denominator universe), <hit_label>, and hit rate %.
    Sorted by DB_n descending; a TOTAL row appended."""
    db = sorted(set(db_organisms))
    hit = set(hit_organisms)
    rank_of = {o: assign_rank(o, rank) for o in set(db) | hit}
    db_counts = collections.Counter(rank_of[o] for o in db)
    hit_counts = collections.Counter(rank_of[o] for o in hit if o in rank_of)
    rows = []
    for r, dbn in db_counts.most_common():
        h = hit_counts.get(r, 0)
        rows.append((r, dbn, h, round(100.0 * h / dbn, 1) if dbn else 0.0))
    df = pd.DataFrame(rows, columns=[rank, "DB_n", hit_label, f"{hit_label}_rate_pct"])
    total = pd.DataFrame([("TOTAL", len(db), len(hit & set(rank_of)),
                           round(100.0 * len(hit) / len(db), 1) if db else 0.0)],
                         columns=df.columns)
    return pd.concat([df, total], ignore_index=True)
