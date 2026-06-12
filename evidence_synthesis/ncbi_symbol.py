"""
Resolve a RefSeq protein accession to its gene symbol via NCBI E-utilities
(protein -> gene -> official symbol), with an on-disk cache. Best-effort: on any
network/parse failure it returns None and the caller falls back to the cleaned title.
Stdlib only (urllib), so no new dependency.
"""
from __future__ import annotations
import json
import os
import re
import time
import urllib.parse
import urllib.request

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_CACHE: dict | None = None
_CACHE_PATH: str | None = None


def clean_title(subject_title: str) -> str:
    """Human-readable protein name from a subject title: drop the leading accession,
    the trailing [Organism], 'isoform X', and RefSeq's [[ ]] subscript markup."""
    t = re.sub(r"^\S+\s+", "", str(subject_title))      # leading accession
    t = re.sub(r"\s*\[[^\]]*\]\s*$", "", t)             # trailing [Organism]
    t = re.sub(r",?\s*isoform\s+\S+", "", t)            # isoform X
    return t.replace("[[", "").replace("]]", "").strip()


def _load_cache(path):
    global _CACHE, _CACHE_PATH
    _CACHE_PATH = path
    if _CACHE is None:
        try:
            _CACHE = json.load(open(path)) if os.path.exists(path) else {}
        except Exception:
            _CACHE = {}
    return _CACHE


def _save_cache():
    if _CACHE_PATH and _CACHE is not None:
        try:
            json.dump(_CACHE, open(_CACHE_PATH, "w"))
        except Exception:
            pass


def _get(url):
    with urllib.request.urlopen(url, timeout=25) as r:
        return json.load(r)


def gene_symbol(accession: str, email: str, cache_path: str,
                tool: str = "genomiss_evidence") -> str | None:
    """Official gene symbol for a protein accession, or None on failure. Cached."""
    cache = _load_cache(cache_path)
    if accession in cache:
        return cache[accession]
    sym = None
    try:
        common = f"&email={urllib.parse.quote(email)}&tool={tool}&retmode=json"
        el = _get(f"{_EUTILS}/elink.fcgi?dbfrom=protein&db=gene&id={accession}{common}")
        dbs = el["linksets"][0].get("linksetdbs", [])
        gid = next((d["links"][0] for d in dbs
                    if d.get("linkname") == "protein_gene" and d.get("links")), None)
        time.sleep(0.34)   # NCBI: <=3 requests/sec without an API key
        if gid:
            es = _get(f"{_EUTILS}/esummary.fcgi?db=gene&id={gid}{common}")
            doc = es["result"][str(gid)]
            sym = doc.get("nomenclaturesymbol") or doc.get("name") or None
            time.sleep(0.34)
    except Exception:
        sym = None
    if sym:                      # cache successes only; let transient failures retry
        cache[accession] = sym
        _save_cache()
    return sym
