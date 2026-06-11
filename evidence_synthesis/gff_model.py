"""
GFF3 parsing -> gene models (gene span/strand + representative-isoform exons).

Same feature handling as GenoMiss's own GFF walk (gene/mRNA/exon, tab-split,
attribute ID/Parent), so gene-ID linkage matches the run. Used for: (a) resolving
the genomic coords of a hit's two LOC genes from the run GFF, and (b) building the
dual-annotation locus figure by pulling all features overlapping a genomic window
from each GFF independently (both GFFs share the assembly -> coordinate overlay,
no LOC<->egapxtmp id map needed).
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional


def parse_attributes(attr: str) -> Dict[str, str]:
    out = {}
    for kv in attr.rstrip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out


@dataclass
class GeneModel:
    gene_id: str           # display id (Name attr, else ID attr)
    chrom: str
    start: int             # 1-based inclusive (GFF)
    end: int
    strand: str
    exons: List[tuple] = field(default_factory=list)   # representative isoform, genomic-sorted (start,end)
    source: str = ""       # annotation label, e.g. "NCBI" / "sheina2025"

    def exon_labels(self) -> List[str]:
        """Exon 1..n in transcription direction (reverse genomic order on '-' strand)."""
        n = len(self.exons)
        order = range(n) if self.strand != "-" else range(n - 1, -1, -1)
        labels = [""] * n
        for num, idx in enumerate(order, start=1):
            labels[idx] = f"Exon {num}"
        return labels


def find_genes(gff_path: str, names: set) -> Dict[str, dict]:
    """Single streaming pass: return {queried_name: gene_record} for each name found.
    Matches a `gene` feature whose Name or ID attribute contains the queried name."""
    want = set(names)
    found: Dict[str, dict] = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = parse_attributes(f[8])
            label = attrs.get("Name") or attrs.get("ID", "")
            for q in list(want):
                if q == label or q in label or q in f[8]:
                    found[q] = dict(chrom=f[0], start=int(f[3]), end=int(f[4]),
                                    strand=f[6], name=label, attrs=attrs)
                    want.discard(q)
            if not want:
                break
    return found


def models_in_window(gff_path: str, chrom: str, wstart: int, wend: int,
                     source: str = "") -> List[GeneModel]:
    """All genes overlapping [wstart, wend] on `chrom`, each with its representative
    isoform's exons (the mRNA with the most exons; ties -> longest span)."""
    genes: Dict[str, dict] = {}          # gene ID -> record
    mrna_parent: Dict[str, str] = {}     # mRNA ID -> gene ID
    exons: Dict[str, List[tuple]] = {}   # mRNA ID -> [(start,end)]
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != chrom:
                continue
            s, e = int(f[3]), int(f[4])
            if e < wstart or s > wend:
                continue
            ftype, attrs = f[2], parse_attributes(f[8])
            if ftype == "gene":
                gid = attrs.get("ID", "")
                genes[gid] = dict(chrom=f[0], start=s, end=e, strand=f[6],
                                  name=attrs.get("Name") or gid)
            elif ftype in ("mRNA", "transcript"):
                mid = attrs.get("ID", "")
                mrna_parent[mid] = attrs.get("Parent", "")
            elif ftype == "exon":
                for parent in attrs.get("Parent", "").split(","):
                    exons.setdefault(parent, []).append((s, e))

    # representative isoform per gene
    by_gene_mrna: Dict[str, List[str]] = {}
    for mid, gid in mrna_parent.items():
        by_gene_mrna.setdefault(gid, []).append(mid)

    models: List[GeneModel] = []
    for gid, g in genes.items():
        mrnas = by_gene_mrna.get(gid, [])
        best, best_exs = None, []
        for mid in mrnas:
            exs = sorted(exons.get(mid, []))
            if not exs:
                continue
            span = exs[-1][1] - exs[0][0]
            if best is None or (len(exs), span) > (len(best_exs), best[1] - best[0]):
                best, best_exs = (exs[0][0], exs[-1][1]), exs
        if not best_exs:  # gene with no parsed exons (e.g. pseudogene); draw the gene box only
            best_exs = [(g["start"], g["end"])]
        models.append(GeneModel(gene_id=g["name"], chrom=g["chrom"], start=g["start"],
                                end=g["end"], strand=g["strand"], exons=best_exs, source=source))
    models.sort(key=lambda m: m.start)
    return models
