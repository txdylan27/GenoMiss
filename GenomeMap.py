from typing import Dict, List, Optional

"""
Object Class GeneNode represents a single Gene, all of its protein isoforms, and any directly neighboring genes downstream of the 
current strand.
"""
class GeneNode:
    def __init__(
            self,
            gene_id: str,
            start_coord: int,
            end_coord: int,
            protein_isoforms: Dict[str, str] = None, # [ID -> sequence)
            neighbors: List['GeneNode'] = None,
            description: str = None
    ):
        self.gene_id: str = gene_id
        self.protein_isoforms: Dict[str, str] = protein_isoforms if protein_isoforms is not None else {}
        self.neighbors: List['GeneNode'] = neighbors if neighbors is not None else []
        self.start_coord: int = start_coord
        self.end_coord: int = end_coord
        self.description: str = description

    def add_neighbor(self, neighbor: 'GeneNode') -> None:
        """Add a neighboring GeneNode."""
        if neighbor not in self.neighbors:
            self.neighbors.append(neighbor)

    def add_protein_isoform(self, key: str, value: list) -> None:
        """Add or update a protein isoform."""
        self.protein_isoforms[key] = value[0]
        self.description = value[1]

    def get_longest_isoform(self) -> Optional[tuple]:
        """Return (isoform_id, sequence) of the longest isoform for a gene, or None if none exist."""
        if not self.protein_isoforms:
            return None
        return max(self.protein_isoforms.items(), key=lambda x: len(x[1]))
    
    def __repr__(self) -> str:
        return f"GeneNode(gene_id='{self.gene_id}', protein_isoforms={self.protein_isoforms}, neighbors={len(self.neighbors)})"

"""
Class GenomeMap keeps track of the head GeneNode of the positive and negative strand of each chromosome for an organism.
"""
class GenomeMap:
    def __init__(self, organism_name: str):
        self.organism_name: str = organism_name
        self.chromosomes: Dict[str, Dict[str, Optional[GeneNode]]] = {}

    def add_chromosome(self, chromosome_id: str) -> None:
        """Add a chromosome with positive and negative strand head nodes."""
        if chromosome_id not in self.chromosomes:
            self.chromosomes[chromosome_id] = {
                '+': None,
                '-': None
            }

    def set_head_node(self, chromosome_id: str, strand: str, head_node: GeneNode) -> None:
        """Set the head node for a specific chromosome strand.

        Args:
            chromosome_id: The chromosome identifier (e.g., 'chr1', '1', 'X')
            strand: Either '+' or '-'
            head_node: The GeneNode to set as the head
        """
        if strand not in ['+', '-']:
            raise ValueError("Strand must be '+' or '-'")

        if chromosome_id not in self.chromosomes:
            self.add_chromosome(chromosome_id)

        self.chromosomes[chromosome_id][strand] = head_node

    def get_head_node(self, chromosome_id: str, strand: str) -> Optional[GeneNode]:
        """Get the head node for a specific chromosome strand."""
        if strand not in ['+', '-']:
            raise ValueError("Strand must be '+' or '-'")

        if chromosome_id in self.chromosomes:
            return self.chromosomes[chromosome_id][strand]
        return None

    def list_chromosomes(self) -> list:
        """Return a list of all chromosome IDs in the genome."""
        return list(self.chromosomes.keys())

    def __repr__(self) -> str:
        return f"Genome(organism='{self.organism_name}', chromosomes={len(self.chromosomes)})"