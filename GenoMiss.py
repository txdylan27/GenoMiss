import argparse
import bisect
import os
import pandas as pd
import shutil
import subprocess
import sys
from tqdm import tqdm
import regex as re
import GenomeMap
import scoring
import output_formatter

# FLAGS AND GLOBAL VARIABLES
ORGANISM_NAME_DETECTED = False
DEBUG_CONTROL = False
tested_genes = set()  # Global set to track fragmented genes being tested (per strand)
fragmented_df = None  # Global to store fragmented genes report
best_diamond_hits = {}  # Global dictionary to store best DIAMOND hit per gene

# Tracking structure for control gene losses
control_losses = {
    'missing_from_map': set(),
    'broken_neighbors': set(),
    'missing_proteins': set(),
    'no_diamond_hits': set(),
    'failed_alignment_overlap': set(),
    'failed_query_coverage': set(),
    'failed_percent_identity': set(),
    'failed_uncharacterized': set(),
    'passed_all_filters': set(),
    'not_tested': set()
}
    
"""
Gene types to include in the genome map analysis
Available gene_biotype values: protein_coding, lncRNA, miRNA, tRNA, rRNA, snRNA, snoRNA, misc_RNA, guide_RNA
Algorithm currently only works with genes that produce actual proteins and are included in the .faa, so we restrict
the possible genes as such. In the future, if we want to concatenate a theorized protein coding gene from a gene currently
annotated as a lncRNA, we would need ot add those included gene types in this set, and also recycle dylan's old methods
for constructing proteins from CDS regions to manually construct the theoretical protein.
"""
INCLUDED_GENE_TYPES = {'protein_coding'}

# Defining all functions for the program
def get_default_num_threads():
    """Function used to identify the number of threads on the user's system and then use half of them for DIAMOND if no threads are inputted."""
    
    total_threads = os.cpu_count() # Attempting to get the thread count for the user's system.
    if total_threads is None: 
        total_threads = 4 # Included in case - for some reason - Python is not able to identify their thread count.
        
    return max(1, total_threads // 2) # Returning the max of either 1 thread or half of the total_threads. The 1 is a fall-back in the rare instance of a system having 0 or only 1 thread.

def construct_faa_string(proteinID, proteinDescription, proteinSequence):
    # Handle empty/None description
    desc = f"{proteinDescription}" if proteinDescription else ""
    return f">{proteinID} {desc}\n{proteinSequence}\n"

def chromosome_processor_unfused_allisoforms(chrom, strand, headNode: GenomeMap.GeneNode, output_prefix, fused_chrom_hits_df, longest_only):
    """
    Traverse all genes in the chromosome strand graph using BFS (Breadth-First Search).
    Handles branching paths from overlapping genes and convergence points.
    """
    tempFAAstring = ""
    visited = set()  # Track visited GeneNode objects to avoid duplicates when paths converge
    queue = [headNode]  # Initialize queue with head node (BFS uses FIFO)

    # Creating sets of product_1 and product_2 to only write out fused gene parts to the faa
    product1_set = set(fused_chrom_hits_df['product_1'])
    product2_set = set(fused_chrom_hits_df['product_2'])

    while queue:
        currentNode = queue.pop(0)  # Dequeue: pop from front (BFS: process level-by-level)

        # Skip if we've already processed this node (handles convergence and duplicate gene names)
        if currentNode in visited:
            continue
        visited.add(currentNode)

        # Handling of longest isoform mode
        if longest_only:
            current_node_isoforms = currentNode.get_longest_isoform()
        else:
            current_node_isoforms = currentNode.protein_isoforms.items()

        # Process all protein isoforms for this gene
        for isoformID, isoformAA in current_node_isoforms:
            # Checking for matches
            if isoformID in product1_set or isoformID in product2_set:
                protString = construct_faa_string(isoformID, currentNode.description, isoformAA)
                tempFAAstring += protString

        # Enqueue all neighbors (handles branching from overlapping genes)
        if currentNode.neighbors:
            queue.extend(currentNode.neighbors)

    temp_fasta_path = f"{output_prefix}/temp_fasta.faa"
    with open(temp_fasta_path, "w") as temp_fasta:
        temp_fasta.write(tempFAAstring)

    return temp_fasta_path

def chromosome_processor_fused_allisoforms(chrom, strand, headNode: GenomeMap.GeneNode, output_prefix, longest_only):
    """
    Traverse all genes and create fused proteins between each gene and its neighbors.
    Uses BFS to handle branching/converging paths from overlapping genes.
    Creates all possible fusions: each isoform from gene1 × each isoform from gene2.
    Returns tuple: (temp_fasta_path, metadata_dataframe)
    """
    tempFAAstring = ""
    fused_metadata = []  # Store metadata for each fusion
    visited_nodes = set()  # Track visited GeneNode objects for graph traversal
    visited_edges = set()  # Track processed edges to avoid duplicate fusions
    queue = [headNode]

    while queue:
        currentNode = queue.pop(0)  # BFS: pop from front

        # Skip if we've already visited this node (handles convergence and duplicate gene names)
        if currentNode in visited_nodes:
            continue
        visited_nodes.add(currentNode)

        # Create fusions with all neighbors
        if currentNode.neighbors:
            for neighbor in currentNode.neighbors:
                # Create unique edge identifier (directed edge)
                edge = (currentNode.gene_id, neighbor.gene_id)

                # Skip if we've already processed this edge
                if edge in visited_edges:
                    continue
                visited_edges.add(edge)

                if DEBUG_CONTROL:
                    # Check if this is a part1 + part2 fusion from GenoFrag
                    if fragmented_df is not None:
                        # Check if currentNode is part1 and neighbor is part2
                        if currentNode.gene_id.endswith('_part1') and neighbor.gene_id.endswith('_part2'):
                            # Extract base gene name (remove _part1 suffix)
                            gene_base = currentNode.gene_id[:-6]  # Remove "_part1"
                            # Verify neighbor matches
                            if neighbor.gene_id == f"{gene_base}_part2":
                                tested_genes.add(gene_base)
                                # Check if both nodes have protein isoforms
                                part1_protein_count = len(currentNode.protein_isoforms)
                                part2_protein_count = len(neighbor.protein_isoforms)
                                if part1_protein_count == 0 or part2_protein_count == 0:
                                    control_losses['missing_proteins'].add(gene_base)
                # Handling of longest isoform mode
                if longest_only:
                    current_node_isoforms = currentNode.get_longest_isoform()
                    neighor_node_isoforms = neighbor.get_longest_isoform()
                else:
                    current_node_isoforms = currentNode.protein_isoforms.items()
                    neighor_node_isoforms = neighbor.protein_isoforms.items()

                # Create all possible isoform combinations (Cartesian product)
                for isoformID1, isoformAA1 in current_node_isoforms:
                    for isoformID2, isoformAA2 in neighor_node_isoforms:
                        # Fuse IDs with underscore
                        fusedID = f"{isoformID1}_{isoformID2}"
                        # Concatenate protein sequences
                        fusedSequence = isoformAA1 + isoformAA2
                        # Concatenate descriptions (or use None if either is missing)
                        fusedDescription = None
                        if currentNode.description and neighbor.description:
                            fusedDescription = f"{currentNode.description} + {neighbor.description}"

                        # Store metadata
                        fused_metadata.append({
                            "fused_protein": fusedID,
                            "fused_product": fusedDescription if fusedDescription else f"{isoformID1} + {isoformID2}",
                            "fused_gene_len": len(fusedSequence),
                            "gene_1": currentNode.gene_id,
                            "product_1": isoformID1,
                            "gene_1_len": len(isoformAA1),
                            "gene_2": neighbor.gene_id,
                            "product_2": isoformID2,
                            "gene_2_len": len(isoformAA2)
                        })

                        # Construct FAA entry
                        protString = construct_faa_string(fusedID, fusedDescription, fusedSequence)
                        tempFAAstring += protString

            # Add all neighbors to queue for continued traversal
            queue.extend(currentNode.neighbors)

    temp_fasta_path = f"{output_prefix}/temp_fasta.faa"
    with open(temp_fasta_path, "w") as temp_fasta:
        temp_fasta.write(tempFAAstring)

    metadata_df = pd.DataFrame(fused_metadata)
    return temp_fasta_path, metadata_df


def unfused_diamond_alignment(chrom, strand_name, diamond_path, diamond_db, genome_faa_path, num_threads, output_prefix, diamond_sensitivity=None, taxon_exclude=None, max_target_seqs=100):
    """Function used to get the original alignment scores to use as a comparison with the fused gene alignments."""

    diamond_output = f"{output_prefix}/{chrom}_{strand_name}_control_diamond_results.tsv"
    diamond_command = [diamond_path, "blastp", "--db", diamond_db, "--query", genome_faa_path, "--out", diamond_output,
                       "--outfmt", "6", "qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "sstart", "send",
                       "pident", "nident", "mismatch", "evalue", "bitscore", "length", "qcovhsp", "scovhsp", "qtitle", "stitle",
                       "--header", "--evalue", "1e-5", "--threads", num_threads,
                       "--max-target-seqs", str(max_target_seqs)]

    # Add sensitivity parameter if specified
    if diamond_sensitivity:
        diamond_command.append(diamond_sensitivity)

    # Add taxon exclusion if specified
    if taxon_exclude:
        diamond_command.extend(["--taxon-exclude", taxon_exclude])

    try:
        result = subprocess.run(diamond_command, check=True, text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        # Check if error is due to missing taxonomy info
        if "taxonomy information" in e.stderr and taxon_exclude:
            print(f"WARNING: Database lacks taxonomy info. Retrying without --taxon-exclude...")
            # Remove taxon-exclude and retry
            diamond_command = [x for x in diamond_command if x != "--taxon-exclude" and x != taxon_exclude]
            result = subprocess.run(diamond_command, check=True, text=True, capture_output=True)
        else:
            print(f"Error during DIAMOND execution for {chrom}, {strand_name} strand.")
            print(f"DIAMOND stderr: {e.stderr}")
            print(f"DIAMOND stdout: {e.stdout}")
            raise

    headers = ["protein", "query_length", "subject_id", "subject_length", "start_of_alignment_in_query", "end_of_alignment_in_query",
               "start_of_alignment_in_subject", "end_of_alignment_in_subject", "percentage_of_identical_matches", "number_of_identical_matches",
               "number_of_mismatches", "expected_value", "bit_score", "alignment_length", "query_coverage", "subject_coverage",
               "query_title", "subject_title"]

    diamond_df = pd.read_csv(f"{output_prefix}/{chrom}_{strand_name}_control_diamond_results.tsv", sep="\t", skiprows = 3, names = headers)
    return diamond_df
    
    
def fused_diamond_alignment(chrom, strand_name, diamond_path, db, temp_fasta_path, fused_metadata_df, num_threads, ident_cutoff, output_prefix, diamond_sensitivity=None, taxon_exclude=None, max_target_seqs=100):
    """Function used to fuse neigboring genes and run the DIAMOND protein alignment script on the fused genes. Returns a dataframe containing unique genes that
    fit the user-inputted filtering criteria."""

    # Creating the diamond output file path and the diamond command that will be ran in the subprocess function.
    diamond_output = f"{output_prefix}/{chrom}_{strand_name}_diamond_results.tsv"
    diamond_command = [diamond_path, "blastp", "--db", db, "--query", temp_fasta_path, "--out", diamond_output,
                       "--outfmt", "6", "qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "sstart", "send",
                       "pident", "nident", "mismatch", "evalue", "bitscore", "length", "qcovhsp", "scovhsp", "qtitle", "stitle",
                       "--header", "--evalue", "1e-5", "--threads", num_threads,
                       "--max-target-seqs", str(max_target_seqs)]

    # Add sensitivity parameter if specified
    if diamond_sensitivity:
        diamond_command.append(diamond_sensitivity)

    # Add taxon exclusion if specified
    if taxon_exclude:
        diamond_command.extend(["--taxon-exclude", taxon_exclude])

    # Running the diamond command via subprocess.
    try:
        result = subprocess.run(diamond_command, check=True, text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        # Check if error is due to missing taxonomy info
        if "taxonomy information" in e.stderr and taxon_exclude:
            print(f"WARNING: Database lacks taxonomy info. Retrying without --taxon-exclude...")
            # Remove taxon-exclude and retry
            diamond_command = [x for x in diamond_command if x != "--taxon-exclude" and x != taxon_exclude]
            result = subprocess.run(diamond_command, check=True, text=True, capture_output=True)
        else:
            print(f"Error during DIAMOND execution for {chrom}, {strand_name} strand.")
            print(f"DIAMOND stderr: {e.stderr}")
            print(f"DIAMOND stdout: {e.stdout}")
            raise

    # As the diamond headers are very obfuscated, included proper header titles.
    headers = ["fused_protein", "query_length", "subject_id", "subject_length", "start_of_alignment_in_query", "end_of_alignment_in_query", 
               "start_of_alignment_in_subject", "end_of_alignment_in_subject", "percentage_of_identical_matches", "number_of_identical_matches", 
               "number_of_mismatches", "expected_value", "bit_score", "alignment_length", "query_coverage", "subject_coverage", 
               "query_title", "subject_title"]
    
    # Converting the .tsv provided by diamond to a dataframe while skipping useless rows (the first 3) and renaming the headers.
    diamond_df = pd.read_csv(f"{output_prefix}/{chrom}_{strand_name}_diamond_results.tsv", sep="\t", skiprows = 3, names = headers)

    if DEBUG_CONTROL and tested_genes:
        # Get expected fusions
        expected_fused_ids = set()
        for _, row in fused_metadata_df.iterrows():
            if pd.notna(row.get('gene_1')) and pd.notna(row.get('gene_2')):
                if row['gene_1'].endswith('_part1') and row['gene_2'].endswith('_part2'):
                    gene_base = row['gene_1'][:-6]
                    if row['gene_2'] == f"{gene_base}_part2" and gene_base in tested_genes:
                        expected_fused_ids.add((gene_base, row['fused_protein']))

        # Check which are missing from DIAMOND
        diamond_fused_ids = set(diamond_df['fused_protein'])
        for gene_base, fused_id in expected_fused_ids:
            if fused_id not in diamond_fused_ids:
                control_losses['no_diamond_hits'].add(gene_base)

    # Merging the diamond results dataframe with the fused protein dataframe. Left join keeps all fusions.
    fused_diamond_df = fused_metadata_df.merge(diamond_df, how='left', on='fused_protein', sort=False)

    # Track best DIAMOND hit for each gene for diagnostic reporting
    if DEBUG_CONTROL and tested_genes:
        # Get all rows for tested genes
        mask = (
            fused_diamond_df['gene_1'].str.endswith('_part1', na=False) &
            fused_diamond_df['gene_2'].str.endswith('_part2', na=False)
        )
        tested_rows = fused_diamond_df[mask].copy()
        tested_rows['gene_base'] = tested_rows['gene_1'].str[:-6]
        tested_rows = tested_rows[tested_rows['gene_base'].isin(tested_genes)]

        # For each tested gene, find its best hit (highest bitscore) and calculate diagnostics
        if not tested_rows.empty:
            for gene in tested_genes:
                gene_rows = tested_rows[tested_rows['gene_base'] == gene]
                if not gene_rows.empty:
                    # Filter out NaN bitscores
                    gene_rows_valid = gene_rows[gene_rows['bit_score'].notna()]
                    if not gene_rows_valid.empty:
                        # Get row with highest bitscore
                        best_hit_idx = gene_rows_valid['bit_score'].idxmax()
                        best_hit_row = gene_rows_valid.loc[best_hit_idx]

                        # Calculate detailed alignment overlap diagnostics
                        gene_1_len = best_hit_row['gene_1_len']
                        align_start = best_hit_row['start_of_alignment_in_query']
                        align_end = best_hit_row['end_of_alignment_in_query']

                        # Check overlap on both sides of the boundary
                        part1_overlap = gene_1_len - align_start  # How far into part1 from the end
                        part2_overlap = align_end - gene_1_len     # How far into part2 from the start

                        # Determine overlap status
                        if align_start < (gene_1_len - 10) and align_end > (gene_1_len + 10):
                            overlap_status = "PASS: >=10 AA on both sides"
                        elif align_start >= gene_1_len:
                            overlap_status = f"FAIL: No part1 coverage (alignment starts at {align_start}, part1 ends at {gene_1_len})"
                        elif align_end <= gene_1_len:
                            overlap_status = f"FAIL: No part2 coverage (alignment ends at {align_end}, part2 starts at {gene_1_len})"
                        else:
                            # Has some overlap but not enough
                            part1_depth = max(0, gene_1_len - align_start)
                            part2_depth = max(0, align_end - gene_1_len)
                            if part1_depth < 10:
                                overlap_status = f"FAIL: Part1 overlap only {part1_depth} AA (need 10)"
                            elif part2_depth < 10:
                                overlap_status = f"FAIL: Part2 overlap only {part2_depth} AA (need 10)"
                            else:
                                overlap_status = f"PASS: Part1={part1_depth} AA, Part2={part2_depth} AA"

                        # Store comprehensive hit information
                        best_diamond_hits[gene] = {
                            'fused_protein': best_hit_row['fused_protein'],
                            'subject_id': best_hit_row['subject_id'],
                            'subject_title': best_hit_row.get('subject_title', ''),
                            'bitscore': best_hit_row['bit_score'],
                            'qcoverage': best_hit_row['query_coverage'],
                            'pident': best_hit_row['percentage_of_identical_matches'],
                            'evalue': best_hit_row['expected_value'],
                            'alignment_start': align_start,
                            'alignment_end': align_end,
                            'gene_1_len': gene_1_len,
                            'gene_2_len': best_hit_row['gene_2_len'],
                            'fused_gene_len': best_hit_row['fused_gene_len'],
                            'overlap_status': overlap_status,
                            'part1_overlap_depth': max(0, gene_1_len - align_start),
                            'part2_overlap_depth': max(0, align_end - gene_1_len)
                        }

    # Only filtering for fused genes that have overlaps with at least 10 AAs, maybe make this user-inputted
    fused_hits_df = fused_diamond_df[
        (fused_diamond_df["start_of_alignment_in_query"] < (fused_diamond_df["gene_1_len"]-10)) &
        (fused_diamond_df["end_of_alignment_in_query"] > (fused_diamond_df["gene_1_len"]+10)) &
        (fused_diamond_df["query_coverage"] >= 50) &
        (fused_diamond_df["percentage_of_identical_matches"] >= 50) &
        (~fused_diamond_df["subject_title"].str.contains("uncharacterized", case=False, na=False))
    ].copy()

    # Track which genes were COMPLETELY filtered out and WHY
    if DEBUG_CONTROL and tested_genes:
        # Step 1: Identify genes that passed (have at least one hit in fused_hits_df)
        mask = (
            fused_hits_df['gene_1'].str.endswith('_part1', na=False) &
            fused_hits_df['gene_2'].str.endswith('_part2', na=False)
        )
        passed_rows = fused_hits_df[mask]

        passed_genes = set()
        if not passed_rows.empty:
            gene_bases = passed_rows['gene_1'].str[:-6]
            passed_genes = set(gene_bases[gene_bases.isin(tested_genes)])

        control_losses['passed_all_filters'].update(passed_genes)

        # Step 2: Identify genes that were completely filtered out
        filtered_out_genes = tested_genes - passed_genes

        # Step 3: For completely filtered-out genes, determine which filters blocked them
        if filtered_out_genes:
            # Get all rows for tested genes from fused_diamond_df
            mask = (
                fused_diamond_df['gene_1'].str.endswith('_part1', na=False) &
                fused_diamond_df['gene_2'].str.endswith('_part2', na=False)
            )
            tested_rows = fused_diamond_df[mask].copy()
            tested_rows['gene_base'] = tested_rows['gene_1'].str[:-6]
            tested_rows = tested_rows[tested_rows['gene_base'].isin(filtered_out_genes)]

            if not tested_rows.empty:
                # For each filtered-out gene, check which filters blocked ALL its hits
                for gene in filtered_out_genes:
                    gene_rows = tested_rows[tested_rows['gene_base'] == gene]

                    if not gene_rows.empty:
                        # Check if ALL hits failed each filter criterion
                        # A gene is counted as failed_X if ALL its hits failed filter X

                        # Filter 1: Alignment overlap
                        passes_overlap = (
                            (gene_rows["start_of_alignment_in_query"] < (gene_rows["gene_1_len"]-10)) &
                            (gene_rows["end_of_alignment_in_query"] > (gene_rows["gene_1_len"]+10))
                        )
                        if not passes_overlap.any():  # ALL hits failed this filter
                            control_losses['failed_alignment_overlap'].add(gene)

                        # Filter 2: Query coverage
                        passes_coverage = (gene_rows["query_coverage"] >= 50)
                        if not passes_coverage.any():  # ALL hits failed this filter
                            control_losses['failed_query_coverage'].add(gene)

                        # Filter 3: Percent identity
                        passes_identity = (gene_rows["percentage_of_identical_matches"] >= 50)
                        if not passes_identity.any():  # ALL hits failed this filter
                            control_losses['failed_percent_identity'].add(gene)

                        # Filter 4: Uncharacterized
                        passes_uncharacterized = ~gene_rows["subject_title"].str.contains("uncharacterized", case=False, na=False)
                        if not passes_uncharacterized.any():  # ALL hits failed this filter
                            control_losses['failed_uncharacterized'].add(gene)
    
    
    os.remove(temp_fasta_path)
    
    return fused_hits_df, fused_diamond_df

def detect_organism_name(gffFile):
    global ORGANISM_NAME_DETECTED # necessary to tell later assignment to change outside global variable
    with open(gffFile) as gff:
        headers_only = [line for line in gff if line.startswith('#')]
        for header in headers_only:
            regex_result = re.search(r'NCBI\s+(.+?)\s+Annotation', header)
            if regex_result:
                organism_name = regex_result.group(1)
                print("Detected Organism Name: " + organism_name)
                ORGANISM_NAME_DETECTED = True
                return organism_name
    print("Organism name ")
    return None #If organism name is not found, ORGANISM_NAME_DETECTED will remain false

def detect_taxon_id(gffFile):
    """Detect NCBI taxonomy ID from GFF file ##species header line."""
    with open(gffFile) as gff:
        for line in gff:
            if line.startswith('##species'):
                # Pattern matches URLs like: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7227
                regex_result = re.search(r'id=(\d+)', line)
                if regex_result:
                    taxon_id = regex_result.group(1)
                    print(f"Detected Taxon ID: {taxon_id}")
                    return taxon_id
            # Stop reading after headers
            if not line.startswith('#'):
                break
    print("WARNING: Taxon ID not found in GFF file. Self-hits will not be excluded.")
    return None

def get_chromosome_names(gffFile):
    # Alternative apporach would be to use pandas to create a dataframe and then use df['seqid'].unique(), but this is
    # faster computationally and more memory efficient.
    chromosomes = set()
    with open(gffFile) as gff:
        for line in gff:
            if not line.startswith('#'):
                seqid = line.split('\t')[0]
                chromosomes.add(seqid)
    return list(chromosomes)

def count_genes_per_chromosome(gffFile):
    gene_counts = {}

    with open(gffFile) as gff:
        for line in gff:
            if not line.startswith('#'):
                fields = line.split('\t')
                seqid = fields[0]
                feature_type = fields[2]

                # Count only protein coding gene features
                if feature_type == 'gene':
                    # Extract gene_biotype and check if it's in the included types
                    attributes = fields[8]
                    attr_dict = {}
                    for attr in attributes.strip().split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            attr_dict[key] = value
                    gene_biotype = attr_dict.get('gene_biotype', '')

                    if gene_biotype in INCLUDED_GENE_TYPES:
                        gene_counts[seqid] = gene_counts.get(seqid, 0) + 1

    return gene_counts

# Returns a set with all of the chromosomes we are interested in
def genecount_filtered_chromosomes(gffFile, minGenePerX):
    x_w_counts = count_genes_per_chromosome(gffFile)
    return {x for x in x_w_counts if x_w_counts[x] > minGenePerX}

# Returns a dictionary of protein ID's to [sequence, description]
def process_faa(faaFile):
    faaDict = {}
    # Made description optional to handle headers like ">PROTEIN_ID" without description
    pattern = r'>(\S+)(?:\s+(.+?)(?:\s+\[([^\]]+)\])?)?$'

    protein_id = None
    description = None
    curProt = ""

    with open(faaFile) as faa:
        for line in faa:
            if line.startswith('>'):
                # Save previous protein before starting new one
                if protein_id:
                    faaDict[protein_id] = [curProt, description]

                # Reset for new protein
                curProt = ""
                match = re.search(pattern, line)
                if match:
                    protein_id = match.group(1)
                    description = match.group(2)  # Can be None if no description
                    #organism = match.group(3)  # Optional: store if needed
                else:
                    print(f"WARNING: Problem parsing line: {line.strip()} in the .faa file!")
                    protein_id = None
                    description = None
            else:
                # Accumulate sequence, stripping newline
                curProt += line.strip()

        # Save last protein when the loop is over
        if protein_id:
            faaDict[protein_id] = [curProt, description]

    return faaDict

def add_edge(upstreamNode: GenomeMap.GeneNode, downstreamNode: GenomeMap.GeneNode):
    upstreamNode.neighbors.append(downstreamNode)

def find_last_nonoverlapping_gene(nodeList, nuNode):
    for node in reversed(nodeList):
        if node.end_coord < nuNode.start_coord:
            return node

def process_nextgene(theMap:GenomeMap.GenomeMap, nodeList: list, open_gene_dict: dict, firstNode: bool, fields: list):
    global nuNode
    curStart = int(fields[3])
    curEnd = int(fields[4])
    curStrandedness = fields[6]
    attributes = fields[8]
    attr_dict = {}
    for attr in attributes.strip().split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attr_dict[key] = value
    gene_id = attr_dict.get('gene') or attr_dict.get('Name') or attr_dict.get('locus_tag')

    if gene_id:
        #name = attr_dict.get('Name', gene_id)
        nuNode = GenomeMap.GeneNode(gene_id, curStart, curEnd)
        if firstNode:
            curNode = nuNode
            theMap.set_head_node(fields[0], curStrandedness, curNode)
            firstNode = False
        else:# Check open genes
            genes_that_got_closed = []
            # Open genes are all genes that don't have a forward connection yet
            for open_gene in open_gene_dict[curStrandedness]:
                if nuNode.start_coord > open_gene.end_coord:
                    add_edge(open_gene, nuNode)
                    genes_that_got_closed.append(open_gene)
                else: # If we were unable to connect one of the open genes, add an additional connection to the last gene without an overlap
                    backconnect_node = find_last_nonoverlapping_gene(nodeList, nuNode)
                    if backconnect_node:
                        add_edge(backconnect_node, nuNode)
                    else:
                        # Edge case: First two genes on chromosome/strand overlap. No non-overlapping gene exists to backconnect to.
                        # Connect to overlapping gene anyway to ensure all non-head nodes have incoming edges.
                        add_edge(open_gene, nuNode)
                        genes_that_got_closed.append(open_gene)
            for closedGene in genes_that_got_closed:
                open_gene_dict[curStrandedness].remove(closedGene)

        # Add current node to the open_gene_dict
        open_gene_dict[curStrandedness].append(nuNode)
    else:
        print("ERROR: Failed to obtain gene_id for: " + "\t".join(fields))
        exit()

    return nuNode

def build_genomemap(organismName, proteomeFile, gffFile, minGenePerX):
    """
    The main method for constructing a strand-specific genome graph, showing which genes are considered "connection"
    candidates to other genes
    """
    chromosomes = genecount_filtered_chromosomes(gffFile, minGenePerX)
    print(f"Chromsomes/contigs that met the minimum gene count threshold of {minGenePerX}:\n{chr(10).join(sorted(chromosomes))}")
    proteins = process_faa(proteomeFile)
    theMap = GenomeMap.GenomeMap(organismName)
    with open(gffFile) as gff:

        curChrom = ""
        curNode = GenomeMap.GeneNode
        first_node_per_strand = {"+": True, "-": True}  # Track first node for each strand separately
        open_genes = {"+": [], "-": []}
        nodeList = {"+": [], "-": []}

        for line in gff:
            if line.startswith('#'):
                continue

            fields = line.split('\t')
            seqid = fields[0] # Chromosome or Scaffold ID
            if seqid not in chromosomes:
                continue
            if seqid != curChrom:
                curChrom = seqid
                first_node_per_strand = {"+": True, "-": True}  # Reset for new chromosome
                nodeList = {"+": [], "-": []}  # Reset node lists for new chromosome
                open_genes = {"+": [], "-": []} # Reset open_genes too, just in case.
                curNode = None

            feature_type = fields[2]
            if feature_type == 'exon':
                continue

            if feature_type == 'gene':
                # Extract gene_biotype and check if it's in the included types
                attributes = fields[8]
                attr_dict = {}
                for attr in attributes.strip().split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value
                gene_biotype = attr_dict.get('gene_biotype', '')

                if gene_biotype in INCLUDED_GENE_TYPES:
                    strand = fields[6]
                    curNode = process_nextgene(theMap, nodeList[strand], open_genes, first_node_per_strand[strand], fields)
                    bisect.insort(nodeList[strand], curNode, key=lambda x: x.end_coord)
                    first_node_per_strand[strand] = False  # Mark that we've seen the first node for this strand
                else:
                    # Skip genes not in INCLUDED_GENE_TYPES (tRNA, lncRNA, miRNA, etc.)
                    curNode = None

            if feature_type == 'CDS':
                # Only add CDS if curNode points to a valid gene (one that passed the gene_biotype filter)
                if curNode is not None:
                    attributes = fields[8]
                    attr_dict = {}
                    for attr in attributes.strip().split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            attr_dict[key] = value
                    protein_id = attr_dict.get('Name')
                    if protein_id and protein_id in proteins:
                        curNode.add_protein_isoform(protein_id, proteins[protein_id])

    return theMap

def compare_bitscores(df_fused, df_control):
    # Create lookup dictionary for quick access to control bitscores
    bitscore_dict = df_control.set_index("protein")["bit_score"].to_dict()

    # Define comparison function for each row
    def is_fused_higher(row):
        b1 = bitscore_dict.get(row["product_1"], None)
        b2 = bitscore_dict.get(row["product_2"], None)
        if b1 is None or b2 is None:
            return None  # skip if missing
        return row["bit_score"] > b1 and row["bit_score"] > b2

    # Apply comparison to each fused protein row
    df_fused["fused_higher"] = df_fused.apply(is_fused_higher, axis=1)

    # Split into two new DataFrames (format preserved)
    df_fused_higher = df_fused[df_fused["fused_higher"] == True].copy()
    df_fused_lower = df_fused[df_fused["fused_higher"] == False].copy()

    return df_fused_higher, df_fused_lower

def calculate_intron_lengths(df_fused, gff_file):
    """
    Calculate theorized intron length between fused gene pairs from GFF file.

    The intron length is the genomic distance between gene_1's end and gene_2's start.
    For genes on the negative strand, we take the absolute distance.

    Args:
        df_fused (pd.DataFrame): Fused hits dataframe with gene_1 and gene_2 columns
        gff_file (str): Path to GFF annotation file

    Returns:
        dict: Mapping of (gene_1, gene_2) tuples to intron length
    """
    # Build a dictionary of gene coordinates from GFF
    gene_coords = {}  # {gene_id: (chrom, start, end, strand)}

    with open(gff_file) as gff:
        for line in gff:
            if line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type == 'gene':
                # Extract gene_biotype
                attributes = fields[8]
                attr_dict = {}
                for attr in attributes.strip().split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value

                gene_biotype = attr_dict.get('gene_biotype', '')
                if gene_biotype in INCLUDED_GENE_TYPES:
                    gene_id = attr_dict.get('gene') or attr_dict.get('Name') or attr_dict.get('locus_tag')
                    if gene_id:
                        chrom = fields[0]
                        start = int(fields[3])
                        end = int(fields[4])
                        strand = fields[6]
                        gene_coords[gene_id] = (chrom, start, end, strand)

    # Calculate intron lengths for each gene pair
    intron_lengths = {}
    unique_pairs = df_fused[['gene_1', 'gene_2']].drop_duplicates()

    for _, row in unique_pairs.iterrows():
        gene_1 = row['gene_1']
        gene_2 = row['gene_2']

        if gene_1 in gene_coords and gene_2 in gene_coords:
            chrom1, start1, end1, strand1 = gene_coords[gene_1]
            chrom2, start2, end2, strand2 = gene_coords[gene_2]

            # Only calculate if on same chromosome and strand
            if chrom1 == chrom2 and strand1 == strand2:
                # Intron length is the gap between genes
                # For + strand: gap = start2 - end1 - 1
                # For - strand: gap = start1 - end2 - 1 (but genes are stored in opposite order)
                # Since our algorithm processes genes in genomic order, we can use:
                gap = start2 - end1 - 1

                # Gap can be negative if genes overlap (which is valid in the graph)
                # Set overlapping genes to have 0 intron length
                intron_length = max(0, gap)
                intron_lengths[(gene_1, gene_2)] = intron_length
            else:
                # Different chromosomes or strands - shouldn't happen but handle it
                intron_lengths[(gene_1, gene_2)] = None
        else:
            # Gene not found in GFF
            intron_lengths[(gene_1, gene_2)] = None

    return intron_lengths

def calculate_scores_for_hits(df_fused, df_control, gff_file=None):
    """
    Calculate composite scores for all fused gene hits.

    Args:
        df_fused (pd.DataFrame): Fused hits dataframe
        df_control (pd.DataFrame): Control hits dataframe
        gff_file (str, optional): Path to GFF file for calculating intron lengths

    Returns:
        pd.DataFrame: Fused hits dataframe with composite_score column added
    """
    # Create lookup dictionary for control bitscores
    control_bitscore_dict = df_control.set_index("protein")["bit_score"].to_dict()

    # Step 1: Calculate intron lengths if GFF file provided
    if gff_file:
        intron_lengths_dict = calculate_intron_lengths(df_fused, gff_file)
        # Add intron length column (in base pairs)
        df_fused["theorized_intron_length_bp"] = df_fused.apply(
            lambda row: intron_lengths_dict.get((row['gene_1'], row['gene_2']), None),
            axis=1
        )
    else:
        df_fused["theorized_intron_length_bp"] = None

    # Step 2: Extract organism names from subject_title and count unique organisms per GENE PAIR
    # Extract organism name from subject_title (format: "protein_id description [Organism name]")
    df_fused["organism"] = df_fused["subject_title"].str.extract(r'\[([^\]]+)\]$')[0]

    # Count unique organisms for each GENE PAIR (not protein isoform pair)
    # This way different isoforms of the same gene pair count as one
    organism_counts = df_fused.groupby(["gene_1", "gene_2"])["organism"].nunique().to_dict()

    # Get max organism count for normalization
    max_organism_count = max(organism_counts.values()) if organism_counts else 1

    # Step 3: Define scoring function that will be applied to each row
    def calculate_row_score(row):
        # Get control bitscores
        control_bs_1 = control_bitscore_dict.get(row["product_1"], None)
        control_bs_2 = control_bitscore_dict.get(row["product_2"], None)

        # Get organism count for this gene pair (not isoform pair)
        org_count = organism_counts.get((row["gene_1"], row["gene_2"]), 1)

        # Calculate composite score
        score = scoring.calculate_composite_score(
            query_coverage=row["query_coverage"],
            fused_bitscore=row["bit_score"],
            control_bitscore_1=control_bs_1,
            control_bitscore_2=control_bs_2,
            start_query=row["start_of_alignment_in_query"],
            end_query=row["end_of_alignment_in_query"],
            gene_1_len=row["gene_1_len"],
            organism_hit_count=org_count,
            max_organism_count=max_organism_count,
            evalue=row["expected_value"],
            percent_identity=row["percentage_of_identical_matches"]
        )

        return score

    # Step 3: Apply scoring function to each row
    df_fused["composite_score"] = df_fused.apply(calculate_row_score, axis=1)

    # Step 4: Add organism_count column (use gene pair lookup)
    df_fused["organism_count"] = df_fused.apply(
        lambda row: organism_counts.get((row["gene_1"], row["gene_2"]), None),
        axis=1
    )

    # Step 5: Add _aa suffix to gene length columns and drop legacy columns
    columns_to_drop = ["organism"]  # Start with temporary organism column

    if "gene_1_len" in df_fused.columns:
        df_fused["gene_1_len_aa"] = df_fused["gene_1_len"]
        columns_to_drop.append("gene_1_len")
    if "gene_2_len" in df_fused.columns:
        df_fused["gene_2_len_aa"] = df_fused["gene_2_len"]
        columns_to_drop.append("gene_2_len")
    if "fused_gene_len" in df_fused.columns:
        df_fused["fused_gene_len_aa"] = df_fused["fused_gene_len"]
        columns_to_drop.append("fused_gene_len")

    # Step 5.5: Add alignment position ranges (both normalized and absolute)
    # Normalized ranges (0-1 scale, rounded to 4 decimal places)
    if "start_of_alignment_in_query" in df_fused.columns and "query_length" in df_fused.columns:
        start_norm = (df_fused["start_of_alignment_in_query"] / df_fused["query_length"]).round(4)
        end_norm = (df_fused["end_of_alignment_in_query"] / df_fused["query_length"]).round(4)
        df_fused["alignment_range_in_query_norm"] = start_norm.astype(str) + "-" + end_norm.astype(str)

    if "start_of_alignment_in_subject" in df_fused.columns and "subject_length" in df_fused.columns:
        start_norm = (df_fused["start_of_alignment_in_subject"] / df_fused["subject_length"]).round(4)
        end_norm = (df_fused["end_of_alignment_in_subject"] / df_fused["subject_length"]).round(4)
        df_fused["alignment_range_in_subject_norm"] = start_norm.astype(str) + "-" + end_norm.astype(str)

    # Absolute ranges (amino acid positions)
    if "start_of_alignment_in_query" in df_fused.columns and "end_of_alignment_in_query" in df_fused.columns:
        df_fused["alignment_range_in_query"] = df_fused["start_of_alignment_in_query"].astype(str) + "-" + df_fused["end_of_alignment_in_query"].astype(str)

    if "start_of_alignment_in_subject" in df_fused.columns and "end_of_alignment_in_subject" in df_fused.columns:
        df_fused["alignment_range_in_subject"] = df_fused["start_of_alignment_in_subject"].astype(str) + "-" + df_fused["end_of_alignment_in_subject"].astype(str)

    # Step 6: Drop individual start/end columns (now replaced by ranges)
    columns_to_drop.extend([
        "start_of_alignment_in_query",
        "end_of_alignment_in_query",
        "start_of_alignment_in_subject",
        "end_of_alignment_in_subject"
    ])

    # Step 7: Sort by score (descending) for easier analysis
    df_fused = df_fused.sort_values(by="composite_score", ascending=False).reset_index(drop=True)

    # Remove temporary columns (organism, legacy length columns, and individual position columns)
    df_fused = df_fused.drop(columns=columns_to_drop)

    return df_fused

def visualize_genome_structure_issues(gene_lookup, fragmented_df, control_losses, output_folder):
    """
    Create visual diagrams for genes with genome map issues.
    Shows the gene parts and 3 neighboring genes on each side in the genome graph.
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    import networkx as nx

    # Helper function to find predecessors (genes that point to this gene)
    def find_predecessors(target_gene_id, gene_lookup, max_depth=3):
        """Find up to max_depth genes that have target as a neighbor."""
        predecessors = []
        for gene_id, node in gene_lookup.items():
            if node.neighbors:
                for neighbor in node.neighbors:
                    if neighbor.gene_id == target_gene_id:
                        predecessors.append(node)
                        break
        return predecessors

    # Helper function to get N genes upstream
    def get_upstream_genes(gene_id, gene_lookup, n=3):
        """Get up to n genes upstream by traversing predecessors."""
        upstream = []
        # Get initial node
        start_node = gene_lookup.get(gene_id)
        if not start_node:
            return upstream

        current_nodes = [start_node]
        visited = set([start_node])  # Track GeneNode objects to handle duplicate gene names

        for level in range(n):
            next_level = []
            for curr_node in current_nodes:
                preds = find_predecessors(curr_node.gene_id, gene_lookup)
                for pred in preds:
                    if pred not in visited:
                        next_level.append(pred)
                        visited.add(pred)
                        upstream.append((pred.gene_id, level + 1))
            current_nodes = next_level
            if not current_nodes:
                break
        return upstream

    # Helper function to get N genes downstream
    def get_downstream_genes(gene_id, gene_lookup, n=3):
        """Get up to n genes downstream by following neighbors."""
        downstream = []
        start_node = gene_lookup.get(gene_id)
        if not start_node:
            return downstream

        current_nodes = [start_node]
        visited = set([start_node])  # Track GeneNode objects to handle duplicate gene names

        for level in range(n):
            next_level = []
            for curr_node in current_nodes:
                if curr_node.neighbors:
                    for neighbor in curr_node.neighbors:
                        if neighbor not in visited:
                            next_level.append(neighbor)
                            visited.add(neighbor)
                            downstream.append((neighbor.gene_id, level + 1))
            current_nodes = next_level
            if not current_nodes:
                break
        return downstream

    # Create output directory
    viz_dir = f"{output_folder}/genome_structure_issues"
    os.makedirs(viz_dir, exist_ok=True)

    # Process missing_from_map genes
    for gene_id in control_losses['missing_from_map']:
        part1_gene = f"{gene_id}_part1"
        part2_gene = f"{gene_id}_part2"

        fig, ax = plt.subplots(figsize=(20, 12))
        G = nx.DiGraph()

        # Check which parts exist
        part1_exists = part1_gene in gene_lookup
        part2_exists = part2_gene in gene_lookup

        # Collect all genes to visualize
        genes_to_show = set()

        # Add part1 and its context
        if part1_exists:
            genes_to_show.add(part1_gene)
            # Get 3 upstream and 3 downstream
            upstream = get_upstream_genes(part1_gene, gene_lookup, n=3)
            downstream = get_downstream_genes(part1_gene, gene_lookup, n=3)
            for g, _ in upstream:
                genes_to_show.add(g)
            for g, _ in downstream:
                genes_to_show.add(g)

        # Add part2 and its context
        if part2_exists:
            genes_to_show.add(part2_gene)
            upstream = get_upstream_genes(part2_gene, gene_lookup, n=3)
            downstream = get_downstream_genes(part2_gene, gene_lookup, n=3)
            for g, _ in upstream:
                genes_to_show.add(g)
            for g, _ in downstream:
                genes_to_show.add(g)

        # Build graph
        for gene_name in genes_to_show:
            if gene_name in gene_lookup:
                node = gene_lookup[gene_name]
                G.add_node(gene_name)
                if node.neighbors:
                    for neighbor in node.neighbors:
                        if neighbor.gene_id in genes_to_show:
                            G.add_edge(gene_name, neighbor.gene_id)

        # Layout the graph
        pos = nx.spring_layout(G, k=2, iterations=50)

        # Draw nodes with colors
        node_colors = []
        for node in G.nodes():
            if node == part1_gene and part1_exists:
                node_colors.append('lightgreen' if part2_exists else 'lightyellow')
            elif node == part2_gene and part2_exists:
                node_colors.append('lightgreen')
            elif node in [part1_gene, part2_gene]:
                node_colors.append('lightcoral')  # Missing
            else:
                node_colors.append('lightblue')  # Context genes

        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=3000, ax=ax, node_shape='s')
        nx.draw_networkx_edges(G, pos, edge_color='gray', arrows=True, arrowsize=20, ax=ax,
                              connectionstyle='arc3,rad=0.1', width=2)

        # Draw labels (shortened for readability)
        labels = {n: n[:20] + '...' if len(n) > 20 else n for n in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels, font_size=8, font_weight='bold', ax=ax)

        # Add text annotations for part1 and part2
        if part1_gene in pos:
            ax.annotate('PART1', xy=pos[part1_gene], xytext=(0, -40),
                       textcoords='offset points', ha='center',
                       bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.7),
                       arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
        if part2_gene in pos:
            ax.annotate('PART2', xy=pos[part2_gene], xytext=(0, -40),
                       textcoords='offset points', ha='center',
                       bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.7),
                       arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

        # Legend
        green_patch = mpatches.Patch(color='lightgreen', label='Part exists in map')
        red_patch = mpatches.Patch(color='lightcoral', label='Part missing from map')
        blue_patch = mpatches.Patch(color='lightblue', label='Context genes (±3)')
        yellow_patch = mpatches.Patch(color='lightyellow', label='Part1 (part2 missing)')
        ax.legend(handles=[green_patch, yellow_patch, red_patch, blue_patch], loc='upper right', fontsize=10)

        ax.set_title(f"Genome Structure: {gene_id}\nIssue: Part(s) missing from genome map\nShowing ±3 genes context",
                    fontsize=14, weight='bold')
        ax.axis('off')

        plt.tight_layout()
        plt.savefig(f"{viz_dir}/{gene_id}_missing_from_map.png", dpi=150, bbox_inches='tight')
        plt.close()

    # Process broken_neighbors genes
    for gene_id in control_losses['broken_neighbors']:
        part1_gene = f"{gene_id}_part1"
        part2_gene = f"{gene_id}_part2"

        fig, ax = plt.subplots(figsize=(20, 12))
        G = nx.DiGraph()

        # Collect all genes to visualize
        genes_to_show = set()

        # Add part1 and its context
        genes_to_show.add(part1_gene)
        upstream = get_upstream_genes(part1_gene, gene_lookup, n=3)
        downstream = get_downstream_genes(part1_gene, gene_lookup, n=3)
        for g, _ in upstream:
            genes_to_show.add(g)
        for g, _ in downstream:
            genes_to_show.add(g)

        # Add part2 and its context
        genes_to_show.add(part2_gene)
        upstream = get_upstream_genes(part2_gene, gene_lookup, n=3)
        downstream = get_downstream_genes(part2_gene, gene_lookup, n=3)
        for g, _ in upstream:
            genes_to_show.add(g)
        for g, _ in downstream:
            genes_to_show.add(g)

        # Build graph
        for gene_name in genes_to_show:
            if gene_name in gene_lookup:
                node = gene_lookup[gene_name]
                G.add_node(gene_name)
                if node.neighbors:
                    for neighbor in node.neighbors:
                        if neighbor.gene_id in genes_to_show:
                            G.add_edge(gene_name, neighbor.gene_id)

        # Layout the graph
        pos = nx.spring_layout(G, k=2, iterations=50)

        # Draw nodes with colors
        node_colors = []
        for node in G.nodes():
            if node == part1_gene:
                node_colors.append('lightyellow')  # Part1 (source)
            elif node == part2_gene:
                node_colors.append('lightcoral')  # Part2 (not a neighbor)
            else:
                node_colors.append('lightblue')  # Context genes

        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=3000, ax=ax, node_shape='s')
        nx.draw_networkx_edges(G, pos, edge_color='gray', arrows=True, arrowsize=20, ax=ax,
                              connectionstyle='arc3,rad=0.1', width=2)

        # Draw labels
        labels = {n: n[:20] + '...' if len(n) > 20 else n for n in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels, font_size=8, font_weight='bold', ax=ax)

        # Add annotations
        if part1_gene in pos:
            ax.annotate('PART1', xy=pos[part1_gene], xytext=(0, -40),
                       textcoords='offset points', ha='center',
                       bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.7),
                       arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
        if part2_gene in pos:
            ax.annotate('PART2\n(NOT NEIGHBOR)', xy=pos[part2_gene], xytext=(0, -50),
                       textcoords='offset points', ha='center',
                       bbox=dict(boxstyle='round,pad=0.5', fc='red', alpha=0.7),
                       arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

        # Draw broken connection if both parts exist
        if part1_gene in pos and part2_gene in pos:
            x1, y1 = pos[part1_gene]
            x2, y2 = pos[part2_gene]
            ax.plot([x1, x2], [y1, y2], 'r--', linewidth=3, alpha=0.5, label='Expected connection (broken)')
            mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2
            ax.text(mid_x, mid_y, 'X', fontsize=30, color='red', weight='bold', ha='center', va='center')

        # Legend
        yellow_patch = mpatches.Patch(color='lightyellow', label='Part1 (source)')
        red_patch = mpatches.Patch(color='lightcoral', label='Part2 (not a neighbor)')
        blue_patch = mpatches.Patch(color='lightblue', label='Context genes (±3)')
        ax.legend(handles=[yellow_patch, red_patch, blue_patch], loc='upper right', fontsize=10)

        ax.set_title(f"Genome Structure: {gene_id}\nIssue: Part2 is not a neighbor of Part1\nShowing ±3 genes context",
                    fontsize=14, weight='bold')
        ax.axis('off')

        plt.tight_layout()
        plt.savefig(f"{viz_dir}/{gene_id}_broken_neighbors.png", dpi=150, bbox_inches='tight')
        plt.close()

    print(f"\nGenome structure visualizations saved to: {viz_dir}/")
    print(f"  - {len(control_losses['missing_from_map'])} missing_from_map diagrams")
    print(f"  - {len(control_losses['broken_neighbors'])} broken_neighbors diagrams")

# Wrapping the main script code in main lets us use the other functions in other scripts without calling the whole thing.
if __name__ == "__main__":
    # Argument method validators
    def valid_file(filepath):
        """Validates that the file exists and is readable"""
        if not os.path.isfile(filepath):
            raise argparse.ArgumentTypeError(f"File '{filepath}' does not exist")
        return filepath

    def valid_faa_file(filepath):
        """Validates .faa proteome file"""
        filepath = valid_file(filepath)
        if not filepath.endswith('.faa'):
            raise argparse.ArgumentTypeError("Proteome file must be a .faa file")
        return filepath

    def valid_gff_file(filepath):
        """Validates .gff annotation file"""
        filepath = valid_file(filepath)
        if not filepath.endswith('.gff'):
            raise argparse.ArgumentTypeError("Annotation file must be a .gff file")
        return filepath

    def valid_dmnd_file(filepath):
        """Validates .dmnd database file"""
        filepath = valid_file(filepath)
        if not filepath.endswith('.dmnd'):
            raise argparse.ArgumentTypeError("Database file must be a .dmnd file")
        return filepath

    def positive_int(value):
        """Validates positive integer"""
        ivalue = int(value)
        if ivalue < 1:
            raise argparse.ArgumentTypeError(f"{value} must be at least 1")
        return ivalue

    def percent_identity(value):
        """Validates percent identity is between 0.0 and 100.0"""
        fvalue = float(value)
        if not 0.0 <= fvalue <= 100.0:
            raise argparse.ArgumentTypeError(f"Identity cutoff must be between 0.0 and 100.0, got {value}")
        return fvalue

    def min_gene_count(value):
        """Validates minimum gene count filter"""
        ivalue = int(value)
        if ivalue < 1:
            raise argparse.ArgumentTypeError(f"Gene count filter must be at least 1, got {value}")
        return ivalue

    def valid_sensitivity(value):
        """Validates diamond sensitivity parameter"""
        valid_modes = ["fast", "mid-sensitive", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"]
        if value not in valid_modes:
            raise argparse.ArgumentTypeError(
                f"Invalid sensitivity mode '{value}'. Must be one of: {', '.join(valid_modes)}"
            )
        return value
    
    # Argparse setup
    parser = argparse.ArgumentParser(prog="Genome Misannotation Checker")

    # REQUIRED PARAMS

    parser.add_argument('-p', "--proteome",
                        help="The filename of the .faa file containing the proteome of the organism of interest. Required input.",
                        required=True,
                        type=valid_faa_file)

    parser.add_argument('-a', '--organism_annotation',
                        help="Input the organism's annotation features in a gff format. This should be a RefSeq annotation. Required input",
                        required=True,
                        type=valid_gff_file)

    parser.add_argument('-db', '--database',
                        help="Input the local reference protein database in a dmnd format. Required input.",
                        required=True,
                        type=valid_dmnd_file)
    parser.add_argument('-o', '--output',
                        help="The output folder path. Required input.",
                        required=True)

    # OPTIONAL PARAMS

    parser.add_argument('-t', '--num_threads',
                        help="Input the number of threads that you would like to use. By default, half of your available threads will be used.",
                        default=get_default_num_threads(),
                        type=positive_int)

    parser.add_argument('-i', '--identity_cutoff',
                        help="Input the percent identity cutoff you would like to use for filtering of the fused gene alignments. The percent identity is the percentage of identical amino acids between two sequences at the same alignment positions. The default is 0.00.",
                        default=0.00,
                        type=percent_identity)

    parser.add_argument('-xf', "--xfilter",
                        help="Filters out all chromosomes and contigs that have less than the specified number of genes (minimum: 1)",
                        type=min_gene_count,
                        default=5)

    parser.add_argument('-ds', '--diamond_sensitivity',
                        help="Sensitivity mode parameter for the Diamond alignment tool. Valid options: fast, mid-sensitive, sensitive, more-sensitive, very-sensitive, ultra-sensitive",
                        type=valid_sensitivity,
                        default=None)
    
    parser.add_argument('-lo', '--longest_only',
                        help="Choose whether each isoform or only the longest isoform for each gene is fused. Default is False meaning that each isoform is fused.",
                        action="store_true")

    parser.add_argument('-n', '--organism_name',
                        help="Scientific name of the organism (e.g., 'Drosophila melanogaster')",
                        type=str,
                        default=None)
    parser.add_argument('-d', '--taxon-id',
                        help="Taxon ID of the organism. Use to filter self-hits during DIAMOND",
                        type=str,
                        default=None)

    parser.add_argument('-c', '--control-report',
                        help="Path to fragmented genes report CSV for synthetic control tracking (enables DEBUG_CONTROL mode)",
                        type=str,
                        default=None)

    args = parser.parse_args()



    annotation = args.organism_annotation
    taxon_id = detect_taxon_id(annotation)

    # Use command-line organism name if provided, otherwise try to detect
    if args.organism_name:
        organism_name = args.organism_name
        print(f"Using provided organism name: {organism_name}")
    else:
        organism_name = detect_organism_name(annotation)
    if organism_name is None:
        print("WARNING: Organism name not detected in the annotation file.")
        exit(1)

    if taxon_id is None:
        print("WARNING: Organism name not detected in the annotation file.")
        if args.taxon_id:
            taxon_id = args.taxon_id
        else:
            print("Please provide organism taxon ID with -n flag (e.g., -n 7010)")
            exit(1)


    # Set max target seqs higher to ensure good hits even after taxon filtering
    max_target_seqs = 200

    diamond_path = shutil.which("diamond")
    if diamond_path is None:
        print("diamond not found in PATH, please check your installation.")
        exit()

    db = args.database
    proteome_file = args.proteome
    num_threads = str(args.num_threads)
    ident_cutoff = float(args.identity_cutoff) * 100
    min_geneperx_threshold = args.xfilter
    gff_file = args.organism_annotation
    output_folder = args.output
    if args.diamond_sensitivity:
        diamond_sensitivity = "--" + args.diamond_sensitivity
    else:
        diamond_sensitivity = args.diamond_sensitivity
    longest_only = args.longest_only
    

    # Create output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Created output folder: {output_folder}")

    # Storing the command-line arguments inputted by the user in a log file for future reference by the user.
    with open(f"{output_folder}/args.log", 'w') as log_file:
        log_file.write("Command-line input:\n")
        log_file.write(" ".join(sys.argv) + "\n")

    # A fasta to store all positive gene parts
    genome_faa_path = f"{output_folder}/genome_wide_positive_hits.faa"
    
    # Ensure the genome-wide file starts empty
    if os.path.exists(genome_faa_path):
        os.remove(genome_faa_path)
        
    # A db to store all positive hits for the entire genome
    fused_hits_genome_df_list = []

    genomeMap = build_genomemap(organism_name, proteome_file, gff_file, min_geneperx_threshold)

    # Enable DEBUG_CONTROL if control report is provided
    if args.control_report:
        if os.path.exists(args.control_report):
            DEBUG_CONTROL = True
            fragmented_df = pd.read_csv(args.control_report)
            print(f"DEBUG_CONTROL enabled. Tracking {len(fragmented_df)} synthetic control genes from: {args.control_report}")
        else:
            print(f"WARNING: Control report file not found at {args.control_report}. Skipping control gene tracking.")
            DEBUG_CONTROL = False

    if DEBUG_CONTROL:

        # Clear debug tracking variables from any previous runs
        tested_genes.clear()
        best_diamond_hits.clear()
        for key in control_losses:
            control_losses[key].clear()
        # Build lookup dictionary for O(1) gene access
        gene_lookup = {}
        for chrom, strands in genomeMap.chromosomes.items():
            for strand_name, head_node in strands.items():
                current = head_node
                while current:
                    gene_lookup[current.gene_id] = current
                    current = current.neighbors[0] if current.neighbors else None

        # Track losses during map construction using optimized lookups
        for gene_id in fragmented_df['gene_name']:
            part1_gene = f"{gene_id}_part1"
            part2_gene = f"{gene_id}_part2"

            # Check if parts exist in map (O(1) lookup)
            if part1_gene not in gene_lookup or part2_gene not in gene_lookup:
                control_losses['missing_from_map'].add(gene_id)
                continue

            # Check if part2 is a neighbor of part1 (O(1) lookup)
            part1_node = gene_lookup[part1_gene]
            if part1_node.neighbors:
                neighbor_ids = [n.gene_id for n in part1_node.neighbors]
                if part2_gene not in neighbor_ids:
                    control_losses['broken_neighbors'].add(gene_id)
            else:
                control_losses['broken_neighbors'].add(gene_id)

    # Initializing a progress bar for tracking the program's status.
    total_steps = (len(genomeMap.chromosomes.keys()) * 4) + 1 
    with tqdm(total=total_steps, unit="step") as pbar:

        # Going into each strand of each chromosome, getting each protein's sequence, fusing neighboring genes, then using DIAMOND to check for misannotations.
        for chrom in genomeMap.chromosomes.keys():
            for strand in ["+", "-"]:
                # Avoiding use of + in file names since it is a no no in the rubric :)
                if strand == "+":
                    strand_name = "plus"
                elif strand == "-":
                    strand_name = "neg"

                output_prefix = f"{output_folder}/{chrom}/{strand_name}"  # Setting up the output file path.
                if not os.path.exists(output_prefix):
                    os.makedirs(output_prefix)

                pbar.set_description(f"Processing {chrom} {strand_name} strand")
                headNode = genomeMap.get_head_node(chrom, strand)
                if headNode is None:
                    print(f"{chrom}'s {strand_name} strand does not contain any protein-coding genes.")
                    pbar.update(2) 
                else:
                    # Process fused proteins first
                    temp_faa_filepath, fused_metadata_df = chromosome_processor_fused_allisoforms(chrom, strand, headNode, output_prefix, longest_only)

                    # Check if there are any fusions to process
                    if os.path.getsize(temp_faa_filepath) > 0 and not fused_metadata_df.empty:
                        pbar.set_description(f"Running DIAMOND on {chrom} {strand_name} strand's fused genes")
                        fused_chrom_hits_df, fused_diamond_df = fused_diamond_alignment(chrom, strand_name, diamond_path, db, temp_faa_filepath,
                                                                                fused_metadata_df, num_threads, ident_cutoff, output_prefix, diamond_sensitivity, taxon_id, max_target_seqs)
                        fused_hits_genome_df_list.append(fused_chrom_hits_df)
                        # fused_diamond_df.to_csv(f"{output_prefix}/test_df.csv")
                    else:
                        pbar.set_description(f"Skipping {chrom} {strand_name} strand (no fusions to process)")
                        # Clean up empty temp file
                        if os.path.exists(temp_faa_filepath):
                            os.remove(temp_faa_filepath)
                    pbar.update(1)

                    if DEBUG_CONTROL:
                        tested_genes.clear()  # Reset for next strand

                    # Process unfused proteins second
                    temp_faa_filepath = chromosome_processor_unfused_allisoforms(chrom, strand, headNode, output_prefix, fused_chrom_hits_df, longest_only)
                    if os.path.exists(temp_faa_filepath) and os.path.getsize(temp_faa_filepath) > 0:
                        with open(temp_faa_filepath, "r") as temp_faa, open(genome_faa_path, "a") as genome_faa:
                            genome_faa.write(temp_faa.read())
                    pbar.update(1)
                    
        #The control is to see if the alignment score increases in the fused gene vs the unfused genes for any overlapping hits in the fused genes
        pbar.set_description("Running control DIAMOND on identified gene parts")
        control_hits = unfused_diamond_alignment(chrom, strand_name, diamond_path, db, genome_faa_path, num_threads, output_prefix, diamond_sensitivity, taxon_id, max_target_seqs)
        pbar.update(1)

    
    
    # Compiling all unique gene alignments to have one dataframe that contains genome-wide results.
    if fused_hits_genome_df_list:
        fused_hits = pd.concat(fused_hits_genome_df_list, ignore_index=True)
        # Organism filtering is handled by --taxon-exclude during BLAST, no need for post-filtering

        # Calculate composite scores for all fused hits
        print("Calculating composite scores for fused gene hits...")
        fused_hits_scored = calculate_scores_for_hits(fused_hits, control_hits, gff_file)
        highest_composite_score = fused_hits_scored.groupby("fused_protein")["composite_score"].idxmax()
        fused_hits_scored_filtered = fused_hits_scored.loc[highest_composite_score].reset_index(drop=True)
        fused_hits_scored_filtered = fused_hits_scored_filtered.sort_values(by="composite_score", ascending=False).reset_index(drop=True)

        # Enrich best_diamond_hits with composite_score and organism_count
        if DEBUG_CONTROL and best_diamond_hits:
            for gene_id, hit_info in best_diamond_hits.items():
                fused_protein = hit_info.get('fused_protein')
                if fused_protein and not pd.isna(fused_protein):
                    # Find this fused_protein in the scored dataframe
                    matching_rows = fused_hits_scored[fused_hits_scored['fused_protein'] == fused_protein]
                    if not matching_rows.empty:
                        # Get the row with highest composite score for this fused_protein
                        best_score_idx = matching_rows['composite_score'].idxmax()
                        best_score_row = matching_rows.loc[best_score_idx]

                        # Add composite score and organism count to the hit info
                        hit_info['composite_score'] = best_score_row.get('composite_score', '')
                        hit_info['organism_count'] = best_score_row.get('organism_count', '')

        # Generate all output formats (CSV, TSV, Excel)
        print("Generating output files...")
        output_formatter.generate_all_outputs(
            df_fused=fused_hits_scored_filtered,
            df_control=control_hits,
            organism_name=organism_name,
            output_folder=output_folder
        )

    else:
        print("No fused genes found across all chromosomes/strands.")
        # Create empty results file with proper headers
        pd.DataFrame(columns=["fused_gene", "fused_product", "fused_gene_len", "gene_1", "product_1",
                              "gene_1_len", "gene_2", "product_2", "gene_2_len"]).to_csv(f"{output_folder}/full_statistics_genome_results.csv")

    if DEBUG_CONTROL and fragmented_df is not None:
        total_genes = fragmented_df['gene_name'].nunique()

        # Build detailed gene-level report
        gene_reports = []

        for gene_id in fragmented_df['gene_name'].unique():
            # Determine filter step
            if gene_id in control_losses['missing_from_map']:
                filter_step = 'missing_from_map'
                notes = 'Part1 or Part2 not found in genome map'
            elif gene_id in control_losses['broken_neighbors']:
                filter_step = 'broken_neighbors'
                notes = 'Part2 is not a neighbor of Part1 in genome graph'
            elif gene_id in control_losses['missing_proteins']:
                filter_step = 'missing_proteins'
                notes = 'Part1 or Part2 has no protein isoforms'
            elif gene_id in control_losses['no_diamond_hits']:
                filter_step = 'no_diamond_hits'
                notes = 'DIAMOND returned no alignments'
            elif gene_id in control_losses['failed_alignment_overlap']:
                filter_step = 'failed_alignment_overlap'
                notes = 'All hits failed to span gene boundary by 10 AA'
            elif gene_id in control_losses['failed_query_coverage']:
                filter_step = 'failed_query_coverage'
                notes = 'All hits had <50% query coverage'
            elif gene_id in control_losses['failed_percent_identity']:
                filter_step = 'failed_percent_identity'
                notes = 'All hits had <50% sequence identity'
            elif gene_id in control_losses['failed_uncharacterized']:
                filter_step = 'failed_uncharacterized'
                notes = 'All hits were to uncharacterized proteins'
            elif gene_id in control_losses['passed_all_filters']:
                # Skip genes that passed - we only want to report problems
                continue
            else:
                filter_step = 'not_tested'
                notes = 'Gene was not tested (not encountered during processing)'
                control_losses['not_tested'].add(gene_id)

            # Get best DIAMOND hit info if available
            best_hit = best_diamond_hits.get(gene_id, {})

            gene_reports.append({
                'gene_name': gene_id,
                'filter_step': filter_step,
                'fused_protein': best_hit.get('fused_protein', ''),
                'composite_score': best_hit.get('composite_score', ''),
                'organism_count': best_hit.get('organism_count', ''),
                'subject_id': best_hit.get('subject_id', ''),
                'subject_title': best_hit.get('subject_title', ''),
                'bitscore': best_hit.get('bitscore', ''),
                'qcoverage': best_hit.get('qcoverage', ''),
                'pident': best_hit.get('pident', ''),
                'evalue': best_hit.get('evalue', ''),
                'alignment_start': best_hit.get('alignment_start', ''),
                'alignment_end': best_hit.get('alignment_end', ''),
                'gene_1_len_aa': best_hit.get('gene_1_len', ''),
                'gene_2_len_aa': best_hit.get('gene_2_len', ''),
                'fused_gene_len_aa': best_hit.get('fused_gene_len', ''),
                'overlap_status': best_hit.get('overlap_status', ''),
                'part1_overlap_aa': best_hit.get('part1_overlap_depth', ''),
                'part2_overlap_aa': best_hit.get('part2_overlap_depth', ''),
                'notes': notes
            })

        # Create DataFrame and save detailed report
        detailed_df = pd.DataFrame(gene_reports)
        detailed_df.to_csv(f"{output_folder}/control_gene_tracking_detailed.csv", index=False)

        # Build summary from detailed report to avoid double-counting genes in multiple categories
        # Add passed_all_filters count separately since it's not in detailed report
        summary_rows = []

        # Count genes per filter_step from detailed report (no overlaps)
        if not detailed_df.empty:
            category_counts = detailed_df.groupby('filter_step')['gene_name'].nunique()
            for category, count in category_counts.items():
                pct = 100 * count / total_genes if total_genes > 0 else 0
                summary_rows.append({
                    'filter_step': category,
                    'genes_lost': count,
                    'percentage': round(pct, 2)
                })

        # Add passed_all_filters from control_losses (these genes are skipped in detailed report)
        passed_count = len(control_losses['passed_all_filters'])
        passed_pct = 100 * passed_count / total_genes if total_genes > 0 else 0
        summary_rows.append({
            'filter_step': 'passed_all_filters',
            'genes_lost': passed_count,
            'percentage': round(passed_pct, 2)
        })

        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_csv(f"{output_folder}/control_gene_tracking_summary.csv", index=False)

        # Print summary table
        print("\n=== SYNTHETIC CONTROL GENE TRACKING SUMMARY ===")
        print(f"Total fragmented genes: {total_genes}")
        print(f"\n{'Filter Step':<30} {'Genes Lost':<12} {'Percentage':<12}")
        print("-" * 55)
        for _, row in summary_df.iterrows():
            print(f"{row['filter_step']:<30} {row['genes_lost']:<12} {row['percentage']:>6.2f}%")
        print("=" * 55)
        print(f"Summary: {output_folder}/control_gene_tracking_summary.csv")
        print(f"Detailed: {output_folder}/control_gene_tracking_detailed.csv")

        # Generate visualizations for genome structure issues
        if control_losses['missing_from_map'] or control_losses['broken_neighbors']:
            visualize_genome_structure_issues(gene_lookup, fragmented_df, control_losses, output_folder)

    print("Your results are ready.")
