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
DEBUG_CONTROL = True
tested_genes = set()  # Global set to track fragmented genes being tested (per strand)
fragmented_df = None  # Global to store fragmented genes report
debug_dfs_saved = False  # Track if we've saved debug dataframes
# Cumulative tracking across all strands
cumulative_tested = set()  # All genes tested across all strands
cumulative_passed = set()  # All genes that passed filters
cumulative_filtered = set()  # All genes that were filtered out
    
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
    visited = set()  # Track visited gene_ids to avoid duplicates when paths converge
    queue = [headNode]  # Initialize queue with head node (BFS uses FIFO)
    
    # Creating sets of product_1 and product_2 to only write out fused gene parts to the faa
    product1_set = set(fused_chrom_hits_df['product_1'])
    product2_set = set(fused_chrom_hits_df['product_2'])

    while queue:
        currentNode = queue.pop(0)  # Dequeue: pop from front (BFS: process level-by-level)

        # Skip if we've already processed this gene (handles convergence)
        if currentNode.gene_id in visited:
            continue
        visited.add(currentNode.gene_id)

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
    visited_nodes = set()  # Track visited nodes for graph traversal
    visited_edges = set()  # Track processed edges to avoid duplicate fusions
    queue = [headNode]

    while queue:
        currentNode = queue.pop(0)  # BFS: pop from front

        # Skip if we've already visited this node
        if currentNode.gene_id in visited_nodes:
            continue
        visited_nodes.add(currentNode.gene_id)

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
                                    print(f"    ⚠ Gene {gene_base}: part1 has {part1_protein_count} proteins, part2 has {part2_protein_count} proteins")
                                    if part1_protein_count == 0:
                                        print(f"      Part1 gene_id: {currentNode.gene_id}")
                                    if part2_protein_count == 0:
                                        print(f"      Part2 gene_id: {neighbor.gene_id}")
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

    if DEBUG_CONTROL:
        # Check if tested genes exist in diamond_df BEFORE merge
        if tested_genes:
            # Get all fused_protein IDs that correspond to tested genes
            expected_fused_ids = set()
            for _, row in fused_metadata_df.iterrows():
                if pd.notna(row.get('gene_1')) and pd.notna(row.get('gene_2')):
                    if row['gene_1'].endswith('_part1') and row['gene_2'].endswith('_part2'):
                        gene_base = row['gene_1'][:-6]
                        if row['gene_2'] == f"{gene_base}_part2" and gene_base in tested_genes:
                            expected_fused_ids.add(row['fused_protein'])

            # Check which ones are in diamond_df
            diamond_fused_ids = set(diamond_df['fused_protein'])
            found_in_diamond = expected_fused_ids & diamond_fused_ids
            missing_from_diamond = expected_fused_ids - diamond_fused_ids

            if missing_from_diamond:
                print(f"  ⚠ BEFORE MERGE: {len(missing_from_diamond)}/{len(expected_fused_ids)} tested gene fusions MISSING from DIAMOND results!")
                print(f"    (DIAMOND never returned hits for these)")
                # Extract gene names for readability
                missing_genes = set()
                for fused_id in list(missing_from_diamond)[:5]:
                    for _, row in fused_metadata_df.iterrows():
                        if row['fused_protein'] == fused_id:
                            missing_genes.add(row['gene_1'][:-6])
                            break
                print(f"    Example genes: {list(missing_genes)}")
            else:
                print(f"  ✓ BEFORE MERGE: All {len(expected_fused_ids)} tested gene fusions found in DIAMOND results")

    if DEBUG_CONTROL:
        print(f"  DEBUG: fused_metadata_df shape: {fused_metadata_df.shape}")
        print(f"  DEBUG: diamond_df shape: {diamond_df.shape}")
        print(f"  DEBUG: Unique fused_proteins in metadata: {fused_metadata_df['fused_protein'].nunique()}")
        print(f"  DEBUG: Unique fused_proteins in diamond: {diamond_df['fused_protein'].nunique()}")
        # Check for exact string matches
        metadata_ids = set(fused_metadata_df['fused_protein'])
        diamond_ids = set(diamond_df['fused_protein'])
        print(f"  DEBUG: Overlap in fused_protein IDs: {len(metadata_ids & diamond_ids)}/{len(metadata_ids)}")

    # Merging the diamond results dataframe with the fused protein dataframe. Left join keeps all fusions.
    fused_diamond_df = fused_metadata_df.merge(diamond_df, how='left', on='fused_protein', sort=False)

    if DEBUG_CONTROL:
        print(f"  DEBUG: fused_diamond_df shape AFTER merge: {fused_diamond_df.shape}")
        print(f"  DEBUG: fused_diamond_df columns: {list(fused_diamond_df.columns)}")
        # Check if gene_1 column exists and is correct
        if 'gene_1' in fused_diamond_df.columns:
            print(f"  DEBUG: gene_1 column exists, non-null count: {fused_diamond_df['gene_1'].notna().sum()}")
        else:
            print(f"  DEBUG: WARNING - gene_1 column NOT FOUND! Available columns: {list(fused_diamond_df.columns)}")

        # Save debug dataframes for first strand only
        global debug_dfs_saved
        if not debug_dfs_saved:
            print(f"  DEBUG: Saving dataframes for analysis...")
            fused_metadata_df.to_csv(f"{output_prefix}/debug_fused_metadata.csv", index=False)
            diamond_df.to_csv(f"{output_prefix}/debug_diamond.csv", index=False)
            fused_diamond_df.to_csv(f"{output_prefix}/debug_fused_diamond.csv", index=False)
            # Also save the tested_genes set
            import json
            with open(f"{output_prefix}/debug_tested_genes.json", 'w') as f:
                json.dump(list(tested_genes), f, indent=2)
            print(f"  DEBUG: Saved to {output_prefix}/")
            debug_dfs_saved = True

    if DEBUG_CONTROL:
        # Check if tested genes made it to DIAMOND results
        if tested_genes:
            # Extract gene pairs from fused_diamond_df
            fused_genes_found = set()
            for _, row in fused_diamond_df.iterrows():
                if pd.notna(row.get('gene_1')) and pd.notna(row.get('gene_2')):
                    if row['gene_1'].endswith('_part1') and row['gene_2'].endswith('_part2'):
                        gene_base = row['gene_1'][:-6]
                        if row['gene_2'] == f"{gene_base}_part2":
                            fused_genes_found.add(gene_base)

            missing_from_diamond = tested_genes - fused_genes_found
            if missing_from_diamond:
                print(f"  WARNING: {len(missing_from_diamond)} tested genes missing from fused_diamond_df!")
                print(f"  Examples: {list(missing_from_diamond)[:5]}")
            else:
                print(f"  ✓ All {len(tested_genes)} tested genes found in fused_diamond_df")

    # Only filtering for fused genes that have overlaps with at least 10 AAs, maybe make this user-inputted
    fused_hits_df = fused_diamond_df[
        (fused_diamond_df["start_of_alignment_in_query"] < (fused_diamond_df["gene_1_len"]-10)) &
        (fused_diamond_df["end_of_alignment_in_query"] > (fused_diamond_df["gene_1_len"]+10)) &
        (fused_diamond_df["query_coverage"] >= 50) &
        (fused_diamond_df["percentage_of_identical_matches"] >= 50) &
        (~fused_diamond_df["subject_title"].str.contains("uncharacterized", case=False, na=False))
    ].copy()
    
    if DEBUG_CONTROL:
        # Check if tested genes survived filtering
        if tested_genes:
            # Extract genes from fused_hits_df
            passed_filter = set()
            for _, row in fused_hits_df.iterrows():
                if pd.notna(row.get('gene_1')) and pd.notna(row.get('gene_2')):
                    if row['gene_1'].endswith('_part1') and row['gene_2'].endswith('_part2'):
                        gene_base = row['gene_1'][:-6]
                        if row['gene_2'] == f"{gene_base}_part2":
                            passed_filter.add(gene_base)

            filtered_out = tested_genes - passed_filter

            # Update cumulative tracking
            cumulative_passed.update(passed_filter)
            cumulative_filtered.update(filtered_out)

            if filtered_out:
                print(f"  ⚠ {len(filtered_out)}/{len(tested_genes)} tested genes FILTERED OUT!")
                print(f"  Examples: {list(filtered_out)[:5]}")
            else:
                print(f"  ✓ All {len(tested_genes)} tested genes PASSED filters")
    
    
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
    if DEBUG_CONTROL:
        # Clear debug tracking variables from any previous runs
        tested_genes.clear()
        cumulative_tested.clear()
        cumulative_passed.clear()
        cumulative_filtered.clear()
        filtered_genes = set()
        GENOFRAG_DIR = '/home/davidbellini/OneDrive/gabbiani/TOOLS/GenoFrag/outputs/simulans_new'
        FRAGMENTED_GFF = os.path.join(GENOFRAG_DIR, 'test_simulans_10pct.gff')
        FRAGMENTED_FAA = os.path.join(GENOFRAG_DIR, 'test_simulans_10pct.faa')
        FRAGMENTED_REPORT = os.path.join(GENOFRAG_DIR, 'test_simulans_10pct_report.csv')
        fragmented_df = pd.read_csv(FRAGMENTED_REPORT)
        all_genes_in_map = set()
        for chrom, strands in genomeMap.chromosomes.items():
            for strand_name, head_node in strands.items():
                current = head_node
                while current:
                    all_genes_in_map.add(current.gene_id)
                    current = current.neighbors[0] if current.neighbors else None

        # Check neighbor relationships
        correct = 0
        broken = 0
        missing = 0

        for idx, row in fragmented_df.iterrows():
            gene_id = row['gene_name']
            part1_gene = f"{gene_id}_part1"
            part2_gene = f"{gene_id}_part2"

            # Check if parts exist in map
            if part1_gene not in all_genes_in_map or part2_gene not in all_genes_in_map:
                missing += 1
                filtered_genes.add(gene_id)
                continue

            # Find part1 node and check if part2 is its neighbor
            found_correct_neighbor = False
            for chrom, strands in genomeMap.chromosomes.items():
                for strand_name, head_node in strands.items():
                    current = head_node
                    while current:
                        if current.gene_id == part1_gene:
                            if current.neighbors:
                                neighbor_ids = [n.gene_id for n in current.neighbors]
                                if part2_gene in neighbor_ids:
                                    found_correct_neighbor = True
                            break
                        current = current.neighbors[0] if current.neighbors else None
                if found_correct_neighbor:
                    break

            if found_correct_neighbor:
                correct += 1
            else:
                broken += 1
                filtered_genes.add(gene_id)

        print(f"\ngenoMap relationships:")
        print(f"  Correct: {correct} ({100 * correct / len(fragmented_df):.1f}%)")
        print(f"  Broken: {broken} ({100 * broken / len(fragmented_df):.1f}%)")
        print(f"  Missing from map: {missing} ({100 * missing / len(fragmented_df):.1f}%)")
        print(f"1st FILTER LOST {len(filtered_genes)} genes")

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
                        if tested_genes:
                            print(f"  {chrom} {strand_name}: Tested {len(tested_genes)} fragmented genes this strand")
                            cumulative_tested.update(tested_genes)
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

    if DEBUG_CONTROL:
        print("\nCUMULATIVE DEBUG SUMMARY:")
        print(f"Total fragmented genes tested: {len(cumulative_tested)}")
        print(f"Genes that PASSED all filters: {len(cumulative_passed)} ({100*len(cumulative_passed)/len(cumulative_tested) if cumulative_tested else 0:.1f}%)")
        print(f"Genes that were FILTERED OUT: {len(cumulative_filtered)} ({100*len(cumulative_filtered)/len(cumulative_tested) if cumulative_tested else 0:.1f}%)")
        if cumulative_filtered:
            print(f"Example filtered genes: {list(cumulative_filtered)[:10]}")
        if fragmented_df is not None:
            total_fragmented = fragmented_df['gene_name'].nunique()
            print(f"Total genes in GenoFrag report: {total_fragmented}")
            print(f"Genes in tested set: {len(cumulative_tested)} ({100*len(cumulative_tested)/total_fragmented:.1f}%)")
            not_tested = total_fragmented - len(cumulative_tested)
            if not_tested > 0:
                print(f"⚠ Genes NOT tested (filtered earlier): {not_tested}")

    print("Your results are ready.")
