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

# FLAGS AND GLOBAL VARIABLES
ORGANISM_NAME_DETECTED = False
    
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

def chromosome_processor_unfused_allisoforms(chrom, strand, headNode: GenomeMap.GeneNode, output_prefix, fused_chrom_hits_df):
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

        # Process all protein isoforms for this gene
        for isoformID, isoformAA in currentNode.protein_isoforms.items():
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

def chromosome_processor_fused_allisoforms(chrom, strand, headNode: GenomeMap.GeneNode, output_prefix):
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

                # Create all possible isoform combinations (Cartesian product)
                for isoformID1, isoformAA1 in currentNode.protein_isoforms.items():
                    for isoformID2, isoformAA2 in neighbor.protein_isoforms.items():
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

def isoform_picker(joined_df):
    """Function that goes through the inputted dataframe and analyzes the length of each gene, keeping only the longest gene of multiple 
    entries of the same gene. This is a means to keep only the longest isoform of a gene with multiple isoforms. Returns a dataframe with 
    only one entry per gene."""
    
    joined_df['sequence_length'] = joined_df['sequence'].apply(len) # Get the length of each gene's sequence.
    joined_df = joined_df.loc[joined_df.groupby(['gene'])['sequence_length'].idxmax()] # Keep only the gene with the longest sequence length.
    joined_df = joined_df.drop(columns=['sequence_length'])
    joined_df = joined_df.sort_values(by=['start', 'end'])
    
    return joined_df

def unfused_diamond_alignment(chrom, strand_name, diamond_path, diamond_db, genome_faa_path, num_threads, output_prefix, organism_name, diamond_sensitivity=None):
    """Function used to get the original alignment scores to use as a comparison with the fused gene alignments."""

    diamond_output = f"{output_prefix}/{chrom}_{strand_name}_control_diamond_results.tsv"
    diamond_command = [diamond_path, "blastp", "--db", diamond_db, "--query", genome_faa_path, "--out", diamond_output,
                       "--outfmt", "6", "qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "sstart", "send",
                       "pident", "nident", "mismatch", "evalue", "bitscore", "length", "qcovhsp", "scovhsp", "qtitle", "stitle",
                       "--header", "--evalue", "1e-5", "--threads", num_threads]

    # Add sensitivity parameter if specified
    if diamond_sensitivity:
        diamond_command.append(diamond_sensitivity)

    try:
        result = subprocess.run(diamond_command, check=True, text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"Error during DIAMOND execution for {chrom}, {strand_name} strand.")
        print(f"DIAMOND stderr: {e.stderr}")
        print(f"DIAMOND stdout: {e.stdout}")
        raise

    headers = ["protein", "query_length", "subject_id", "subject_length", "start_of_alignment_in_query", "end_of_alignment_in_query", 
               "start_of_alignment_in_subject", "end_of_alignment_in_subject", "percentage_of_identical_matches", "number_of_identical_matches", 
               "number_of_mismatches", "expected_value", "bit_score", "alignment_length", "query_coverage", "subject_coverage", 
               "query_title", "subject_title"]
    
    diamond_df = pd.read_csv(f"{output_prefix}/{chrom}_{strand_name}_control_diamond_results.tsv", sep="\t", skiprows = 3, names = headers)
    
    # Filtering out matches to own organism
    diamond_df_no_org = diamond_df[~diamond_df["subject_title"].str.contains(organism_name, regex=False, na=False)]
    return diamond_df_no_org
    # keeping entry with the highest qcov and scov
    # diamond_df['coverage_sum'] = diamond_df['query_coverage'] + diamond_df['subject_coverage']
    # diamond_df.to_csv(f"{output_prefix}/{chrom}_{strand_name}_control_diamond_results.tsv", sep="\t")
    # unique_diamond_df = diamond_df.loc[diamond_df.groupby('gene')['coverage_sum'].idxmax()].drop(columns='coverage_sum')
    # unique_diamond_df = unique_diamond_df.loc[unique_diamond_df.groupby('gene')['percentage_of_identical_matches'].idxmax()]
    # unique_diamond_df.to_csv(f"{output_prefix}/{chrom}_{strand_name}_control_unique_diamond_results.tsv", sep="\t", header=True, index=True)
    
def fused_diamond_alignment(chrom, strand_name, diamond_path, db, temp_fasta_path, fused_metadata_df, num_threads, ident_cutoff, output_prefix, diamond_sensitivity=None):
    """Function used to fuse neigboring genes and run the DIAMOND protein alignment script on the fused genes. Returns a dataframe containing unique genes that
    fit the user-inputted filtering criteria."""

    # Creating the diamond output file path and the diamond command that will be ran in the subprocess function.
    diamond_output = f"{output_prefix}/{chrom}_{strand_name}_diamond_results.tsv"
    diamond_command = [diamond_path, "blastp", "--db", db, "--query", temp_fasta_path, "--out", diamond_output,
                       "--outfmt", "6", "qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "sstart", "send",
                       "pident", "nident", "mismatch", "evalue", "bitscore", "length", "qcovhsp", "scovhsp", "qtitle", "stitle",
                       "--header", "--evalue", "1e-5", "--threads", num_threads]

    # Add sensitivity parameter if specified
    if diamond_sensitivity:
        diamond_command.append(diamond_sensitivity)

    # Running the diamond command via subprocess.
    try:
        result = subprocess.run(diamond_command, check=True, text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
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
    # diamond_df["coverage_sum"] = diamond_df["query_coverage"] + diamond_df["subject_coverage"]
    # diamond_df.to_csv(f"{output_prefix}/{chrom}_{strand_name}_diamond_results.tsv", sep="\t")
    # Keep only the row with the maximum coverage sum for each fused gene
    # diamond_df = diamond_df.loc[diamond_df.groupby('fused_gene')['coverage_sum'].idxmax()].drop(columns="coverage_sum")

    # Merging the diamond results dataframe with the fused protein dataframe. Had to use a weird gimmicky solution I found to keep the original order preserved.
    fused_diamond_df = fused_metadata_df.merge(fused_metadata_df.merge(diamond_df, how='outer', on='fused_protein', sort=False))
    # Splitting up the subject_title column if not doing a one organism alignment
    # diamond_df[["subject_XP", "subject_product", "subject_organism"]] = diamond_df["subject_title"].str.extract(r"((?:XP|NP)_[\d\.]+)\s(.+)\s\[(.+)\]")
    # diamond_df = diamond_df.drop(columns=["subject_title"])
    # fused_diamond_df.to_csv(f"{output_prefix}/{chrom}_{strand_name}_max_diamond_results.tsv", sep="\t")
    
    # Only filtering for fused genes that have overlaps with at least 10 AAs, maybe make this user-inputted
    fused_hits_df = fused_diamond_df[(fused_diamond_df["start_of_alignment_in_query"] < (fused_diamond_df["gene_1_len"]-10)) &
                                     (fused_diamond_df["end_of_alignment_in_query"] > (fused_diamond_df["gene_1_len"]+10))].copy()
    
    # Applying the identity cutoff filter inputted by the user as well as filtering to only keep hits with overlaps between the two gene parts.
    # filtered_fused_diamond_df = fused_diamond_df[(fused_diamond_df["start_of_alignment_in_query"] < fused_diamond_df["gene_1_len"]) &
                                                 # (fused_diamond_df["end_of_alignment_in_query"] > fused_diamond_df["gene_1_len"]) & 
                                                 # (fused_diamond_df["percentage_of_identical_matches"] > ident_cutoff)
                                                 # ].copy()
    
    # Creating overlap scores for one's that we know have an overlap for easier viewing for the user. These are the average of each gene's coverage within the alignment.
    # filtered_fused_diamond_df.loc[:, "coverage_gene_1"] = ((filtered_fused_diamond_df["gene_1_len"] - filtered_fused_diamond_df["start_of_alignment_in_query"]) / filtered_fused_diamond_df["gene_1_len"])
    # filtered_fused_diamond_df.loc[:, "coverage_gene_2"] = ((filtered_fused_diamond_df["end_of_alignment_in_query"] - (filtered_fused_diamond_df["gene_1_len"] + 1)) / filtered_fused_diamond_df["gene_2_len"])
    # filtered_fused_diamond_df.loc[:, "overlap_score"] = ((filtered_fused_diamond_df["coverage_gene_1"] + filtered_fused_diamond_df["coverage_gene_2"]) / 2)                                                       
    
    # Load control results
    # control_file = f"{output_prefix}/{chrom}_{strand_name}_control_unique_diamond_results.tsv"
    # control_df = pd.read_csv(control_file, sep="\t")
    
    # For each fused gene, compare qcovhsp and scovhsp with the corresponding control genes
    # filtered_rows = []
    # for row in filtered_fused_diamond_df.itertuples():
        # gene1 = row.gene_1
        # gene2 = row.gene_2
        # qcov = row.query_coverage
        # scov = row.subject_coverage

        # gene1_control = control_df[control_df["gene"] == gene1]
        # gene2_control = control_df[control_df["gene"] == gene2]

        # if not gene1_control.empty and not gene2_control.empty:
            # gene1_qcov = gene1_control["query_coverage"].max()
            # gene2_qcov = gene2_control["query_coverage"].max()
            # gene1_scov = gene1_control["subject_coverage"].max()
            # gene2_scov = gene2_control["subject_coverage"].max()

            # Check if fused qcov/scov exceeds either gene1 or gene2
            # if (qcov > max(gene1_qcov, gene2_qcov)) or (scov > max(gene1_scov, gene2_scov)):
                # filtered_rows.append(row._asdict())

    # Convert filtered rows back to a DataFrame, preserving structure even if empty
    # if filtered_rows:
        # filtered_df = pd.DataFrame(filtered_rows)
    # else:
        # Create empty dataframe with same columns as filtered_fused_diamond_df
        # filtered_df = filtered_fused_diamond_df.iloc[0:0].copy()
    
    # Merge with overlap stuff
    # filtered_df = filtered_fused_diamond_df.merge(filtered_df, how="inner", on="fused_gene")
    
    # filtered_df.to_csv(f"{output_prefix}/{chrom}_{strand_name}_filtered_fused_diamond.csv")
    
    
    # Creating a dataframe that takes the greatest match for each alignment after filtering for easy viewing.
    # unique_filtered_fused_diamond_df = filtered_fused_diamond_df.groupby('fused_gene').first().reset_index()
    # unique_filtered_fused_diamond_df.to_csv(f"{output_prefix}/{chrom}_{strand_name}_unique_filtered_fused_diamond.csv")
    
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
    pattern = r'>(\S+)\s+(.+?)(?:\s+\[([^\]]+)\])?$'

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
                    description = match.group(2)
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
    print(f"Chromsomes/contigs that met the minimum gene count threshold of {min_geneperx_threshold}:\n{chr(10).join(sorted(chromosomes))}")
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
        valid_modes = ["--fast", "--mid-sensitive", "--sensitive", "--more-sensitive", "--very-sensitive", "--ultra-sensitive"]
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
                        help="Sensitivity mode parameter for the Diamond alignment tool. Valid options: --fast, --mid-sensitive, --sensitive, --more-sensitive, --very-sensitive, --ultra-sensitive",
                        type=valid_sensitivity,
                        default=None)

    args = parser.parse_args()



    annotation = args.organism_annotation
    organism_name = detect_organism_name(annotation)

    # TODO: Implement -n parameter
    if organism_name is None:
        print("ERROR: Organism name not detected. Please enter it manually using the '-n' or '--name' parameter.")
        exit()

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
    diamond_sensitivity = args.diamond_sensitivity

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
                    temp_faa_filepath, fused_metadata_df = chromosome_processor_fused_allisoforms(chrom, strand, headNode, output_prefix)

                    # Check if there are any fusions to process
                    if os.path.getsize(temp_faa_filepath) > 0 and not fused_metadata_df.empty:
                        pbar.set_description(f"Running DIAMOND on {chrom} {strand_name} strand's fused genes")
                        fused_chrom_hits_df, fused_diamond_df = fused_diamond_alignment(chrom, strand_name, diamond_path, db, temp_faa_filepath,
                                                                                fused_metadata_df, num_threads, ident_cutoff, output_prefix, diamond_sensitivity)
                        fused_hits_genome_df_list.append(fused_chrom_hits_df)
                        fused_diamond_df.to_csv(f"{output_prefix}/test_df.csv")
                    else:
                        pbar.set_description(f"Skipping {chrom} {strand_name} strand (no fusions to process)")
                        # Clean up empty temp file
                        if os.path.exists(temp_faa_filepath):
                            os.remove(temp_faa_filepath)
                    pbar.update(1)
                    
                    # Process unfused proteins second
                    temp_faa_filepath = chromosome_processor_unfused_allisoforms(chrom, strand, headNode, output_prefix, fused_chrom_hits_df)
                    if os.path.exists(temp_faa_filepath) and os.path.getsize(temp_faa_filepath) > 0:
                        with open(temp_faa_filepath, "r") as temp_faa, open(genome_faa_path, "a") as genome_faa:
                            genome_faa.write(temp_faa.read())
                    pbar.update(1)
                    
        #The control is to see if the alignment score increases in the fused gene vs the unfused genes for any overlapping hits in the fused genes
        pbar.set_description("Running control DIAMOND on identified gene parts")
        control_hits = unfused_diamond_alignment(chrom, strand_name, diamond_path, db, genome_faa_path, num_threads, output_prefix, organism_name, diamond_sensitivity)
        pbar.update(1)

    
    
    # Compiling all unique gene alignments to have one dataframe that contains genome-wide results.
    if fused_hits_genome_df_list:
        fused_hits = pd.concat(fused_hits_genome_df_list, ignore_index=True)
        max_eval_fuse = fused_hits.groupby("fused_protein")["expected_value"].idxmax()
        fused_hits_filter = fused_hits.loc[max_eval_fuse].reset_index(drop=True)
        max_eval_control = control_hits.groupby("protein")["expected_value"].idxmax()
        control_hits_filter = control_hits.loc[max_eval_control].reset_index(drop=True)
        # unique_genome_df = unique_genome_df[unique_genome_df["subject_organism"].str.lower() != name] # Filtering out hits to own organism this line problematic for some reason
        # unique_genome_df = unique_genome_df.sort_values(by="overlap_score", ascending=False)
        # unique_genome_df = unique_genome_df[(unique_genome_df["query_coverage"] > 70) &
                                            # (unique_genome_df["subject_coverage"] > 50)]
        fused_hits_filter.to_csv(f"{output_folder}/fused_hits.csv")
        control_hits_filter.to_csv(f"{output_folder}/control_hits.csv")
    else:
        print("No fused genes found across all chromosomes/strands.")
        # Create empty results file with proper headers
        pd.DataFrame(columns=["fused_gene", "fused_product", "fused_gene_len", "gene_1", "product_1",
                              "gene_1_len", "gene_2", "product_2", "gene_2_len"]).to_csv(f"{output_folder}/full_statistics_genome_results.csv")

    # Making it look pretty
    # final_output_df = unique_genome_df[["fused_gene", "fused_product", "subject_XP", "subject_product", "subject_organism", "fused_gene_len", "gene_1_len", "gene_2_len", "subject_length", "alignment_length", "coverage_gene_1", "coverage_gene_2", "overlap_score", "percentage_of_identical_matches"]]
    # final_output_df = unique_genome_df[["fused_gene", "fused_product", "subject_title", "fused_gene_len", "gene_1_len", "gene_2_len", "subject_length", "alignment_length", "coverage_gene_1", "coverage_gene_2", "overlap_score", "percentage_of_identical_matches"]]
    # final_output_df = final_output_df.sort_values(by="overlap_score", ascending=False)
    #final_output_df = unique_genome_df[["fused_gene", "fused_product", "subject_title", ""]]
    #final_output_df.to_csv(f"{output_folder}/final_genome_results.csv")
    print("Your results are ready.")