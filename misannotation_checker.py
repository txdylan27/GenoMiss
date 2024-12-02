import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import os
import pyranges as pr
import pandas as pd
import shutil
import subprocess
import sys
import time
import tkinter as tk
from tkinter import ttk
from tqdm import tqdm

# Defining all functions for the program
def get_default_num_threads():
    """Function used to identify the number of threads on the user's system and then use half of them for DIAMOND if no threads are inputted."""
    
    total_threads = os.cpu_count() # Attempting to get the thread count for the user's system.
    if total_threads is None: 
        total_threads = 4 # Included in case - for some reason - Python is not able to identify their thread count.
        
    return max(1, total_threads // 2) # Returning the max of either 1 thread or half of the total_threads. The 1 is a fall-back in the rare instance of a system having 0 or only 1 thread.
    
def genome_loader(genome, annotation):
    """Function used to load in the provided genome and gff of the organism. Returns an indexed version of the genome and the gff as a dataframe."""
    
    print("Loading in genome and annotation...")
    
    idx_genome = SeqIO.index(genome, "fasta") # This is allowing for the indexing of the genome by chromosome without the loading of the entire genome.
    gff = pr.read_gff3(annotation, as_df=True) # Reorganizing the gff into a dataframe for easy manipulation.
    
    return idx_genome, gff

def chromosome_parser(gff):
    """Function used to identify the chromosomes within the genome and gff based on the user's input. This can be aided by using NCBI's entry of the user's
    genome. Returns a filtered chromosome list."""
    
    unfiltered_chrom_list = gff["Chromosome"].unique().tolist() # Getting a list of all of the unique chromosome identifiers present in the annotation file.
    while True:
        try:
            chrom_input = input("Input the chromsome prefixes - comma-separated - listed by RefSeq on the NCBI entry for your genome: ")
            if not chrom_input:
                raise ValueError("Input cannot be empty. Please enter at least one chromosome prefix.")

            prefix_list = [prefix.strip().upper() for prefix in chrom_input.split(",")] # Create a list of the prefixes provided by the user.
            
            chrom_list = [chrom for chrom in unfiltered_chrom_list if any(prefix in chrom for prefix in prefix_list)] # Refine the initial chromosome list based on the user's inputted prefixes.
            if not chrom_list:
                print("No chromosomes were found within the organism's annotation with the provided prefixes. Please try again.")
                continue
            
            return chrom_list
        
        # Error handling
        except ValueError as ve:
            print(f"Input error: {ve}")
        except Exception as e:
            print(f"Unexpected error: {e}")

def chromosome_selector(chrom_list):
    """Function that uses a GUI tool to allow the user to manually check the identified chromosomes before translation and fusion of the genome. Returns
    a manually selected and filtered chromosome list."""
    
    # Creating a list that will store the selected chromosomes by the checkboxes.
    selected_chrom = []
    
    def on_submit():
        """Function that extracts the selected chromosomes and then closes the GUI."""
        selected_chrom.extend([chrom for chrom, var in checkboxes.items() if var.get()]) # Adds the chromosome to the list if it is checked once the submit button is pressed.
        root.destroy() # Close the GUI window
    
    # Initialize GUI window
    root = tk.Tk()
    root.title("Chromosome Selection")
    ttk.Label(root, text="""This is a list of chromosomes identified in the annotation file based on the inputted prefixes. Check which chromosomes you
              would like to keep or remove, such as the mitochondrial chromsome, and then click submit when finished.""").pack(pady=10)
    checkboxes = {}
    for chrom in chrom_list:
        var = tk.BooleanVar(value=True) # Set each chromosome to default as checked.
        checkbox = ttk.Checkbutton(root, text=chrom, variable=var)
        checkbox.pack(anchor="w", padx=10, pady=2)
        checkboxes[chrom] = var
    submit_button = ttk.Button(root, text="Submit", command=on_submit)
    submit_button.pack(pady=10)
    root.mainloop()
    
    return selected_chrom

def gff_manipulator(gff, chrom_list):
    """Function containing all manipulations (filterings, sortings, etc.) of the annotation dataframe. Returns a filtered annotation dataframe."""
    
    filtered_gff = gff[gff["Chromosome"].isin(chrom_list)] # Creating a df that only contains annotations belonging to genes within the selected chromosomes.
    filtered_sorted_gff = filtered_gff.sort_values(by=["Start"]) # Making sure that the genes are sorted by position.
    filtered_sorted_gff = filtered_sorted_gff[filtered_sorted_gff["Feature"] == "CDS"] # Since we are doing protein alignment, only need the coding strands of each gene.
    
    return filtered_sorted_gff

def chromosome_processor(chrom, strand, strand_name, genome, gff, flag):
    """Function that compiles multiple functions to return a dataframe containing all translated proteins for each gene on each chromosome for each strand."""
    
    chrom_seq = genome[chrom].format("fasta") # Going into the indexed genome and retrieving the value (chromosome sequence) relating to the key (chromosome name).
    
    # Manipulating the chromosome sequence by discarding the sequence title line, replacing new line characters, and making the character case uniform.
    chrom_seq = chrom_seq.split("\n", 1)
    chrom_seq = chrom_seq[1].replace("\n", "")
    chrom_seq = chrom_seq.upper()
    
    # Manipulating gff dataframe to only have entries related to the current chromosome and strand.
    c_gff = gff[gff["Chromosome"] == chrom]
    if strand == "+":
        cs_gff = c_gff[c_gff["Strand"] == "+"]
    elif strand == "-":
        cs_gff = c_gff[c_gff["Strand"] == "-"]
    
    processed_gene = gene_processor(chrom_seq, cs_gff, strand) # Using CDS coordinates to build gene sequences while maintaining positional order.
    processed_gene_df = pd.DataFrame(processed_gene) # Turning the processed_gene list into a dataframe.
    total_prot_df = gene_joiner(processed_gene_df) # Combining DNA sequences that may have been fragmented.
    refined_prot_df = isoform_picker(total_prot_df) # Keeping the longest isoforms and discarding the rest.
    translated_df = sequence_translator(refined_prot_df) # Translating the kept sequences.
      
    translated_df.to_csv(f"{output_prefix}/{chrom}_{strand_name}_gene2protein.csv", index=False)
    
    return translated_df

def gene_processor(chrom_seq, cs_gff, strand):
    """Function that uses the CDS coordinates within the gff to retrieve all gene coding sequences on a specific strand of a chromosome. Returns
    a list of dictionaries containing information about each gene and their full CDS sequence based on the coordinates."""
    
    # Establishing variables.
    results = [] # Using a list to store dictionaries that contain each gene's protein information for faster processing.
    CDS_coords = [] # Using a list to store each coordinate pair for a gene's coding strand exons.
    previous_product = None
    previous_gene = None
    isoform_info = {} # Using a dictionary to store the names and CDS coordinates of all isoforms of a gene if they exist.
    
    # Going through each gene entry and processing them.
    for entry in cs_gff.itertuples(index=False): # Using itertuples as a memory-efficient and fast way to access each row of the dataframe.   
            
        # Evaluating if the current gene entry is a continuation of the previous gene's coding strand exons.
        if entry.gene == previous_gene:
            exon_coords = [entry.Start, entry.End] # If it is, record the exon's genomic coordinates.
            if "isoform" in str(entry.product).lower(): # Check for isoforms. If they exist, record their information.
                isoform_name = entry.product
                if isoform_name in isoform_info: # If this is the first instance of this isoform, add it to the isoform storage dictionary.
                    isoform_info[isoform_name].append(exon_coords)
                else: # If this is another exon of the isoform, just add the exon coordinates to the previous entry.
                    isoform_info[isoform_name] = [exon_coords]
            else: # If there is no isoform found, just add the genomic coordinates to the CDS dictionary.
                CDS_coords.append(exon_coords)
        
        # If the current gene entry is not a continuation of the previous gene entry, write out the gene's entire CDS sequence. The not statement is to exclude the first gene entry.
        elif previous_gene is not None and entry.gene != previous_gene:
            if isoform_info: # If the isoform dictionary exists, then write out the contents of it.
                for isoform, i_coords in isoform_info.items():
                    gene_seq = ""
                    beginning = i_coords[0][0] # Setting a variable for the beginning of the isoform, dictated by the first coordinate of the first coordinate pair.
                    ending = i_coords[-1][-1] # Setting a variable for the ending of the isoform, dictated by the last coordinate of the last coordinate pair.
                    # Extracting the genomic sequence of the isoform using the stored coordinates.
                    for coords in i_coords:
                        start, end = coords
                        gene_seq += chrom_seq[start:end]
                    gene_seq = Seq(gene_seq) # Making the extracted gene sequence a Seq object for easy manipulation (reverse complement).
                    # If the negative strand is being evaluated, take the reverse complement.
                    if strand == "-":
                        gene_seq = gene_seq.reverse_complement()
                    results.append({'gene': previous_gene,'product_name': isoform, "start": beginning, "end": ending, 'sequence': str(gene_seq)}) # Append the extracted isoform information to the results list.
            if CDS_coords: # If the CDS dictionary exists (which it should whenever there are no isoforms), then write out the contents of it. Logic is the same as for the isoform extraction.
                gene_seq = ""
                beginning = CDS_coords[0][0] 
                ending = CDS_coords[-1][-1] 
                for coords in CDS_coords:
                    start, end = coords
                    gene_seq += chrom_seq[start:end]
                gene_seq = Seq(gene_seq)
                if strand == "-":
                    gene_seq = gene_seq.reverse_complement()
                results.append({'gene': previous_gene,'product_name': previous_product, "start": beginning, "end": ending, 'sequence': str(gene_seq)})
        
            # Resetting variables since the current entry is a new gene.
            CDS_coords = []
            isoform_info = {}
        
            # Recording information of the current gene entry since it is the start of a new gene.
            exon_coords = [entry.Start, entry.End]
            if "isoform" in str(entry.product).lower():
                isoform_name = entry.product
                if isoform_name in isoform_info:
                    isoform_info[isoform_name].append(exon_coords)
                else:
                    isoform_info[isoform_name] = [exon_coords]
            else:
                CDS_coords.append(exon_coords)
            
        # Establishing variables for loops.
        previous_gene = entry.gene
        previous_product = entry.product

    # Loop will break when it goes past the last entry, need to account for the last gene and add it's information to the results.
    if isoform_info:
        for isoform, i_coords in isoform_info.items():
            gene_seq = ""
            beginning = i_coords[0][0]
            ending = i_coords[-1][-1]
            for coords in i_coords:
                start, end = coords
                gene_seq += chrom_seq[start:end]
            gene_seq = Seq(gene_seq)
            if strand == "-":
                gene_seq = gene_seq.reverse_complement()
            results.append({'gene': previous_gene,'product_name': isoform, "start": beginning, "end": ending, 'sequence': str(gene_seq)})
    if CDS_coords: 
        gene_seq = ""
        beginning = CDS_coords[0][0] 
        ending = CDS_coords[-1][-1] 
        for coords in CDS_coords:
            start, end = coords
            gene_seq += chrom_seq[start:end]
        gene_seq = Seq(gene_seq)
        if strand == "-":
            gene_seq = gene_seq.reverse_complement()
        results.append({'gene': previous_gene,'product_name': previous_product, "start": beginning, "end": ending, 'sequence': str(gene_seq)})
        
    return results

def gene_joiner(processed_df):
    """Function that combines gene sequences that may have been fragmented due to overlapping gene regions leading to breaks in gene entries. Returns
    a dataframe containing full-length gene sequences."""
    
    grouped_dict = processed_df.groupby("gene") # This creates a dictionary where the keys are genes and the values are their indices within the dataframe.
    joined_results = [] # Using a list to store the results to build the dataframe after.
    
    for gene, indices in grouped_dict:
        indices = indices.sort_values(by='start') # Sorting gene entries by start to maintain sequence order.
        
        # Checking to see if the gene has isoforms or not.
        base_genes = indices[~indices['product_name'].astype(str).str.contains('isoform X', case=False, na=False)]
        isoforms = indices[indices['product_name'].astype(str).str.contains('isoform X', case=False, na=False)]
        
        # For genes with no isoforms, see if there are multiple entries for the same gene name. If there are, combine their sequences.
        combined_base_genes = {}
        for row in base_genes.itertuples(index=False):
            gene_name = row.gene
            if gene_name in combined_base_genes:
                combined_base_genes[gene_name]['sequence'] += row.sequence
                combined_base_genes[gene_name]['end'] = max(combined_base_genes[gene_name]['end'], row.end)
            else:
                combined_base_genes[gene_name] = {'gene': gene_name, 'product_name': row.product_name, 'start': row.start, 'end': row.end, 'sequence': row.sequence}
        joined_results.extend(combined_base_genes.values())
        
        # For genes with isoforms, see if there are multiple entries for the same product name. If there are, combine their sequences.
        combined_isoforms = {}
        for row in isoforms.itertuples(index=False):
            product_name = row.product_name
            if product_name in combined_isoforms:
                combined_isoforms[product_name]['sequence'] += row.sequence
                combined_isoforms[product_name]['end'] = max(combined_isoforms[product_name]['end'], row.end)
            else:
                combined_isoforms[product_name] = {'gene': row.gene, 'product_name': product_name, 'start': row.start, 'end': row.end, 'sequence': row.sequence}
        joined_results.extend(combined_isoforms.values())
        
    joined_results = pd.DataFrame(joined_results)
    joined_results = joined_results.sort_values(by=['start', 'end'])
    
    return joined_results
        
def isoform_picker(joined_df):
    """Function that goes through the inputted dataframe and analyzes the length of each gene, keeping only the longest gene of multiple 
    entries of the same gene. This is a means to keep only the longest isoform of a gene with multiple isoforms. Returns a dataframe with 
    only one entry per gene."""
    
    joined_df['sequence_length'] = joined_df['sequence'].apply(len) # Get the length of each gene's sequence.
    joined_df = joined_df.loc[joined_df.groupby(['gene'])['sequence_length'].idxmax()] # Keep only the gene with the longest sequence length.
    joined_df = joined_df.drop(columns=['sequence_length'])
    joined_df = joined_df.sort_values(by=['start', 'end'])
    
    return joined_df

def sequence_translator(refined_df):
    """Function that translates the kept gene CDS sequences. To take into account rare, incomplete gene sequences due to the annotation process,
    padding is added when applicable to keep triplets. Returns a dataframe containing the translated sequences."""
    
    refined_df['protein_sequence'] = refined_df['sequence'].apply(lambda seq: str(Seq(seq + "N" * ((3 - len(seq) % 3) % 3)).translate(to_stop=True)))
    translated_df = refined_df.drop(columns=['sequence'])
    translated_df = translated_df.rename(columns={'protein_sequence': 'sequence'})

    return translated_df
    
def diamond_alignment(prot_output_df, chrom, strand, strand_name, cwd, diamond_path, db_path, num_threads, ident_cutoff, len_buffer, output_prefix):
    """Function used to fuse neigboring genes and run the DIAMOND protein alignment script on the fused genes. Returns a dataframe containing unique genes that
    fit the user-inputted filtering criteria."""   
    
    temp_fasta_path = f"{output_prefix}/temp_fasta.faa" # The input for DIAMOND is a fasta file. Therefore, need to convert the protein dataframe into a fasta format.
    fused_prot_df_list = [] # This list will be used to store information about the fused proteins.
    
    with open(temp_fasta_path, "w") as temp_fasta:
        first = True
        for prot in prot_output_df.itertuples(index=False): # Going into each row (protein) of the protein dataframe.
            if first == True: # The program works by fusing proteins upstream. So the current protein in the loop will be fused with the protein upstream of it. Since there is nothing upstream of the first protein on the chromosome, it is skipped.
                previous_prot = prot # Keeping track of the previous protein.
                first = False
                continue
            else:
                current_prot = prot # Keeping track of the current protein.
                fused_prot_gene = f"{previous_prot.gene}+{current_prot.gene}" # Fusing gene names.
                fused_prot_product = f"{previous_prot.product_name}+{current_prot.product_name}" # Fusing product names.
                fused_prot_seq = previous_prot.sequence + current_prot.sequence # Fusing the protein sequences.
                fused_prot_df_list.append({"fused_gene": fused_prot_gene, "fused_product": fused_prot_product, "fused_gene_len": len(fused_prot_seq),
                                        "gene_1": previous_prot.gene, "product_1": previous_prot.product_name, "gene_1_len": len(previous_prot.sequence),
                                        "gene_2": current_prot.gene, "product_2": current_prot.product_name, "gene_2_len": len(current_prot.sequence)}) # Adding a dictionary of the information to the list.
                temp_fasta.write(f">{fused_prot_gene} {fused_prot_product}\n{fused_prot_seq}\n") # Writing to the fasta.
                previous_prot = prot # Keeping track of the previous protein as we move to the next protein in the dataframe.
                
    fused_prot_df = pd.DataFrame(fused_prot_df_list)
    
    # Creating the diamond output file path and the diamond command that will be ran in the subprocess function. 
    diamond_output = f"{output_prefix}/{chrom}_{strand_name}_diamond_results.tsv"
    diamond_command = [diamond_path, "blastp", "--db", db_path, "--query", temp_fasta_path, "--out", diamond_output, 
                       "--outfmt", "6", "qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "sstart", "send",
                       "pident", "nident", "mismatch", "evalue", "bitscore", "length", "qtitle", "stitle",
                       "--header", "--evalue", "1e-5", "--threads", num_threads]
    
    # Running the diamond command via subprocess.
    try:
        subprocess.run(diamond_command, check=True, text=True, capture_output=True)
    except subprocess.CalledProcessError:
        print(f"Error during DIAMOND execution for {chrom}, {strand_name} strand.")
        raise
    
    # As the diamond headers are very obfuscated, included proper header titles.
    headers = ["fused_gene", "query_length", "subject_id", "subject_length", "start_of_alignment_in_query", "end_of_alignment_in_query", 
               "start_of_alignment_in_subject", "end_of_alignment_in_subject", "percentage_of_identical_matches", "number_of_identical_matches", 
               "number_of_mismatches", "expected_value", "bit_score", "alignment_length", "query_title", "subject_title"]
    
    # Converting the .tsv provided by diamond to a dataframe while skipping useless rows (the first 3) and renaming the headers.
    diamond_df = pd.read_csv(f"{output_prefix}/{chrom}_{strand_name}_diamond_results.tsv", sep="\t", skiprows = 3, names = headers)
    
    # Merging the diamond results dataframe with the fused protein dataframe. Had to use a weird gimmicky solution I found to keep the original order preserved.
    fused_diamond_df = fused_prot_df.merge(fused_prot_df.merge(diamond_df, how='outer', on='fused_gene', sort=False))
    fused_diamond_df.to_csv(f"{output_prefix}/{chrom}_{strand_name}_diamond_results.tsv", sep="\t")
    
    # Applying the filters inputted by the user.
    filtered_fused_diamond_df = fused_diamond_df[((fused_diamond_df["subject_length"] > (fused_diamond_df["gene_1_len"] + len_buffer)) | (fused_diamond_df["subject_length"] > (fused_diamond_df["gene_2_len"] + len_buffer)))
    & (fused_diamond_df["percentage_of_identical_matches"] > ident_cutoff) & ((fused_diamond_df["alignment_length"] > (fused_diamond_df["gene_1_len"] + len_buffer)) & (fused_diamond_df["alignment_length"] > (fused_diamond_df["gene_2_len"] + len_buffer)))] 
    filtered_fused_diamond_df.to_csv(f"{output_prefix}/{chrom}_{strand_name}_full_filtered_fused_diamond.csv")
    
    # Creating a dataframe that takes the highest match for each alignment after filtering for easy viewing.
    unique_filtered_fused_diamond_df = filtered_fused_diamond_df.groupby('fused_gene').first().reset_index()
    unique_filtered_fused_diamond_df.to_csv(f"{output_prefix}/{chrom}_{strand_name}_unique_filtered_fused_diamond.csv")
    
    os.remove(temp_fasta_path)
    
    return unique_filtered_fused_diamond_df
    

# Establishing command-line arguments
parser = argparse.ArgumentParser(
                    prog="Genome Misannotation Checker",
                    description="""This program takes an organism's annotation and genome files, extracts the amino acid sequences of all 
                    proteins in a sequential manner from the genome by chromosome and strand, then fuses neighboring gene sequences. These fused 
                    sequences are then aligned with a local reference protein database using the command-line protein alignment tool DIAMOND. 
                    Based on user-inputted criteria, potential misannotated genes are returned to the user in a CSV format. This program only works
                    with genomes containing assembled chromosomes. This program is not compatible with scaffolded genomes.""")
parser.add_argument('-wd', '--project_directory', help="""Input the file path to the project directory. Ideally this directory will contain the organism's 
                    genome, gff, as well as the protein database you wish to use with DIAMOND. If an input is not provided, the current directory 
                    will be used.""", default = os.getcwd())
parser.add_argument('-g', '--organism_genome', help="""Input the organism's genome in a fasta format. This should be a RefSeq genome. This input is required.""", required=True)
parser.add_argument('-a', '--organism_annotation', help="""Input the organism's annotation features in a gff format. This should be a RefSeq annotation. This input is required.""", required=True)
parser.add_argument('-db', '--database_path', help="""Input the file path to where the local reference protein database is stored. If an input is 
                    not provided, the directory provided in -wd will be used.""")
parser.add_argument('-t', '--num_threads', help="""Input the number of threads that you would like to use. By default, half of your available threads
                    will be used.""", default = get_default_num_threads())
parser.add_argument('-i', '--identity_cutoff', help="""Input the percent identity cutoff you would like to use for filtering of the fused gene alignments.
                    The percent identity is the percentage of identical amino acids between two sequences at the same alignment positions. The default is 0.75.""",
                    default = 0.75)
parser.add_argument('-l', '--length_buffer', help="""Input the length buffer you would like to use for filtering of the fused gene alignments. The length
                    buffer refers to the value added to the subject length and alignment length as a means to avoid protein alignments to the indivdual gene parts
                    of the fused gene. The default is 50 amino acids.""", default = 100)

args = parser.parse_args()

# Establishing variables from args.
cwd = args.project_directory
if cwd.endswith("/"):
    cwd = cwd[:-1]
genome = args.organism_genome
annotation = args.organism_annotation
diamond_path = shutil.which("diamond")
if args.database_path is None:
    db_path = cwd
else:
    db_path = args.database_path
num_threads = str(args.num_threads)
ident_cutoff = float(args.identity_cutoff) * 100
len_buffer = float(args.length_buffer)

# Storing the command-line arguments inputted by the user in a log file for future reference by the user.
if not os.path.exists("output"):
    os.makedirs("output")
with open("output/args.log", 'w') as log_file:
    log_file.write("Command-line input:\n")
    log_file.write(" ".join(sys.argv) + "\n")

# Loading in the organism's genome and annotation.
os.chdir(cwd)
start_time = time.time()
idx_genome, gff = genome_loader(genome, annotation)
end_time = time.time()
print(f"It took {end_time - start_time:0.2f} seconds to load the provided genome.")

# Determining the chromosome format within the genome, returning a list of them which has been curated by the user, and then filtering the gff df based on it.
filtered_chrom_list = chromosome_parser(gff)
selected_chrom_list = chromosome_selector(filtered_chrom_list)

# Manipulating the df for iteration
filtered_gff = gff_manipulator(gff, selected_chrom_list)

# A list to store each chromosome strand's top gene alignments.
fused_hits_genome_df_list = []

# Initializing a progress bar for tracking the program's status.
total_steps = len(selected_chrom_list) * 4 # For each chromosome, there are two strands. For each strand, there are two steps (processing and protein alignment).
with tqdm(total=total_steps, unit="step") as pbar:
    
    # Going into each strand of each chromosome, getting each protein's sequence, fusing neighboring genes, then using DIAMOND to check for misannotations.
    for chrom in selected_chrom_list:
        for strand in ["+", "-"]:
            # Avoiding use of + in file names since it is a no no in the rubric :)
            if strand == "+":
                strand_name = "plus"
            elif strand == "-":
                strand_name = "neg"
            
            output_prefix = f"output/{chrom}/{strand_name}"  # Setting up the output file path.
            if not os.path.exists(output_prefix):
                os.makedirs(output_prefix)
            
            pbar.set_description(f"Processing {chrom} {strand_name} strand")
            prot_output_df = chromosome_processor(chrom, strand, strand_name, idx_genome, filtered_gff, output_prefix)
            pbar.update(1)
            
            pbar.set_description(f"Running DIAMOND on {chrom} {strand_name} strand's fused genes")
            unique_chrom_strand_df = diamond_alignment(prot_output_df, chrom, strand, strand_name, cwd, diamond_path, db_path, 
                                                      num_threads, ident_cutoff, len_buffer, output_prefix)
            pbar.update(1)
            fused_hits_genome_df_list.append(unique_chrom_strand_df)

# Compiling all unique gene alignments to have one dataframe that contains genome-wide results.
unique_genome_df = pd.concat(fused_hits_genome_df_list, ignore_index=True)
unique_genome_df.to_csv("output/unique_genome_fused_results.csv")
