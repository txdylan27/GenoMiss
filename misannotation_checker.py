import pyranges as pr
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import time
import os
import subprocess

# reading in gff and genome
# need to make it so that these are all user input
start_time = time.time()
print("Loading in the gff and genome....")
os.chdir("/home/dylan/Files/Spyder/Intro_Prog_Project/apis_mellifera")
gff = pr.read_gff3("apis_mellifera_genomic.gff", as_df = True)
genome = SeqIO.index("GCF_003254395.2_Amel_HAv3.1_genomic.fna", "fasta")
cwd = os.getcwd()
end_time = time.time()
print(f"Finished loading. Time taken: {end_time - start_time:.2f} seconds")

# Extracting list of chromosomes from the gff
# NC are chromosomes, NW/NT are scaffolds
chrom_list = gff["Chromosome"].unique().tolist()
chrom_list.sort() # doesn't really need to be sorted
# chrom_list = [chrom for chrom in chrom_list if "NC" in chrom or "NT" in chrom] # temporary until we start testing other random genomes
# chrom_list = [chrom for chrom in chrom_list if "NW" in chrom] this varies, maybe take a chromosome annotation type input by the user then use whatever idek
chrom_list = [chrom for chrom in chrom_list if "NC" in chrom]
print(chrom_list)
# what are the NW chromosomes, it seems only the NC chromosomes are the actual chromosomes reported via RefSeq, maybe make it so that it outputs identified chromosomes to the user and asks if they are correct, if they are not can input the correct ones (or select via GUI)
feature_list = gff["Feature"].unique().tolist() # can also make it so you select which features to include
# need to make it so that they can select if they dont want the mitochondrial chromosome (this can be included in above feature)
# make this so it fits with the chrom list
gff_filtered = gff[
    gff["Chromosome"].str.contains("NC") |
    gff["Chromosome"].str.contains("NT") |
    gff["Chromosome"].str.contains("NW")
]# this needs to be based on the chromosome list
gff_sorted = gff_filtered.sort_values(by=["Start"])
gff_sorted_filtered = gff_sorted.query("Feature == 'exon' | Feature == 'gene' | Feature == 'CDS'") # for ease, only doing genes and exons right now
gff_sorted_filtered = gff_sorted_filtered[["Chromosome", "Feature", "gbkey", "gene", "product", "transcript_id", "Start", "End", "Strand"]]
gff_sorted_filtered = gff_sorted_filtered[~gff_sorted_filtered['gene'].str.contains('tRNA', case=False, na=False)]
gff_sorted_filtered = gff_sorted_filtered[~((gff_sorted_filtered["Feature"] == "exon") & (gff_sorted_filtered["product"].isna()))]
# put all filterings (except for chrom) up here, can honestly probably get rid of a lot of the filterings as this was for my own visual aid during the debugging

# can make this more compact by splitting it into functions and putting the DIAMOND stuff in the same chrom loop, can probably get rid of mRNA stuff as well since we are doing protein blast

# Going through gff by chromosome and strand, then going into the genome and grabbing the CDS and mRNA sequence for a gene. I then translate the CDS regions. I then combine genes that might have been fragmented due to genomic overlap.
# After, I sort each entry, split the proteins and mRNAs, choose the longest isoform/variant and output the dataframe for BLASTing.
# cycling through each chromosome from the chrom list
for chrom in chrom_list:
    # start timing how long it takes to go through the chromosome
    start_time = time.time()
    print("-----------------------------")
    print(f"Going through {chrom}...")
    print("-----------------------------")
    # getting the genome sequence correlating to the chromosome
    chrom_seq = genome[chrom].format("fasta")
    # getting rid of first line of chrom sequence, replacing all of the new line characters, and making the case uniform
    seq = chrom_seq.split("\n", 1)
    seq = seq[1].replace("\n", "")
    seq = seq.upper()
    # cycling through strands
    for strand in ["plus", "neg"]:
        chrom_strand_out = f"{cwd}/output/{chrom}/{strand}"
        if not os.path.exists(chrom_strand_out):
            os.makedirs(chrom_strand_out)
        # setting up strand-specific chromosome-specific list to store all entries then put that into a df after
        results_list = []
        # making sure that the gene list only consists of correct chromosome and strand
        gff_chrom = gff_sorted_filtered[gff_sorted_filtered["Chromosome"] == chrom]
        if strand == "plus":
            gff_chrom = gff_chrom[gff_chrom["Strand"] == "+"]
        elif strand == "neg":
            gff_chrom = gff_chrom[gff_chrom["Strand"] == "-"]
        
        exon_coords = []
        CDS_coords = []
        previous_gene = None
        exon_check = False
        CDS_check = False
        product_name = None
        variant_dict = {}
        isoform_dict = {}
        beginning = None
        ending = None

        # going through each gene entry
        for entry in gff_chrom.itertuples(index=False):
            # if this entry is an exon and is an exon of the previous gene, then add its coordinates to the coords list
            if entry.gene == previous_gene and entry.Feature == "exon":
                exon_check = True
                e_coord = [entry.Start, entry.End]
                # variant handling
                if "transcript variant" in str(entry.product).lower():
                    variant_name = entry.product
                    if variant_name in variant_dict:
                        variant_dict[variant_name].append(e_coord)
                    else:
                        variant_dict[variant_name] = [e_coord]
                else:
                    exon_coords.append(e_coord)
            # if this entry is a CDS and is a CDS of the previous gene, add its coordinates to the CDS list
            elif entry.gene == previous_gene and entry.Feature == "CDS":
                CDS_check = True
                p_coord = [entry.Start, entry.End]
                # isoform handling
                if "isoform" in str(entry.product).lower():
                    isoform_name = entry.product
                    if isoform_name in isoform_dict:
                        isoform_dict[isoform_name].append(p_coord)
                    else:
                        isoform_dict[isoform_name] = [p_coord]
                else:
                    CDS_coords.append(p_coord)
    
            # if this entry is the same gene but is not an exon, do not add its coordinates but continue going through the entries
            elif entry.gene == previous_gene and entry.Feature == "gene":
                continue
            
            # if this entry is not the same gene and not an exon, write out the gene sequence based on the coordinates
            else:
                if previous_gene is not None:
    
                    # check for variants
                    if variant_dict:
                        for variant, v_coords in variant_dict.items():
                            gene_dna_seq = ""
                            beginning = v_coords[0][0]
                            ending = v_coords[-1][-1]
                            for coords in v_coords:
                                start, end = coords
                                gene_dna_seq += seq[start:end]
                            gene_dna_seq = Seq(gene_dna_seq)
                            if strand == "neg":
                                gene_dna_seq = gene_dna_seq.reverse_complement()
                                results_list.append({'gene': previous_gene,'product_name': variant,'sequence_type': 'mRNA', "start": beginning, "end": ending, 'sequence': str(gene_dna_seq)})
                            else:
                                results_list.append({'gene': previous_gene,'product_name': variant,'sequence_type': 'mRNA', "start": beginning, "end": ending, 'sequence': str(gene_dna_seq)})
    
                    # check for isoforms
                    if isoform_dict:
                        for isoform, i_coords in isoform_dict.items():
                            gene_prot_seq = ""
                            beginning = i_coords[0][0]
                            ending = i_coords[-1][-1]
                            for coords in i_coords:
                                start, end = coords
                                gene_prot_seq += seq[start:end]
                            gene_prot_seq = Seq(gene_prot_seq)
                            if strand == "neg":
                                gene_prot_seq = gene_prot_seq.reverse_complement()
                                gene_prot_seq = gene_prot_seq.translate(to_stop=True)
                                results_list.append({'gene': previous_gene,'product_name': isoform,'sequence_type': 'protein', "start": beginning, "end": ending, 'sequence': str(gene_prot_seq)})
                            else:
                                gene_prot_seq = gene_prot_seq.translate(to_stop=True)
                                results_list.append({'gene': previous_gene,'product_name': isoform,'sequence_type': 'protein', "start": beginning, "end": ending, 'sequence': str(gene_prot_seq)})
                    
                    # no variants or isoforms
                    if CDS_coords:
                        # if the gene has CDS, translate
                        gene_prot_seq = ""
                        beginning = CDS_coords[0][0]
                        ending = CDS_coords[-1][-1]
                        for coords in CDS_coords:
                            start, end = coords
                            gene_prot_seq += seq[start:end]
                        gene_prot_seq = Seq(gene_prot_seq)
                        if strand == "neg":
                            gene_prot_seq = gene_prot_seq.reverse_complement()
                            gene_prot_seq = gene_prot_seq.translate(to_stop=True)
                            results_list.append({'gene': previous_gene,'product_name': product_name,'sequence_type': 'protein', "start": beginning, "end": ending, 'sequence': str(gene_prot_seq)})
                        else:
                            gene_prot_seq = gene_prot_seq.translate(to_stop=True)
                            results_list.append({'gene': previous_gene,'product_name': product_name,'sequence_type': 'protein', "start": beginning, "end": ending, 'sequence': str(gene_prot_seq)})
                            
                    # if no CDS, just output DNA seq
                    elif exon_coords:
                            gene_dna_seq = ""
                            beginning = exon_coords[0][0]
                            ending = exon_coords[-1][-1]
                            # go into genomic sequence and piece together the gene based on stored exon coordinates
                            for coords in exon_coords:
                                start, end = coords
                                gene_dna_seq += seq[start:end]
                            # convert the gene sequence into a Seq object so that we can easily manipulate it
                            gene_dna_seq = Seq(gene_dna_seq)
                            # if the gene is on the negative strand, get the reverse complement
                            if strand == "neg":
                                gene_dna_seq = gene_dna_seq.reverse_complement()
                                results_list.append({'gene': previous_gene,'product_name': product_name,'sequence_type': 'mRNA', "start": beginning, "end": ending, 'sequence': str(gene_dna_seq)})
                            else:
                                results_list.append({'gene': previous_gene,'product_name': product_name,'sequence_type': 'mRNA', "start": beginning, "end": ending, 'sequence': str(gene_dna_seq)})
                
                # reset the variables
                exon_coords = []
                CDS_coords = []
                exon_check = False
                CDS_check = False
                variant_dict = {}
                isoform_dict = {}
                beginning = None
                ending = None

                if entry.Feature == "exon":
                    exon_check = True
                    e_coord = [entry.Start, entry.End]
                    # variant handling
                    if "transcript variant" in str(entry.product).lower():
                        variant_name = entry.product
                        if variant_name in variant_dict:
                            variant_dict[variant_name].append(e_coord)
                        else:
                            variant_dict[variant_name] = [e_coord]
                    else:
                        exon_coords.append(e_coord)
                elif entry.Feature == "CDS":
                    CDS_check = True
                    p_coord = [entry.Start, entry.End]
                    # isoform handling
                    if "isoform" in str(entry.product).lower():
                        isoform_name = entry.product
                        if isoform_name in isoform_dict:
                            isoform_dict[isoform_name].append(p_coord)
                        else:
                            isoform_dict[isoform_name] = [p_coord]
                    else:
                        CDS_coords.append(p_coord)
    
            # store the gene and product names
            previous_gene = entry.gene
            product_name = entry.product
    
        # same as above, just doing for the last entry since it will kick us out of the loop
        if previous_gene is not None:
            # check for variants
            if variant_dict:
                for variant, v_coords in variant_dict.items():
                    gene_dna_seq = ""
                    beginning = v_coords[0][0]
                    ending = v_coords[-1][-1]
                    for coords in v_coords:
                        start, end = coords
                        gene_dna_seq += seq[start:end]
                    gene_dna_seq = Seq(gene_dna_seq)
                    if strand == "neg":
                        gene_dna_seq = gene_dna_seq.reverse_complement()
                        results_list.append({'gene': previous_gene,'product_name': variant,'sequence_type': 'mRNA', "start": beginning, "end": ending, 'sequence': str(gene_dna_seq)})
                    else:
                        results_list.append({'gene': previous_gene,'product_name': variant,'sequence_type': 'mRNA', "start": beginning, "end": ending, 'sequence': str(gene_dna_seq)})
    
            # check for isoforms
            if isoform_dict:
                for isoform, i_coords in isoform_dict.items():
                    gene_prot_seq = ""
                    beginning = i_coords[0][0]
                    ending = i_coords[-1][-1]
                    for coords in i_coords:
                        start, end = coords
                        gene_prot_seq += seq[start:end]
                    gene_prot_seq = Seq(gene_prot_seq)
                    if strand == "neg":
                        gene_prot_seq = gene_prot_seq.reverse_complement()
                        gene_prot_seq = gene_prot_seq.translate(to_stop=True)
                        results_list.append({'gene': previous_gene,'product_name': isoform,'sequence_type': 'protein', "start": beginning, "end": ending, 'sequence': str(gene_prot_seq)})
                    else:
                        gene_prot_seq = gene_prot_seq.translate(to_stop=True)
                        results_list.append({'gene': previous_gene,'product_name': isoform,'sequence_type': 'protein', "start": beginning, "end": ending, 'sequence': str(gene_prot_seq)})
            if CDS_coords:
                # if the gene has CDS, translate
                gene_prot_seq = ""
                beginning = CDS_coords[0][0]
                ending = CDS_coords[-1][-1]
                for coords in CDS_coords:
                    start, end = coords
                    gene_prot_seq += seq[start:end]
                gene_prot_seq = Seq(gene_prot_seq)
                if strand == "neg":
                    gene_prot_seq = gene_prot_seq.reverse_complement()
                    gene_prot_seq = gene_prot_seq.translate(to_stop=True)
                    results_list.append({'gene': previous_gene,'product_name': product_name,'sequence_type': 'protein', "start": beginning, "end": ending, 'sequence': str(gene_prot_seq)})
                else:
                    gene_prot_seq = gene_prot_seq.translate(to_stop=True)
                    results_list.append({'gene': previous_gene,'product_name': product_name,'sequence_type': 'protein', "start": beginning, "end": ending, 'sequence': str(gene_prot_seq)})
                    
            # if no CDS, just output DNA seq, but make sure it exists first
            elif exon_coords:
                gene_dna_seq = ""
                beginning = exon_coords[0][0]
                ending = exon_coords[-1][-1]
                # go into genomic sequence and piece together the gene based on stored exon coordinates
                for coords in exon_coords:
                    start, end = coords
                    gene_dna_seq += seq[start:end]
                # convert the gene sequence into a Seq object so that we can easily manipulate it
                gene_dna_seq = Seq(gene_dna_seq)

                # if the gene is on the negative strand, get the reverse complement
                if strand == "neg":
                    gene_dna_seq = gene_dna_seq.reverse_complement()
                    results_list.append({'gene': previous_gene,'product_name': product_name,'sequence_type': 'mRNA', "start": beginning, "end": ending, 'sequence': str(gene_dna_seq)})
                else:
                    results_list.append({'gene': previous_gene,'product_name': product_name,'sequence_type': 'mRNA', "start": beginning, "end": ending, 'sequence': str(gene_dna_seq)})
        

        # combine sequences that might have been fragmented from overlapping genomic regions
        results = pd.DataFrame(results_list)
        grouped = results.groupby('gene')
        combined_results = []
        
        # going through each gene group, separating isoforms, variants, and unique genes 
        for gene, group in grouped:
            # sorting group by start to maintain sequence order
            group = group.sort_values(by='start')
            
            # want to combine a gene if it was fragmented at the gene level, protein level, and transcript level
            base_gene = group[~group['product_name'].str.contains('isoform X|transcript variant X', case=False)]
            isoforms = group[group['product_name'].str.contains('isoform X', case=False)]
            variants = group[group['product_name'].str.contains('transcript variant X', case=False)]
            
            # combine sequences for base gene if fragmented
            combined_dict = {}
            for row in base_gene.itertuples(index=False):
                gene_name = row.gene
                if gene_name in combined_dict:
                    combined_dict[gene_name]['sequence'] += row.sequence
                    combined_dict[gene_name]['end'] = max(combined_dict[gene_name]['end'], row.end)
                else:
                    combined_dict[gene_name] = {'gene': gene_name, 'product_name': row.product_name, 'sequence_type': row.sequence_type, 'start': row.start, 'end': row.end, 'sequence': row.sequence}
            combined_results.extend(combined_dict.values())
            
            # combine sequences for isoforms if fragmented
            combined_dict = {}
            for row in isoforms.itertuples(index=False):
                product_name = row.product_name
                if product_name in combined_dict:
                    combined_dict[product_name]['sequence'] += row.sequence
                    combined_dict[product_name]['end'] = max(combined_dict[product_name]['end'], row.end)
                else:
                    combined_dict[product_name] = {'gene': row.gene, 'product_name': product_name, 'sequence_type': row.sequence_type, 'start': row.start, 'end': row.end, 'sequence': row.sequence}
            combined_results.extend(combined_dict.values())
            
            # combine sequences for variants if fragmented
            combined_dict = {}
            for row in variants.itertuples(index=False):
                product_name = row.product_name
                if product_name in combined_dict:
                    combined_dict[product_name]['sequence'] += row.sequence
                    combined_dict[product_name]['end'] = max(combined_dict[product_name]['end'], row.end)
                else:
                    combined_dict[product_name] = {'gene': row.gene, 'product_name': product_name, 'sequence_type': row.sequence_type, 'start': row.start, 'end': row.end, 'sequence': row.sequence}
            combined_results.extend(combined_dict.values())
        
        combined_results_df = pd.DataFrame(combined_results)
        combined_results_df = combined_results_df.sort_values(by=['start', 'end'])
        combined_results_df.to_csv(f"{chrom_strand_out}/{chrom}_{strand}_total_output.csv", index=False)
        combined_results_df = combined_results_df.drop_duplicates()

        # making separate protein and mRNA dataframes
        prot_df = combined_results_df[combined_results_df["sequence_type"] == "protein"]
        mRNA_df = combined_results_df[combined_results_df["sequence_type"] == "mRNA"]

        # keeping only the longest isoform or variant (intend to change this such that it only uses the longest variant/isoform if no Kallisto data is provided)
        prot_df['sequence_length'] = prot_df['sequence'].apply(len)
        mRNA_df['sequence_length'] = mRNA_df['sequence'].apply(len)
        
        prot_df = prot_df.loc[prot_df.groupby(['gene'])['sequence_length'].idxmax()]
        mRNA_df = mRNA_df.loc[mRNA_df.groupby(['gene'])['sequence_length'].idxmax()]
        
        prot_df = prot_df.drop(columns=['sequence_length'])
        mRNA_df = mRNA_df.drop(columns=['sequence_length'])
        
        prot_df = prot_df.sort_values(by=['start', 'end'])
        mRNA_df = mRNA_df.sort_values(by=['start', 'end'])
        
        prot_df.to_csv(f"{chrom_strand_out}/{chrom}_{strand}_prot_output.csv", index=False)
        mRNA_df.to_csv(f"{chrom_strand_out}/{chrom}_{strand}_mRNA_output.csv", index=False)
        results_list.clear()
        combined_results.clear()

        
    # end timing
    end_time = time.time()
    print("-----------------------------")
    print(f"Finished processing {chrom}. Time taken: {end_time - start_time:.2f} seconds")
    print("-----------------------------")


# running fusion diamond (blast) with upstream genes
# make these user-inputted
diamond_path = "/media/dylan/maindata_2/protein_tools/diamond"
db_path = "/media/dylan/maindata_2/protein_tools/insecta_refseq_protein_db.dmnd"
num_threads = "12"
ident_cutoff = 94
# length_tolerance = 0.1
len_buffer = 100
unique_genome_df_list = []

for chrom in chrom_list:
    # start timing how long it takes to BLAST a single chromosome
    start_time = time.time()
    print("-----------------------------")
    print(f"Running DIAMOND on Fusion Genes in {chrom}...")
    print("-----------------------------")
    for strand in ["plus", "neg"]:
        chrom_strand_prot_file = f"{cwd}/output/{chrom}/{strand}/{chrom}_{strand}_prot_output.csv"
        prot_df = pd.read_csv(chrom_strand_prot_file)
        temp_fasta_path = "temp_fasta.faa"
        fused_prot_list = []
        with open(temp_fasta_path, "w") as temp_fasta:
            idx = 0
            for protein in prot_df.itertuples(index=False):
                if idx == 0:
                    previous_prot = protein
                    idx += 1
                    continue
                else:
                    current_prot = protein
                    fused_prot_gene = f"{previous_prot.gene}+{current_prot.gene}"
                    fused_prot_product = f"{previous_prot.product_name}+{current_prot.product_name}"
                    fused_prot_seq = previous_prot.sequence + current_prot.sequence
                    fused_prot_list.append({"fused_gene": fused_prot_gene, "fused_product": fused_prot_product, "fused_gene_len": len(fused_prot_seq),
                                            "gene_1": previous_prot.gene, "product_1": previous_prot.product_name, "gene_1_len": len(previous_prot.sequence),
                                            "gene_2": current_prot.gene, "product_2": current_prot.product_name, "gene_2_len": len(current_prot.sequence)})
                    temp_fasta.write(f">{fused_prot_gene} {fused_prot_product}\n{fused_prot_seq}\n")
                    previous_prot = protein
                    idx += 1
        fused_prot_df = pd.DataFrame(fused_prot_list)
        diamond_output = f"{cwd}/output/{chrom}/{strand}/{chrom}_{strand}_diamond_results.tsv"
        diamond_command = [diamond_path, "blastp", "--db", db_path, "--query", temp_fasta_path, "--out", diamond_output, 
                           "--outfmt", "6", "qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "sstart", "send",
                           "pident", "nident", "mismatch", "evalue", "bitscore", "length", "qtitle", "stitle",
                           "--header", "--evalue", "1e-5", "--threads", num_threads]
        try:
            diamond_result = subprocess.run(diamond_command, check=True, text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print(f"Error during DIAMOND execution for {chrom}, {strand} strand.")
        headers = ["fused_gene", "query_length", "subject_id", "subject_length", "start_of_alignment_in_query", "end_of_alignment_in_query", 
                   "start_of_alignment_in_subject", "end_of_alignment_in_subject", "percentage_of_identical_matches", "number_of_identical_matches", 
                   "number_of_mismatches", "expected_value", "bit_score", "alignment_length", "query_title", "subject_title"]
        diamond_df = pd.read_csv(f"{cwd}/output/{chrom}/{strand}/{chrom}_{strand}_diamond_results.tsv", sep = "\t", skiprows = 3, names = headers)
        # fused_diamond_df = pd.merge(fused_prot_df, diamond_df, on = "fused_gene", how = "outer")
        fused_diamond_df = fused_prot_df.merge(fused_prot_df.merge(diamond_df, how='outer', on='fused_gene', sort=False)) # weird gimmicky solution I found to original order not being maintained after merging
        fused_diamond_df.to_csv(f"{cwd}/output/{chrom}/{strand}/{chrom}_{strand}_fused_diamond.csv")
        # get rid of exact matches to either gene part of the fused gene
        # filtered_diamond_df = fused_diamond_df[((fused_diamond_df["subject_length"] != fused_diamond_df["gene_1_len"]) &
        # (fused_diamond_df["subject_length"] != fused_diamond_df["gene_2_len"]))]
        # check that the remaining entries meet the length tolerance (falls apart with longer sequences, hmmmm) maybe do a hard cutoff - this was for percentage cutoff
        # filtered_diamond_df["length_difference"] = abs(filtered_diamond_df["query_length"] - filtered_diamond_df["subject_length"])
        # filtered_diamond_df["length_within_tolerance"] = filtered_diamond_df["length_difference"] / filtered_diamond_df["query_length"] <= length_tolerance
        # filtered_diamond_df = filtered_diamond_df[filtered_diamond_df["length_within_tolerance"]]
        # this is hard length cutoff, added alignment length cutoff as well, could maybe do a dynamic cutoff that changes based on length of gene size (> 50% of gene 1 len or 50% of gene 2 len)
        filtered_fused_diamond_df = fused_diamond_df[((fused_diamond_df["subject_length"] > (fused_diamond_df["gene_1_len"] + len_buffer)) | (fused_diamond_df["subject_length"] > (fused_diamond_df["gene_2_len"] + len_buffer)))
        & (fused_diamond_df["percentage_of_identical_matches"] > ident_cutoff) & ((fused_diamond_df["alignment_length"] > (fused_diamond_df["gene_1_len"] + len_buffer)) & (fused_diamond_df["alignment_length"] > (fused_diamond_df["gene_2_len"] + len_buffer)))] 
        filtered_fused_diamond_df.to_csv(f"{cwd}/output/{chrom}/{strand}/{chrom}_{strand}_full_filtered_fused_diamond.csv")
        # choosing unique based on highest bit score
        unique_filtered_fused_diamond_df = filtered_fused_diamond_df.groupby('fused_gene').first().reset_index()
        unique_filtered_fused_diamond_df.to_csv(f"{cwd}/output/{chrom}/{strand}/{chrom}_{strand}_unique_filtered_fused_diamond.csv")
        unique_genome_df_list.append(unique_filtered_fused_diamond_df)
        
    end_time = time.time()
    print(f"DIAMOND for Fusion Genes in {chrom} was completed in {end_time - start_time:.2f} seconds")

unique_genome_df = pd.concat(unique_genome_df_list, ignore_index=True)
unique_genome_df.to_csv(f"{cwd}/output/unique_genome_fused_results.csv")

# need to tighten parameters, make it so that there is genome-wide df as well that appends and prints after, check alignment quality in addition to length cutoff

