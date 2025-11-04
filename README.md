# GenoMiss

![Python](https://img.shields.io/badge/python-3.12.7-blue.svg)
![GitHub issues](https://img.shields.io/github/issues/txdylan27/GenoMiss)
![GitHub last commit](https://img.shields.io/github/last-commit/txdylan27/GenoMiss)
![License](https://img.shields.io/github/license/txdylan27/GenoMiss)

GenoMiss is a computational tool for detecting potential gene misannotations in protein-coding genomes. GenoMiss identifies cases where adjacent genes may actually represent fragments of a single misannotated gene by fusing neighboring genes and comparing their alignment scores against a reference protein database.

## Overview

Gene misannotations can occur when sequencing errors, assembly issues, or annotation pipeline limitations cause a single gene to be incorrectly split into multiple adjacent genes. This tool addresses this problem by:

1. **Constructing a genome graph** representing the spatial relationships between protein-coding genes
2. **Creating fusion proteins** by concatenating adjacent gene sequences
3. **Aligning fused and unfused proteins** against a reference database using DIAMOND
4. **Scoring alignment quality** using a composite scoring system that evaluates multiple factors
5. **Generating comprehensive reports** in multiple formats (CSV, TSV, Excel)

## Features

- **Graph-based genome representation**: Handles overlapping genes and complex genomic architectures using strand-specific directed graphs
- **Comprehensive scoring system**: Evaluates fused genes based on query coverage, bit score improvement, alignment equilibrium, organism hit count, and e-value significance
- **Multi-format output**: Generates CSV, TSV, and Excel reports with professional formatting and documentation
- **Parallel processing**: Utilizes multiple CPU threads for faster DIAMOND alignments
- **Flexible filtering**: Customizable thresholds for percent identity, chromosome/contig filtering, and alignment sensitivity
- **Detailed statistics**: Includes intron length calculations, organism diversity metrics, and confidence categorization

## Requirements

### Software Dependencies

- **Python 3.7+** with the following packages:
  - `pandas` - Data manipulation and analysis
  - `tqdm` - Progress bar visualization
  - `regex` - Enhanced regular expression support
  - `openpyxl` (optional) - Excel file generation

- **DIAMOND** - Fast protein sequence aligner
  - Must be installed and available in your system PATH
  - Download from: https://github.com/bbuchfink/diamond

### Input Files

1. **Proteome file** (`.faa`) - FASTA file containing all protein sequences for the organism
2. **Genome annotation** (`.gff`) - GFF3 format annotation file (RefSeq recommended)
3. **Reference database** (`.dmnd`) - Pre-built DIAMOND database of reference proteins

## Installation

1. Clone the repository:
```bash
git clone https://github.com/txdylan27/GenoMiss.git
cd GenoMiss
```

2. Install Python dependencies:
```bash
pip install pandas tqdm regex openpyxl
```

3. Install DIAMOND (if not already installed):
```bash
# On Ubuntu/Debian
sudo apt-get install diamond-aligner

# Or download from GitHub releases
wget http://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
sudo mv diamond /usr/local/bin/
```

## Usage

### Basic Usage

```bash
python GenoMiss.py \
  -p <proteome.faa> \
  -a <annotation.gff> \
  -db <reference.dmnd> \
  -o <output_directory>
```

### Example

```bash
python GenoMiss.py \
  -p apis_mellifera_protein.faa \
  -a apis_mellifera_annotation.gff \
  -db insecta_refseq_protein_db.dmnd \
  -o honeybee_results \
  -t 16 \
  -ds very-sensitive
```

### Command-Line Options

#### Required Arguments

| Option | Description |
|--------|-------------|
| `-p`, `--proteome` | Path to the `.faa` proteome file containing all protein sequences |
| `-a`, `--organism_annotation` | Path to the `.gff` genome annotation file (RefSeq format recommended) |
| `-db`, `--database` | Path to the `.dmnd` DIAMOND reference protein database |
| `-o`, `--output` | Output directory path for results and intermediate files |

#### Optional Arguments

| Option | Default | Description |
|--------|---------|-------------|
| `-t`, `--num_threads` | Half of available CPU cores | Number of threads for DIAMOND alignment (positive integer) |
| `-i`, `--identity_cutoff` | 0.00 | Percent identity cutoff for filtering fused gene alignments (0.0-100.0) |
| `-xf`, `--xfilter` | 5 | Minimum number of genes required per chromosome/contig for analysis (filters out small contigs) |
| `-ds`, `--diamond_sensitivity` | None | DIAMOND sensitivity mode: `fast`, `mid-sensitive`, `sensitive`, `more-sensitive`, `very-sensitive`, or `ultra-sensitive` |

## Output Files

The tool generates multiple output files in the specified output directory:

### Primary Output Files

| File | Format | Description |
|------|--------|-------------|
| `misannotation_results.xlsx` | Excel | Multi-sheet workbook with formatted results, including Top Hits, All Fused Hits, Control Hits, and Scoring Methodology |
| `fused_hits.csv` | CSV | All fused gene hits sorted by composite score with comprehensive metadata |
| `high_confidence_hits.csv` | CSV | Subset of fused hits with composite score ≥ 70 |
| `control_hits.csv` | CSV | Alignment results for individual unfused gene parts |
| `fused_hits.tsv` | TSV | Tab-separated version of all fused hits (terminal-friendly) |
| `control_hits.tsv` | TSV | Tab-separated version of control hits |

### Intermediate Files

| File | Description |
|------|-------------|
| `args.log` | Command-line arguments used for the analysis |
| `genome_wide_positive_hits.faa` | FASTA file containing gene parts that were components of positive fused hits |
| `<chrom>/<strand>/` | Per-chromosome/strand intermediate files and DIAMOND results |

## Scoring Methodology

The tool uses a composite scoring system (0-100 scale) that combines five components:

| Component | Weight | Description |
|-----------|--------|-------------|
| **Query Coverage** | 50% | Percentage of the fused protein sequence covered by the alignment |
| **Bit Score Improvement** | 25% | Relative improvement of fused protein bit score compared to individual gene part bit scores |
| **Overlap Equilibrium** | 5% | Balance of alignment between the two gene parts (favors alignments spanning both genes evenly) |
| **Organism Count** | 15% | Number of different organisms with hits to the fused gene (higher diversity increases confidence) |
| **E-value** | 5% | Statistical significance of the alignment (lower e-values increase confidence) |

### Score Interpretation

- **High Confidence (≥70)**: Strong evidence of misannotation; fused gene shows substantial improvement
- **Medium Confidence (40-69)**: Moderate evidence; warrants manual inspection
- **Low Confidence (<40)**: Weak evidence; may represent false positives or edge cases

## Algorithm Overview

1. **Genome Map Construction**
   - Parse GFF annotation to identify protein-coding genes
   - Build strand-specific directed graphs representing gene adjacency
   - Handle overlapping genes using branching paths

2. **Fusion Generation**
   - Traverse genome graph for each chromosome and strand
   - Create all possible isoform combinations between adjacent genes
   - Store metadata (gene IDs, product IDs, sequence lengths)

3. **DIAMOND Alignment**
   - Align fused proteins against reference database
   - Align individual gene parts (control group)
   - Filter hits based on overlap threshold (10 amino acids on each side of fusion junction)

4. **Score Calculation**
   - Calculate composite scores for all fused hits
   - Compute intron lengths from GFF coordinates
   - Count unique organisms per gene pair

5. **Output Generation**
   - Reorganize columns for readability
   - Generate summary statistics
   - Export in multiple formats with professional formatting

## Project Structure

```
GenoMiss/
├── GenoMiss.py                 # Main program and CLI interface
├── GenomeMap.py                # Genome graph data structures
├── scoring.py                  # Composite scoring system
├── output_formatter.py         # Multi-format report generation
├── input_files/                # Example input files (if provided)
└── README.md                   # This file
```

## Performance Considerations

- **Chromosome filtering** (`-xf`): Increase threshold to skip small contigs and reduce runtime
- **Thread count** (`-t`): More threads accelerate DIAMOND but increase memory usage
- **Sensitivity mode** (`-ds`): Higher sensitivity increases accuracy but significantly increases runtime
  - `fast`: ~10x faster than default, suitable for large-scale screening
  - `very-sensitive`: ~10x slower than default, recommended for publication-quality results

## Citation

If you use this tool in your research, please cite:

```
[Citation information to be added]
```

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Authors

- Dylan Ulloa
- David Bellini

## Acknowledgments

- DIAMOND alignment tool by Benjamin Buchfink

## Contact

For questions, issues, or feature requests, please open an issue on GitHub or contact the authors.
