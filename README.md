# ASV Processing Pipeline

## Overview
This pipeline is designed to process Amplicon Sequence Variant (ASV) data generated from environmental samples, typically after initial bioinformatics steps like denoising (e.g., DADA2 output). It handles filtering of unassigned ASVs, merging paired-end reads, combining ASV sequences, generating a comprehensive ASV count table, and producing summary statistics including ASV read count histograms and unique taxonomy counts.

The pipeline is structured to iterate through multiple plates and genetic markers, ensuring consistent processing and organized output for each dataset.

## Features
* **Automated Processing:** Loops through specified plate and marker directories to process each dataset.
* **Unassigned ASV Removal:** Identifies and filters out ASVs without clear taxonomic assignment (based on configurable keywords), preventing them from appearing in final tables.
* **Read Merging:** Utilizes BBmerge to combine paired-end reads and handles unmerged reads appropriately.
* **Comprehensive ASV Table Generation:** Creates a final tab-separated value (TSV) table that includes ASV read counts per sample, their full taxonomic path, and their sequence.
* **Dynamic Metadata:** Automatically adds 'Plate' and 'Marker' columns to the final ASV table.
* **Summary Reports:** Generates two key summary files for each marker:
    * A histogram of ASV read counts, categorized into predefined bins.
    * A count of unique taxonomic paths.

## Prerequisites

### Software
* **Python 3.6+**: The pipeline is written in Python 3 and uses f-strings, which require Python 3.6 or newer.
* **BBMap Suite**: Specifically, `bbmerge.sh` is required for merging paired-end reads. Ensure its path is correctly configured in `run_pipeline_loop.py`.

### Python Libraries
Ensure you have the following Python libraries installed in your environment (preferably a `conda` or `venv` environment):
* `pandas`
* `numpy`
* `biopython`

You can install them using pip:
pip install pandas numpy biopython

# Setup

Organize your scripts:
Place all the Python scripts (i.e., run_pipeline_loop.py, create_asv_count_table.py, identify_unassigned_asvs.py, filter_fasta_by_asv_list.py, combine_merged_and_unmerged_forwards.py, and asv_summary_generator.py) in the same directory.

# Configure run_pipeline_loop.py:
Open run_pipeline_loop.py and modify the following variables in the "--- Configuration for the Loop Script ---" section:

BASE_BIO_DIR: Set this to the root directory containing your "PlateX" folders (e.g., /projects/RachelMeyer/BioSCAPE/).

BBMERGE_PATH: Crucially, update this to the absolute path of your bbmerge.sh executable (e.g., /home/rameyer/miniconda2/bin/bbmerge.sh).

# Ensure input file structure:
The pipeline expects your data to be organized hierarchically:

BASE_BIO_DIR/
├── PlateA/
│   ├── 16S_Bacteria/
│   │   └── paired/
│   │       ├── [unique_prefix]-16S_Bacteria-paired_F.asv
│   │       ├── [unique_prefix]-16S_Bacteria-paired.txt
│   │       ├── [unique_prefix]-16S_Bacteria-paired_F.fasta
│   │       └── [unique_prefix]-16S_Bacteria-paired_R.fasta
│   └── ITS1_Fungi/
│       └── paired/
│           ├── ...
└── PlateB/
    └── ...

Where [unique_prefix] is an identifier specific to your DADA2 output files (e.g., clve5i5mn0003l50gam95h4jx).

# Usage
Activate your Python environment:
It is highly recommended to use a dedicated Python environment (like a conda environment) where all prerequisites are installed.

conda activate my_bio_env
Or source activate my_bio_env (for older conda versions)

# Run the pipeline:
Navigate to the directory where you saved run_pipeline_loop.py and execute it:
python run_pipeline_loop.py

The script will then print its progress to the console as it processes each plate and marker.

# Output
For each paired directory within BASE_BIO_DIR/PlateX/MarkerY/, the pipeline will generate the following files:

unassigned_asvs_to_remove.txt: A list of ASV names identified as 'unassigned' based on keywords defined in identify_unassigned_asvs.py. These ASVs will be excluded from further processing.

filtered_forward_asvs.fasta: Forward ASV sequences after removing unassigned ASVs.

filtered_reverse_asvs.fasta: Reverse ASV sequences after removing unassigned ASVs.

merged_asvs.fasta: ASV sequences successfully merged by BBmerge.

unmerged_pairs.fasta: Unmerged forward ASV sequences (from filtered_forward_asvs.fasta).

final_combined_asvs.fasta: A FASTA file containing all ASV sequences that were either successfully merged or were unmerged forwards (from unmerged_pairs.fasta), after filtering.

[PlateName]_[MarkerName]_final_asv_table.tsv: The main output ASV table. This tab-separated file contains ASV names, sequences, taxonomic assignments (without "Unassigned" entries), and read counts for each sample.

[MarkerName]_asv_read_count_histogram.csv: A CSV file summarizing the distribution of ASV total read counts into predefined bins.

[MarkerName]_unique_taxonomy_counts.csv: A CSV file listing each unique full taxonomic path found in the final ASV table and the count of ASVs assigned to it.

# Script Breakdown
run_pipeline_loop.py: The main orchestration script. It traverses plate and marker directories and calls the other modular scripts in the correct sequence.

identify_unassigned_asvs.py: Reads the DADA2 .txt statistics file and creates a list of ASV names to be considered 'unassigned' based on taxonomic keywords or single-domain assignments.

filter_fasta_by_asv_list.py: Takes the list of unassigned ASVs and filters corresponding ASV sequences from the raw forward and reverse FASTA files.

combine_merged_and_unmerged_forwards.py: Combines the ASV sequences that were successfully merged by BBmerge with the forward reads that could not be merged, creating a single comprehensive FASTA file of all remaining ASVs.

create_asv_count_table.py: Generates the final ASV abundance table, integrating read counts from the .asv file, taxonomic information from the .txt file, and ASV sequences from the combined FASTA. Crucially, it now excludes ASVs that lack a taxonomic entry.

asv_summary_generator.py: The new script that reads the _final_asv_table.tsv and generates the ASV read count histogram and the unique taxonomy counts.




