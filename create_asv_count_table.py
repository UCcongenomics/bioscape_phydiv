# create_asv_count_table.py
import pandas as pd
from Bio import SeqIO
import os

def create_final_asv_table(counts_file, taxonomy_file, combined_fasta, output_file,
                           plate_name=None, marker_name=None):
    """
    Creates a combined ASV table with read counts per sample, taxonomic identity, ASV sequence,
    and optionally adds columns for plate and marker names.
    ASVs without a matching taxonomic entry will now be excluded from the final table.

    Args:
        counts_file (str): Path to the ASV counts file (.asv).
        taxonomy_file (str): Path to the taxonomic information file (.txt).
        combined_fasta (str): Path to the FASTA file containing the final ASV names and sequences to include.
        output_file (str): Path for the new tab-delimited output table.
        plate_name (str, optional): Name of the plate directory (e.g., 'PlateA'). Defaults to None.
        marker_name (str, optional): Name of the marker directory (e.g., 'ITS2_Plants'). Defaults to None.
    """
    print(f"--- Starting final ASV table creation for {plate_name}/{marker_name} ---")

    # --- 1. Load Taxonomy Data ---
    print(f"Loading taxonomic data from: {taxonomy_file}")
    try:
        if not os.path.exists(taxonomy_file):
            print(f"Error: Taxonomy file '{taxonomy_file}' not found. Skipping table creation.")
            return

        taxonomy_df = pd.read_csv(taxonomy_file, sep='\t', usecols=['readname', 'taxonomic_path'])
        taxonomy_dict = taxonomy_df.set_index('readname')['taxonomic_path'].to_dict()
        print(f"Loaded {len(taxonomy_dict)} taxonomic entries.")
    except Exception as e:
        print(f"An error occurred while loading taxonomy data from '{taxonomy_file}': {e}")
        return

    # --- 2. Load ASV Counts Data ---
    print(f"Loading ASV counts data from: {counts_file}")
    try:
        if not os.path.exists(counts_file):
            print(f"Error: ASV counts file '{counts_file}' not found. Skipping table creation.")
            return

        asv_counts_df = pd.read_csv(counts_file, sep='\t')

        # Rename the ASV ID column for consistency
        if 'ITS2_Plants_seq_number' in asv_counts_df.columns:
            asv_counts_df.rename(columns={'ITS2_Plants_seq_number': 'ASV_Name'}, inplace=True)
        else:
            potential_asv_col = asv_counts_df.columns[0]
            if not potential_asv_col.startswith("K") and not potential_asv_col.startswith("S"):
                 asv_counts_df.rename(columns={potential_asv_col: 'ASV_Name'}, inplace=True)
            else:
                 print(f"Warning: Could not automatically identify ASV Name column. First column is '{potential_asv_col}'. Please verify input .asv file structure if results are unexpected.")
                 asv_counts_df.rename(columns={potential_asv_col: 'ASV_Name'}, inplace=True)


        # Define columns to be explicitly dropped from the ASV counts file if they exist
        columns_to_drop = [
            'score', 'forward_mismatch', 'reverse_mismatch',
            'divergence', 'chi2', 'tree_number', 'node_number',
            'forward_length', 'reverse_length', 'taxonomic_path' # Ensure taxonomic_path from .asv is dropped if present
        ]
        asv_counts_df.drop(columns=columns_to_drop, errors='ignore', inplace=True)

        metadata_cols_remaining = ['ASV_Name', 'forward_sequence'] # Removed 'taxonomic_path' as it's being mapped from .txt
        metadata_cols_existing = [col for col in metadata_cols_remaining if col in asv_counts_df.columns]

        sample_cols = [col for col in asv_counts_df.columns if col not in metadata_cols_existing]

        asv_counts_subset_df = asv_counts_df[['ASV_Name'] + sample_cols].copy()
        asv_counts_subset_df['Taxonomy'] = asv_counts_subset_df['ASV_Name'].map(taxonomy_dict)

        # --- MODIFICATION START ---
        # Instead of labeling 'Unassigned', we will drop ASVs that did not have a matching taxonomic entry.
        initial_asv_count = len(asv_counts_subset_df)
        asv_counts_subset_df.dropna(subset=['Taxonomy'], inplace=True)
        dropped_asv_count = initial_asv_count - len(asv_counts_subset_df)

        if dropped_asv_count > 0:
            print(f"Removed {dropped_asv_count} ASVs from the counts data due to missing taxonomic entries.")
        else:
            print("No ASVs were removed from the counts data due to missing taxonomic entries.")
        # --- MODIFICATION END ---

        asv_counts_subset_df.set_index('ASV_Name', inplace=True)
        print(f"Loaded {len(asv_counts_subset_df)} ASVs with counts (after taxonomy check).")

    except Exception as e:
        print(f"An error occurred while loading ASV counts data from '{counts_file}': {e}")
        return

    # --- 3. Get ASV Names and Sequences from final_combined_asvs.fasta ---
    print(f"Collecting ASV names and sequences from: {combined_fasta}")
    included_asv_data = {}
    try:
        if not os.path.exists(combined_fasta):
            print(f"Error: Combined FASTA file '{combined_fasta}' not found. Skipping table creation.")
            return

        for record in SeqIO.parse(combined_fasta, "fasta"):
            included_asv_data[record.id] = str(record.seq)
        print(f"Found {len(included_asv_data)} ASVs with sequences in the combined FASTA file.")
    except Exception as e:
        print(f"An error occurred while reading combined FASTA '{combined_fasta}': {e}")
        return

    # --- 4. Filter and Prepare Final Table ---
    print("Filtering and preparing final table with sequence column...")

    final_table_df = asv_counts_subset_df.loc[
        asv_counts_subset_df.index.intersection(included_asv_data.keys())
    ].copy()

    final_table_df['Sequence'] = final_table_df.index.map(included_asv_data)
    final_table_df.reset_index(inplace=True)

    # --- Add Plate and Marker Columns ---
    if plate_name:
        final_table_df['Plate'] = plate_name
    if marker_name:
        final_table_df['Marker'] = marker_name

    # Reorder columns: Plate, Marker, ASV_Name, Sequence, Taxonomy, then all sample columns
    desired_cols = []
    if plate_name: desired_cols.append('Plate')
    if marker_name: desired_cols.append('Marker')
    desired_cols.extend(['ASV_Name', 'Sequence', 'Taxonomy'])
    desired_cols.extend(sample_cols)

    for col in desired_cols:
        if col not in final_table_df.columns:
            final_table_df[col] = None

    final_table_df = final_table_df[desired_cols]


    print(f"Final table contains {len(final_table_df)} ASVs matching the combined FASTA.")

    # --- 5. Write to Tab-Delimited File ---
    try:
        final_table_df.to_csv(output_file, sep='\t', index=False)
        print(f"Final ASV table saved to: {output_file}")
    except Exception as e:
        print(f"An error occurred while writing the final table to '{output_file}': {e}")

# (The __main__ block is still commented out or can be deleted)
