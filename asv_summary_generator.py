# asv_summary_generator.py
import pandas as pd
import os
import numpy as np

def generate_asv_summary(file_path, marker_name="UnknownMarker"):
    """
    Generates a histogram of ASV read counts and a summary of unique taxonomies
    from a final ASV table. The output files are named using the marker_name.

    Args:
        file_path (str): Path to the ASV table file (e.g., PlateA_ITS1_Fungi_final_asv_table.tsv).
        marker_name (str): The name of the marker (e.g., '16S_Bacteria', '18S_Euk') for naming output files.
    """
    print(f"\n--- Generating ASV Summary for {marker_name} from: {file_path} ---")

    if not os.path.exists(file_path):
        print(f"Error: The file '{file_path}' was not found. Skipping summary for {marker_name}.")
        return

    try:
        df = pd.read_csv(file_path, sep='\t')
        print(f"Successfully loaded {len(df)} rows for {marker_name}.")

        # --- Part 1: Generate Histogram of ASV Read Counts ---
        print(f"--- Generating Histogram of ASV Read Counts for {marker_name} ---")
        meta_cols = ['Plate', 'Marker', 'ASV_Name', 'Sequence', 'Taxonomy']
        sample_cols = [col for col in df.columns if col not in meta_cols]

        if not sample_cols:
            print(f"No sample (read count) columns identified for {marker_name}. Cannot calculate total reads for histogram.")
        else:
            for col in sample_cols:
                df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)

            df['Total_Reads'] = df[sample_cols].sum(axis=1)

            bins = [0, 1, 5, 10, 50, 100, np.inf]
            labels = ['1', '2-5', '6-10', '11-50', '51-100', '>100']

            histogram_table = pd.cut(df['Total_Reads'], bins=bins, labels=labels, include_lowest=True, right=True).value_counts().reindex(labels, fill_value=0)
            histogram_df = histogram_table.reset_index()
            histogram_df.columns = ['Read_Count_Range', 'Number_of_ASVs']

            print(f"\nHistogram Table of ASV Read Counts for {marker_name}:")
            print(histogram_df.to_string(index=False))

            # Output path for histogram CSV, placed in the same directory as the input ASV table
            output_histogram_csv_path = os.path.join(os.path.dirname(file_path), f"{marker_name}_asv_read_count_histogram.csv")
            histogram_df.to_csv(output_histogram_csv_path, index=False)
            print(f"Histogram table saved to: {output_histogram_csv_path}")

        # --- Part 2: Count Unique Taxonomy ---
        print(f"\n--- Counting Unique Taxonomy for {marker_name} ---")
        if 'Taxonomy' not in df.columns:
            print(f"Error: 'Taxonomy' column not found for {marker_name}. Cannot count unique taxonomies.")
        else:
            taxonomy_counts = df['Taxonomy'].value_counts()
            taxonomy_counts_df = taxonomy_counts.reset_index()
            taxonomy_counts_df.columns = ['Taxonomy_Path', 'Count']

            print(f"\nCounts of Unique Taxonomy Paths (Top 10) for {marker_name}:")
            print(taxonomy_counts_df.head(10).to_string(index=False))

            # Output path for taxonomy counts CSV, placed in the same directory as the input ASV table
            output_taxonomy_csv_path = os.path.join(os.path.dirname(file_path), f"{marker_name}_unique_taxonomy_counts.csv")
            taxonomy_counts_df.to_csv(output_taxonomy_csv_path, index=False)
            print(f"Unique taxonomy counts saved to: {output_taxonomy_csv_path}")

    except Exception as e:
        print(f"An error occurred while processing the file for {marker_name}: {e}")

# Example of how you might integrate this into your main loop script:
# (This part is commented out as it's for illustration, integrate it into your existing pipeline script)
#
# if __name__ == "__main__":
#     base_dir = "/projects/RachelMeyer/BioSCAPE/PlateA" # Adjust to your base directory
#     # List all your markers here that you process in your pipeline
#     markers = ["16S_Bacteria", "18S_Euk", "ITS1_Fungi"]
#
#     for marker in markers:
#         # Assuming your final ASV tables are named consistently, e.g., PlateA_16S_Bacteria_final_asv_table.tsv
#         final_asv_table_filename = f"PlateA_{marker}_final_asv_table.tsv"
#         final_asv_table_path = os.path.join(base_dir, marker, "paired", final_asv_table_filename)
#
#         # Call the function for each marker's final ASV table
#         generate_asv_summary(final_asv_table_path, marker)
