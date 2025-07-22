# combine_merged_and_unmerged_forwards.py
from Bio import SeqIO
import os

def combine_asv_fastas(merged_fasta_file, unmerged_fasta_file, output_combined_fasta):
    """
    Combines successfully merged ASVs with the forward reads of unmerged pairs.

    Args:
        merged_fasta_file (str): Path to the FASTA file containing successfully merged ASVs.
        unmerged_fasta_file (str): Path to the FASTA file containing unmerged read pairs.
        output_combined_fasta (str): Path for the new combined output FASTA file.
    """
    final_asvs = {} # Use a dictionary to store unique ASVs and avoid duplicates

    # Process merged ASVs
    print(f"Processing merged ASVs from: {merged_fasta_file}")
    try:
        if not os.path.exists(merged_fasta_file):
            print(f"Error: Merged FASTA file '{merged_fasta_file}' not found. Skipping processing.")
            return # Exit function if crucial input is missing

        for record in SeqIO.parse(merged_fasta_file, "fasta"):
            if record.id not in final_asvs:
                final_asvs[record.id] = record
        print(f"Added {len(final_asvs)} merged ASVs.")
    except Exception as e:
        print(f"An error occurred while processing merged FASTA '{merged_fasta_file}': {e}")
        return

    # Process unmerged forward reads
    print(f"Processing unmerged reads from: {unmerged_fasta_file}")
    try:
        if not os.path.exists(unmerged_fasta_file):
            print(f"Error: Unmerged FASTA file '{unmerged_fasta_file}' not found. Skipping processing.")
            return # Exit function if crucial input is missing

        unmerged_forward_count = 0
        for record in SeqIO.parse(unmerged_fasta_file, "fasta"):
            # Assuming forward reads have '_F_' in their ID for unmerged pairs,
            # but more robustly, we just need the first sequence of the pair.
            # If the file contains concatenated F/R sequences, we might need
            # to split them. For now, assuming distinct F records.
            # Let's check if it's explicitly a forward read if BBmerge splits them like that.
            # Or, if BBmerge's `outunmerged` just gives one FASTA per pair (first is F).
            # Based on standard usage, `outunmerged` often gives R1 and R2 records.
            # We want the R1 record (forward).
            if '_R_' not in record.id and '_2:' not in record.id: # Simple check for forward ID
                if record.id not in final_asvs:
                    final_asvs[record.id] = record
                    unmerged_forward_count += 1
        print(f"Added {unmerged_forward_count} unique unmerged forward ASVs.")
    except Exception as e:
        print(f"An error occurred while processing unmerged FASTA '{unmerged_fasta_file}': {e}")
        return

    # Write the combined unique ASVs to a new FASTA file
    print(f"Writing combined unique ASVs to: {output_combined_fasta}")
    try:
        with open(output_combined_fasta, "w") as output_handle:
            SeqIO.write(final_asvs.values(), output_handle, "fasta")
        print(f"Successfully created {output_combined_fasta} with {len(final_asvs)} unique ASVs.")
    except Exception as e:
        print(f"An error occurred while writing the combined FASTA: {e}")

# The __main__ block is now commented out or can be deleted,
# as this script will be called by the wrapper.
# if __name__ == "__main__":
#     # --- Configuration ---
#     MERGED_FASTA_FILE = "merged_ITS2_asvs.fasta"
#     UNMERGED_FASTA_FILE = "unmerged_ITS2_pairs.fasta"
#     OUTPUT_COMBINED_FASTA = "final_combined_asvs.fasta"
#
#     combine_asv_fastas(MERGED_FASTA_FILE, UNMERGED_FASTA_FILE, OUTPUT_COMBINED_FASTA)
