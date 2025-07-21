# combine_merged_and_unmerged_forwards.py
from Bio import SeqIO
import os

def combine_merged_and_unmerged_forwards(merged_fasta_path, unmerged_pairs_fasta_path, final_output_fasta_path):
    """
    Combines successfully merged ASVs with the forward ASVs from unmerged pairs.

    Args:
        merged_fasta_path (str): Path to the FASTA file containing successfully merged ASVs (from bbmerge's 'out=').
        unmerged_pairs_fasta_path (str): Path to the FASTA file containing unmerged paired ASVs (from bbmerge's 'outunmerged=').
        final_output_fasta_path (str): Path for the new combined FASTA file.
    """
    final_records = []
    seen_ids = set() # To ensure uniqueness if an ID could somehow appear in both files (unlikely with mix=f)

    print(f"Collecting successfully merged ASVs from: {merged_fasta_path}")
    try:
        for record in SeqIO.parse(merged_fasta_path, "fasta"):
            if record.id not in seen_ids:
                final_records.append(record)
                seen_ids.add(record.id)
        print(f"Collected {len(final_records)} merged ASVs.")
    except FileNotFoundError:
        print(f"Error: Merged FASTA file '{merged_fasta_path}' not found. Exiting.")
        return
    except Exception as e:
        print(f"An error occurred while reading merged FASTA: {e}")
        return

    print(f"Collecting forward ASVs from unmerged pairs in: {unmerged_pairs_fasta_path}")
    num_unmerged_f_added = 0
    try:
        # This file should contain interleaved forward and reverse reads for unmerged pairs.
        # We only want the forward read.
        # We assume forward reads have '_F_' in their ID (e.g., ITS2_Plants_paired_F_XYZ).
        for record in SeqIO.parse(unmerged_pairs_fasta_path, "fasta"):
            if "_F_" in record.id and record.id not in seen_ids:
                final_records.append(record)
                seen_ids.add(record.id)
                num_unmerged_f_added += 1
            # We implicitly skip reverse reads (those with '_R_') and any other unexpected IDs.

        print(f"Added {num_unmerged_f_added} forward ASVs from unmerged pairs.")
    except FileNotFoundError:
        print(f"Warning: Unmerged pairs FASTA file '{unmerged_pairs_fasta_path}' not found. Skipping addition of unmerged forwards.")
    except Exception as e:
        print(f"An error occurred while reading unmerged pairs FASTA: {e}")

    print(f"Total sequences in final output: {len(final_records)}")

    try:
        SeqIO.write(final_records, final_output_fasta_path, "fasta")
        print(f"Combined ASVs saved to: {final_output_fasta_path}")
    except Exception as e:
        print(f"An error occurred while writing final FASTA: {e}")


if __name__ == "__main__":
    # --- Configuration ---
    # !!! IMPORTANT: Ensure these paths are correct and match your bbmerge output !!!
    MERGED_ASVS_FILE = "merged_ITS2_asvs.fasta"            # Output from bbmerge's 'out=' parameter
    UNMERGED_PAIRS_FILE = "unmerged_ITS2_pairs.fasta"      # Output from bbmerge's 'outunmerged=' parameter
    FINAL_COMBINED_ASVS_FILE = "final_combined_asvs.fasta" # The new file this script will create

    combine_merged_and_unmerged_forwards(MERGED_ASVS_FILE, UNMERGED_PAIRS_FILE, FINAL_COMBINED_ASVS_FILE)
