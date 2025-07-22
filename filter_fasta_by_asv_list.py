# filter_fasta_by_asv_list.py
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys # Import sys for stderr printing

def filter_fasta_files(asv_list_path, forward_fasta_path, reverse_fasta_path,
                       output_forward_fasta_path, output_reverse_fasta_path):
    """
    Filters FASTA files, removing sequences whose IDs are in a given list.

    Args:
        asv_list_path (str): Path to the text file containing ASV names to remove.
        forward_fasta_path (str): Path to the input forward ASV FASTA file.
        reverse_fasta_path (str): Path to the input reverse ASV FASTA file.
        output_forward_fasta_path (str): Path for the output filtered forward FASTA file.
        output_reverse_fasta_path (str): Path for the output filtered reverse FASTA file.
    """
    asvs_to_remove = set()
    try:
        with open(asv_list_path, 'r') as f:
            for line in f:
                stripped_line = line.strip()
                asvs_to_remove.add(stripped_line)
        print(f"Loaded {len(asvs_to_remove)} ASVs to remove from {asv_list_path}", file=sys.stderr)
        # Debug: Print a few ASVs to remove to verify format
        print(f"Sample ASVs to remove: {list(asvs_to_remove)[:5]}", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: ASV list file '{asv_list_path}' not found. Cannot filter FASTA files.", file=sys.stderr)
        return
    except Exception as e:
        print(f"An error occurred while reading the ASV list: {e}", file=sys.stderr)
        return

    # Filter Forward FASTA
    filtered_forward_records = []
    try:
        total_forward_read = 0
        removed_forward_read = 0
        for record in SeqIO.parse(forward_fasta_path, "fasta"):
            total_forward_read += 1
            # Debug: Print current ASV ID and check if it's in the set
            if total_forward_read % 1000 == 0: # Print for every 1000th ASV to avoid too much output
                print(f"Processing forward ASV: {record.id}", file=sys.stderr)
                # For a specific ASV like 18S_Euk_paired_F_8338, you can add a direct check:
                # if record.id == "18S_Euk_paired_F_8338":
                #     print(f"DEBUG: Found 18S_Euk_paired_F_8338. Is in asvs_to_remove? {record.id in asvs_to_remove}", file=sys.stderr)

            if record.id in asvs_to_remove: # Changed to 'in' for clarity in the condition
                removed_forward_read += 1
            else:
                filtered_forward_records.append(record)

        SeqIO.write(filtered_forward_records, output_forward_fasta_path, "fasta")
        print(f"Processed {total_forward_read} forward ASVs. Removed {removed_forward_read}. Saved remaining {len(filtered_forward_records)} to {output_forward_fasta_path}", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: Forward FASTA file '{forward_fasta_path}' not found.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred while filtering forward FASTA: {e}", file=sys.stderr)

    # Filter Reverse FASTA
    filtered_reverse_records = []
    try:
        total_reverse_read = 0
        removed_reverse_read = 0
        for record in SeqIO.parse(reverse_fasta_path, "fasta"):
            total_reverse_read += 1
            derived_forward_id = record.id.replace('_R_', '_F_') # Match reverse ASVs to forward IDs in list
            # Debug: Print current reverse ASV ID and derived forward ID
            if total_reverse_read % 1000 == 0: # Print for every 1000th ASV
                print(f"Processing reverse ASV: {record.id}, Derived forward ID: {derived_forward_id}", file=sys.stderr)
                # if derived_forward_id == "18S_Euk_paired_F_8338":
                #     print(f"DEBUG: Found reverse for 18S_Euk_paired_F_8338. Is in asvs_to_remove? {derived_forward_id in asvs_to_remove}", file=sys.stderr)

            if derived_forward_id in asvs_to_remove: # Check if the derived ID is in the set to remove
                removed_reverse_read += 1
            else:
                filtered_reverse_records.append(record)

        SeqIO.write(filtered_reverse_records, output_reverse_fasta_path, "fasta")
        print(f"Processed {total_reverse_read} reverse ASVs. Removed {removed_reverse_read}. Saved remaining {len(filtered_reverse_records)} to {output_reverse_fasta_path}", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: Reverse FASTA file '{reverse_fasta_path}' not found.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred while filtering reverse FASTA: {e}", file=sys.stderr)

if __name__ == "__main__":
    # --- Configuration ---
    # !!! IMPORTANT: Replace these with your actual file paths !!!
    ASV_LIST_FILE = "unassigned_asvs_to_remove.txt" # Output from the first script
    FORWARD_ASV_FASTA = "clve5i5mn0003l50gam95h4jx-ITS2_Plants-paired_F.fasta"   # Your original forward ASV FASTA
    REVERSE_ASV_FASTA = "clve5i5mn0003l50gam95h4jx-ITS2_Plants-paired_R.fasta"   # Your original reverse ASV FASTA

    # Names for the new filtered output files
    OUTPUT_FORWARD_FASTA = "filtered_forward_asvs.fasta"
    OUTPUT_REVERSE_FASTA = "filtered_reverse_asvs.fasta"
    # ---------------------

    filter_fasta_files(ASV_LIST_FILE, FORWARD_ASV_FASTA, REVERSE_ASV_FASTA,
                       OUTPUT_FORWARD_FASTA, OUTPUT_REVERSE_FASTA)
