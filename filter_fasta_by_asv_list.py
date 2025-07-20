# filter_fasta_by_asv_list.py
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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
                asvs_to_remove.add(line.strip())
        print(f"Loaded {len(asvs_to_remove)} ASVs to remove from {asv_list_path}")
    except FileNotFoundError:
        print(f"Error: ASV list file '{asv_list_path}' not found. Cannot filter FASTA files.")
        return
    except Exception as e:
        print(f"An error occurred while reading the ASV list: {e}")
        return

    # Filter Forward FASTA
    filtered_forward_records = []
    try:
        total_forward_read = 0
        removed_forward_read = 0
        for record in SeqIO.parse(forward_fasta_path, "fasta"):
            total_forward_read += 1
            if record.id not in asvs_to_remove:
                filtered_forward_records.append(record)
            else:
                removed_forward_read += 1
        SeqIO.write(filtered_forward_records, output_forward_fasta_path, "fasta")
        print(f"Processed {total_forward_read} forward ASVs. Removed {removed_forward_read}. Saved remaining {len(filtered_forward_records)} to {output_forward_fasta_path}")
    except FileNotFoundError:
        print(f"Error: Forward FASTA file '{forward_fasta_path}' not found.")
    except Exception as e:
        print(f"An error occurred while filtering forward FASTA: {e}")

    # Filter Reverse FASTA
    filtered_reverse_records = []
    try:
        total_reverse_read = 0
        removed_reverse_read = 0
        for record in SeqIO.parse(reverse_fasta_path, "fasta"):
            total_reverse_read += 1
            if record.id not in asvs_to_remove:
                filtered_reverse_records.append(record)
            else:
                removed_reverse_read += 1
        SeqIO.write(filtered_reverse_records, output_reverse_fasta_path, "fasta")
        print(f"Processed {total_reverse_read} reverse ASVs. Removed {removed_reverse_read}. Saved remaining {len(filtered_reverse_records)} to {output_reverse_fasta_path}")
    except FileNotFoundError:
        print(f"Error: Reverse FASTA file '{reverse_fasta_path}' not found.")
    except Exception as e:
        print(f"An error occurred while filtering reverse FASTA: {e}")

if __name__ == "__main__":
    # --- Configuration ---
    # !!! IMPORTANT: Replace these with your actual file paths !!!
    ASV_LIST_FILE = "unassigned_asvs_to_remove.txt" # Output from the first script
    FORWARD_ASV_FASTA = "clve5i5mn0003l50gam95h4jx-ITS2_Plants-paired_F.fasta"   # Your original forward ASV FASTA
    REVERSE_ASV_FASTA = "clve5i5mn0003l50gam95h4jx-ITS2_Plants-paired_R.fasta"   # Your original reverse ASV FASTA
    
    # Names for the new filtered output files
    OUTPUT_FORWARD_FASTA = "filtered_ITS2_forward_asvs.fasta"
    OUTPUT_REVERSE_FASTA = "filtered_ITS2_reverse_asvs.fasta"
    # ---------------------

    filter_fasta_files(
        ASV_LIST_FILE,
        FORWARD_ASV_FASTA,
        REVERSE_ASV_FASTA,
        OUTPUT_FORWARD_FASTA,
        OUTPUT_REVERSE_FASTA
    )
