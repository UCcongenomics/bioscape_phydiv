# run_pipeline_loop.py
import os
import glob
import subprocess
import sys

# Import the functions from your modularized scripts
# Make sure these script files (filter_fasta_by_asv_list.py, etc.) are in the same directory
# as this run_pipeline_loop.py script, or modify the import paths.
try:
    from filter_fasta_by_asv_list import filter_fasta_files
    from combine_merged_and_unmerged_forwards import combine_asv_fastas
    from create_asv_count_table import create_final_asv_table
    from identify_unassigned_asvs import identify_unassigned_asvs
    from asv_summary_generator import generate_asv_summary # <--- NEW IMPORT: For histogram and taxonomy summary
except ImportError as e:
    print(f"Error importing one of the pipeline scripts: {e}")
    print("Please ensure all pipeline Python scripts are in the same directory as 'run_pipeline_loop.py'.")
    sys.exit(1)


# --- Configuration for the Loop Script ---
BASE_BIO_DIR = "/projects/RachelMeyer/BioSCAPE/"
# Path to your BBMap bbmerge.sh script
BBMERGE_PATH = "/home/rameyer/miniconda2/bin/bbmerge.sh" # !!! IMPORTANT: CHANGE THIS TO YOUR ACTUAL BBmerge.sh PATH !!!

# Suffixes for input files to find the unique prefix (e.g., ITS2_Plants)
ASV_COUNTS_SUFFIX = "_F.asv"
TAXONOMY_SUFFIX = ".txt"
FORWARD_FASTA_SUFFIX = "_F.fasta"
REVERSE_FASTA_SUFFIX = "_R.fasta"

def find_unique_prefix(directory, marker_name, suffix):
    """
    Finds a unique file prefix (e.g., 'clve5i5mn0003l50gam95h4jx-ITS2_Plants-paired')
    in a given directory based on a file suffix and marker name.
    The marker_name is used to filter files.
    """
    # Look for files like 'SOME_ID-MARKER_NAME-paired_F.asv'
    pattern = os.path.join(directory, f"*-{marker_name}-paired{suffix}")
    files = glob.glob(pattern)
    if not files:
        # Fallback: if specific marker name not in prefix, look for any *paired* file
        # This fallback is less specific but might be needed if the unique_prefix
        # doesn't strictly contain the 'marker_name' string.
        pattern_fallback = os.path.join(directory, f"*-paired{suffix}")
        files_fallback = glob.glob(pattern_fallback)
        if not files_fallback:
            return None, f"No files matching pattern '{pattern}' or '{pattern_fallback}' found."
        files = files_fallback # Use fallback if specific pattern fails

    # Extract the prefix before the suffix
    # Example: 'clve5i5mn0003l50gam95h4jx-ITS2_Plants-paired_F.asv' -> 'clve5i5mn0003l50gam95h4jx-ITS2_Plants-paired'
    first_file_name = os.path.basename(files[0])
    prefix = first_file_name.replace(suffix, '')
    return prefix, None # Return prefix and no error

def run_pipeline_for_dataset(plate_name, marker_name, data_dir, bbmerge_path):
    """
    Runs the full ASV processing pipeline for a single marker dataset.

    Args:
        plate_name (str): Name of the plate directory (e.g., 'PlateA').
        marker_name (str): Name of the marker directory (e.g., 'ITS2_Plants').
        data_dir (str): Full path to the 'paired' data directory for this marker.
        bbmerge_path (str): Path to the bbmerge.sh executable.
    """
    print(f"\n--- Processing Dataset: Plate='{plate_name}', Marker='{marker_name}' ---")

    if not os.path.isdir(data_dir):
        print(f"Skipping {plate_name}/{marker_name}: Data directory '{data_dir}' not found.")
        return

    # Discover the unique prefix (e.g., 'clve5i5mn0003l50gam95h4jx-ITS2_Plants-paired')
    unique_prefix, error_msg = find_unique_prefix(data_dir, marker_name, ASV_COUNTS_SUFFIX)
    if not unique_prefix:
        print(f"Warning: {error_msg} in {data_dir}. Skipping this dataset.")
        return

    print(f"Identified unique file prefix: {unique_prefix}")

    # Define all input and output paths for the current marker
    original_forward_fasta = os.path.join(data_dir, f"{unique_prefix}{FORWARD_FASTA_SUFFIX}")
    original_reverse_fasta = os.path.join(data_dir, f"{unique_prefix}{REVERSE_FASTA_SUFFIX}")
    asv_counts_file = os.path.join(data_dir, f"{unique_prefix}{ASV_COUNTS_SUFFIX}")
    taxonomy_file = os.path.join(data_dir, f"{unique_prefix}{TAXONOMY_SUFFIX}")

    # Intermediate files
    unassigned_asvs_to_remove = os.path.join(data_dir, "unassigned_asvs_to_remove.txt") # <--- Path for the generated list
    filtered_forward_fasta = os.path.join(data_dir, "filtered_forward_asvs.fasta")
    filtered_reverse_fasta = os.path.join(data_dir, "filtered_reverse_asvs.fasta")

    merged_fasta = os.path.join(data_dir, "merged_asvs.fasta")
    unmerged_pairs_fasta = os.path.join(data_dir, "unmerged_pairs.fasta")
    final_combined_fasta = os.path.join(data_dir, "final_combined_asvs.fasta")

    # Final output file name (INCLUDING PLATE AND MARKER NAMES)
    final_output_table = os.path.join(data_dir, f"{plate_name}_{marker_name}_final_asv_table.tsv")

    # --- NEW STEP: 0. Identify Unassigned ASVs ---
    print("\n--- Running Step 0: Identifying Unassigned ASVs ---")
    # The taxonomy_file (e.g., clve5i5mn0003l50gam95h4jx-ITS2_Plants-paired.txt)
    # is the stats file with 'readname' and 'taxonomic_path' columns.
    identify_unassigned_asvs(taxonomy_file, unassigned_asvs_to_remove)
    # Note: identify_unassigned_asvs.py will create the 'unassigned_asvs_to_remove.txt' file.
    # If no unassigned ASVs are found, the file will be empty, which is fine for the next step.
    # If the stats_file_path is not found, identify_unassigned_asvs will print an error and return.

    # --- Step 1: Filter FASTA by ASV List ---
    print("\n--- Running Step 1: Filtering FASTA files ---")
    # This step will now use the generated unassigned_asvs_to_remove.txt
    filter_fasta_files(unassigned_asvs_to_remove, original_forward_fasta, original_reverse_fasta,
                       filtered_forward_fasta, filtered_reverse_fasta)

    # --- Step 2: Run BBmerge ---
    print("\n--- Running Step 2: Merging Paired-End Reads with BBmerge ---")
    if not os.path.exists(filtered_forward_fasta) or not os.path.exists(filtered_reverse_fasta):
        print("Skipping BBmerge: Filtered FASTA files are missing. Check Step 1 logs.")
    elif not os.path.exists(bbmerge_path):
        print(f"Error: BBmerge script not found at {bbmerge_path}. Please update BBMERGE_PATH. Skipping BBmerge step.")
    else:
        bbmerge_command = [
            bbmerge_path,
            f"in1={filtered_forward_fasta}",
            f"in2={filtered_reverse_fasta}",
            f"out={merged_fasta}",
            f"outunmerged={unmerged_pairs_fasta}",
            "trimq=0",
            "minoverlap=10",
            "mix=f"
        ]
        try:
            print(f"Executing BBmerge in {data_dir}: {' '.join(bbmerge_command)}")
            subprocess.run(bbmerge_command, check=True, cwd=data_dir,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            print("BBmerge completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running BBmerge for {marker_name}: {e}")
            print(f"Command output (stdout):\n{e.stdout}")
            print(f"Command error (stderr):\n{e.stderr}")
        except FileNotFoundError:
            print(f"BBmerge command not found. Ensure {bbmerge_path} is correct and executable.")
        except Exception as e:
            print(f"An unexpected error occurred during BBmerge for {marker_name}: {e}")

    # --- Step 3: Combine Merged ASVs and Unmerged Forwards ---
    print("\n--- Running Step 3: Combining Merged and Unmerged ASVs ---")
    combine_asv_fastas(merged_fasta, unmerged_pairs_fasta, final_combined_fasta)

    # --- Step 4: Create Final ASV Abundance Table ---
    print("\n--- Running Step 4: Creating Final ASV Abundance Table ---")
    create_final_asv_table(asv_counts_file, taxonomy_file, final_combined_fasta, final_output_table,
                           plate_name, marker_name)

    # --- NEW STEP: 5. Generate Histogram and Taxonomy Summary ---
    print("\n--- Running Step 5: Generating ASV Summary (Histogram & Taxonomy) ---")
    generate_asv_summary(
        file_path=final_output_table,
        marker_name=marker_name
    )
    # --- END NEW STEP ---

    print(f"--- Finished processing {plate_name}/{marker_name} ---")


if __name__ == "__main__":
    if not os.path.exists(BBMERGE_PATH):
        print(f"CRITICAL ERROR: BBmerge path not set or invalid: {BBMERGE_PATH}")
        print("Please edit 'run_pipeline_loop.py' and set the correct BBMERGE_PATH.")
        sys.exit(1)

    print(f"Starting pipeline loop from base directory: {BASE_BIO_DIR}")

    plate_folders = [d for d in os.listdir(BASE_BIO_DIR) if os.path.isdir(os.path.join(BASE_BIO_DIR, d)) and d.startswith("Plate")]
    plate_folders.sort()

    if not plate_folders:
        print(f"No 'PlateX' folders found in {BASE_BIO_DIR}. Exiting.")
        sys.exit(0)

    for plate_name in plate_folders:
        plate_path = os.path.join(BASE_BIO_DIR, plate_name)
        print(f"\n===== Entering Plate: {plate_name} =====")

        marker_folders = [d for d in os.listdir(plate_path) if os.path.isdir(os.path.join(plate_path, d))]
        marker_folders.sort()

        if not marker_folders:
            print(f"No marker folders found in {plate_path}. Skipping this plate.")
            continue

        for marker_name in marker_folders:
            data_dir = os.path.join(plate_path, marker_name, "paired")
            run_pipeline_for_dataset(plate_name, marker_name, data_dir, BBMERGE_PATH)

    print("\n--- Pipeline loop completed for all plates and markers ---")
