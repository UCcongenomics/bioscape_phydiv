# identify_unassigned_asvs.py
import csv
import os # Added for path checking

def identify_unassigned_asvs(stats_file_path, output_asv_list_path):
    """
    Identifies 'unassigned' ASVs from a DADA2 statistics file based on taxonomic path
    and writes their names to a list file.

    Args:
        stats_file_path (str): Path to the ASV statistics file (e.g., ITS2_Plants-paired.txt).
        output_asv_list_path (str): Path to save the list of unassigned ASV names.
    """
    unassigned_asv_names = []
    # Common keywords for unassigned or poorly classified taxonomy
    unassigned_keywords = ["unassigned", "unclassified", "no_hit", "root", "environmental_sample"]
    # Common single domain paths that might indicate broad or unspecific classification
    single_domain_paths = ["eukaryota", "bacteria", "archaea", "viruses"]

    print(f"Identifying unassigned ASVs from: {stats_file_path}")

    try:
        if not os.path.exists(stats_file_path):
            print(f"Error: Input stats file '{stats_file_path}' not found. Cannot identify unassigned ASVs.")
            return

        with open(stats_file_path, mode='r', newline='') as infile:
            reader = csv.DictReader(infile, delimiter='\t')
            if 'readname' not in reader.fieldnames or 'taxonomic_path' not in reader.fieldnames:
                print("Error: Input file must contain 'readname' and 'taxonomic_path' columns. Skipping.")
                return

            for row in reader:
                asv_name = row['readname'].strip()
                tax_path = row['taxonomic_path'].strip().lower() # Convert to lowercase for case-insensitive check

                is_unassigned = False

                # Check for explicit unassigned keywords
                for keyword in unassigned_keywords:
                    if keyword in tax_path:
                        is_unassigned = True
                        break

                # Check for single domain assignment if not already marked unassigned
                # This logic assumes single component paths are undesirable (e.g., just "bacteria")
                if not is_unassigned:
                    path_components = [pc.strip() for pc in tax_path.split(';') if pc.strip()] # Split and clean
                    if len(path_components) == 1 and path_components[0] in single_domain_paths:
                        is_unassigned = True

                if is_unassigned:
                    unassigned_asv_names.append(asv_name)

        # Ensure the output directory exists
        output_dir = os.path.dirname(output_asv_list_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        with open(output_asv_list_path, mode='w') as outfile:
            for asv_name in unassigned_asv_names:
                outfile.write(f"{asv_name}\n")

        print(f"Identified {len(unassigned_asv_names)} unassigned ASVs. List saved to: {output_asv_list_path}")

    except FileNotFoundError:
        print(f"Error: The file '{stats_file_path}' was not found. Cannot identify unassigned ASVs.")
    except Exception as e:
        print(f"An unexpected error occurred during identification of unassigned ASVs: {e}")

# (The __main__ block is now commented out or can be deleted)
# if __name__ == "__main__":
#     # --- Configuration ---
#     STATS_FILE = "clve5i5mn0003l50gam95h4jx-ITS2_Plants-paired.txt"
#     OUTPUT_ASV_LIST = "unassigned_asvs_to_remove.txt"
#     identify_unassigned_asvs(STATS_FILE, OUTPUT_ASV_LIST)
