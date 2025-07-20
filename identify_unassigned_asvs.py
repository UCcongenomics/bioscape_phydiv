#/projects/RachelMeyer/BioSCAPE/PlateA_ITS2_Plants/paired/clve5i5mn0003l50gam95h4jx-ITS2_Plants-paired.txt

# identify_unassigned_asvs.py
import csv

def identify_unassigned_asvs(stats_file_path, output_asv_list_path):
    """
    Identifies 'unassigned' ASVs from a DADA2 statistics file based on taxonomic path.

    Args:
        stats_file_path (str): Path to the ASV statistics file (e.g., ITS2_Plants-paired.txt).
        output_asv_list_path (str): Path to save the list of unassigned ASV names.
    """
    unassigned_asv_names = []
    unassigned_keywords = ["unassigned", "unclassified", "no_hit", "root"] # Common keywords for unassigned
    single_domain_paths = ["eukaryota", "bacteria", "archaea", "viruses"] # Common single domain paths

    try:
        with open(stats_file_path, mode='r', newline='') as infile:
            reader = csv.DictReader(infile, delimiter='\t')
            if 'readname' not in reader.fieldnames or 'taxonomic_path' not in reader.fieldnames:
                print("Error: Input file must contain 'readname' and 'taxonomic_path' columns.")
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
                if not is_unassigned:
                    path_components = tax_path.split(';')
                    if len(path_components) == 1 and path_components[0] in single_domain_paths:
                        is_unassigned = True

                if is_unassigned:
                    unassigned_asv_names.append(asv_name)

        with open(output_asv_list_path, mode='w') as outfile:
            for asv_name in unassigned_asv_names:
                outfile.write(f"{asv_name}\n")

        print(f"Identified {len(unassigned_asv_names)} unassigned ASVs.")
        print(f"List of unassigned ASV names saved to: {output_asv_list_path}")

    except FileNotFoundError:
        print(f"Error: The file '{stats_file_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # --- Configuration ---
    # !!! IMPORTANT: Replace this with the actual path to your statistics file !!!
    STATS_FILE = "/projects/RachelMeyer/BioSCAPE/PlateA_ITS2_Plants/paired/clve5i5mn0003l50gam95h4jx-ITS2_Plants-paired.txt"
    OUTPUT_ASV_LIST = "unassigned_asvs_to_remove.txt"
    # ---------------------

    identify_unassigned_asvs(STATS_FILE, OUTPUT_ASV_LIST)
