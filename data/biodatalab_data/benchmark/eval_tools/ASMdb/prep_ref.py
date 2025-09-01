# from Bio import Entrez
# import json
# # Entrez.email = "jiaxianyan@mail.ustc.edu.cn"

# Entrez.email = 'youmail@mail.com'
# # 1. Searching the NCBI GEO database

# def search_geo(search_term, end_date="2019/12/31", retmax=10000):
#     """Action: Search the GEO database for relevant studies."""
#     print(f"Searching GEO with term: '{search_term}' until {end_date}")
    
#     handle = Entrez.esearch(db="gds", term=search_term, retmax=retmax)
#     record = Entrez.read(handle)
#     handle.close()
    
#     gse_ids = record["IdList"]
#     print(f"Found {len(gse_ids)} potential GEO database ")
#     print(record.keys())
#     return gse_ids


    
# SEARCH_TERM = "whole-genome bisulfite sequencing[All Fields]) AND (\"0001/01/01\"[PDAT] : \"2019/10/31\"[PDAT])"
# gse_ids = search_geo(SEARCH_TERM)

# try:
#     handle = Entrez.esummary(db="gds", id=",".join(gse_ids))
#     summaries = Entrez.read(handle)
#     handle.close()
# except Exception as e:
#     print(f"An error occurred while fetching summaries: {e}")

# Accession_list = []
# for i, summary in enumerate(summaries):
#     # print(summary)
#     gse_accession = summary["Accession"]
#     # print(gse_accession)
#     # title = summary["title"]
#     # input(title)
#     Accession_list.append(gse_accession)

# with open('data/biodatalab_data/benchmark/results/ASMdb/bs_seq_2019_accessions_ref.json', 'w') as f:
#     json.dump(Accession_list, f)


# import json
# import os
# from pathlib import Path
# import GEOparse

# def get_gsm_sample_ids_from_gse(gse_id, output_file):
#     """
#     Fetch all GSM sample IDs from a given GSE series and save to JSON file.
    
#     Args:
#         gse_id (str): The GSE series ID (e.g., 'GSE285258')
#         output_file (str): Path to output JSON file
#     """
#     try:
#         print(f"Fetching data for {gse_id}...")
        
#         # Download and parse the GSE data
#         gse = GEOparse.get_GEO(geo=gse_id, destdir="./temp_geo")
        
#         # Extract GSM sample IDs
#         gsm_ids = list(gse.gsms.keys())
        
#         print(f"Found {len(gsm_ids)} GSM samples in {gse_id}")
        
#         # Create output directory if it doesn't exist
#         output_path = Path(output_file)
#         output_path.parent.mkdir(parents=True, exist_ok=True)
        
#         # Save GSM IDs to JSON file
#         with open(output_file, 'w') as f:
#             json.dump({
#                 "gse_id": gse_id,
#                 "gsm_count": len(gsm_ids),
#                 "gsm_sample_ids": gsm_ids
#             }, f, indent=2)
        
#         print(f"GSM sample IDs saved to: {output_file}")
#         return gsm_ids
        
#     except Exception as e:
#         print(f"Error fetching data for {gse_id}: {str(e)}")
#         return None

# # Main execution
# if __name__ == "__main__":
#     gse_id = "GSE285258"
#     output_file = "workdir/ncbi_gsm_sample_ids.json"
    
#     # Get GSM sample IDs and save to file
#     gsm_ids = get_gsm_sample_ids_from_gse(gse_id, output_file)
    
#     if gsm_ids:
#         print(f"Successfully retrieved {len(gsm_ids)} GSM sample IDs")
#         print("First 5 GSM IDs:")
#         for gsm_id in gsm_ids[:5]:
#             print(f"  - {gsm_id}")
#     else:
#         print("Failed to retrieve GSM sample IDs")

# import GEOparse
# from Bio import Entrez
# import os
# import json
# import re

# # --- Configuration ---
# GSM_ID = "GSM864033"
# WORKDIR = "workdir"
# OUTPUT_FILE = os.path.join(WORKDIR, "ncbi_sra_run_ids.json")

# # IMPORTANT: Always tell NCBI who you are
# # Replace with your actual email address
# Entrez.email = "your.email@example.com"


# def get_sra_project_from_gsm(gsm):
#     """Extracts the SRA project ID (SRP, ERP, DRP) from GSM metadata."""
#     if 'relation' not in gsm.metadata:
#         print(f"No 'relation' field found in metadata for {gsm.name}.")
#         return None

#     for relation in gsm.metadata['relation']:
#         if relation.startswith("SRA:"):
#             # Use regex to find the SRA Project/Experiment ID
#             match = re.search(r'(SRP|ERP|DRP|SRX)\d+', relation)
#             if match:
#                 sra_id = match.group(0)
#                 print(f"Found SRA ID: {sra_id}")
#                 return sra_id
    
#     print(f"No SRA link found in the 'relation' field for {gsm.name}.")
#     return None


# def get_sra_runs_from_project(sra_id):
#     """Fetches all SRA Run IDs (SRR) for a given SRA Project/Experiment ID."""
#     print(f"Querying NCBI SRA for runs associated with {sra_id}...")
#     try:
#         # Use esearch to find links from the project to run objects
#         handle = Entrez.esearch(db="sra", term=sra_id, retmax="1000")
#         record = Entrez.read(handle)
#         handle.close()
        
#         if not record['IdList']:
#             print(f"No SRA records found for term {sra_id}.")
#             return []

#         # Use esummary to get details for the found IDs
#         id_list = record['IdList']
#         handle = Entrez.esummary(db="sra", id=",".join(id_list))
#         summary_records = Entrez.read(handle)
#         handle.close()
        
#         run_ids = []
#         # The structure of the summary can be complex. We need to parse it carefully.
#         # Each record in summary_records corresponds to an SRA Experiment (SRX)
#         for record in summary_records:
#             # The Run information is nested inside an XML string called 'Runs'
#             # We can parse the run accession from the XML attributes
#             run_xml = record.get('Runs', '')
#             # Find all occurrences of acc="SRR..." in the XML string
#             runs_in_record = re.findall(r'acc="(SRR\d+)"', run_xml)
#             if runs_in_record:
#                 run_ids.extend(runs_in_record)

#         print(f"Found {len(run_ids)} SRA Run(s): {run_ids}")
#         return run_ids

#     except Exception as e:
#         print(f"An error occurred while querying NCBI: {e}")
#         return []

# def main():
#     """Main function to execute the workflow."""
#     # 1. Create the working directory if it doesn't exist
#     os.makedirs(WORKDIR, exist_ok=True)
#     print(f"Working directory '{WORKDIR}' is ready.")

#     # 2. Get GEO sample data
#     print(f"Fetching GEO metadata for {GSM_ID}...")
#     try:
#         gsm = GEOparse.get_GEO(geo=GSM_ID, silent=True)
#     except Exception as e:
#         print(f"Could not download data for {GSM_ID}. Error: {e}")
#         return

#     # 3. Extract SRA project ID from the metadata
#     sra_project_id = get_sra_project_from_gsm(gsm)

#     if not sra_project_id:
#         print("Could not proceed without an SRA Project ID.")
#         return
        
#     # 4. Get all SRA Run IDs from the SRA Project ID
#     sra_run_ids = get_sra_runs_from_project(sra_project_id)

#     # 5. Save the list of SRA Run IDs to a JSON file
#     if sra_run_ids:
#         print(f"Saving {len(sra_run_ids)} SRA Run ID(s) to {OUTPUT_FILE}...")
#         try:
#             with open(OUTPUT_FILE, 'w') as f:
#                 json.dump(sra_run_ids, f, indent=4)
#             print("Successfully saved the file.")
#             print("\n--- Content of the file ---")
#             print(json.dumps(sra_run_ids, indent=4))
#             print("---------------------------")
#         except IOError as e:
#             print(f"Error writing to file {OUTPUT_FILE}. Error: {e}")
#     else:
#         print("No SRA Run IDs were found to save.")

# if __name__ == "__main__":
#     main()


import subprocess
import os

# Input files
input_files = [
    'data/biodatalab_data/benchmark/results/ASMdb/raw_data/SRR7207779_1.fastq',
    'data/biodatalab_data/benchmark/results/ASMdb/raw_data/SRR7207779_2.fastq'
]

# Assuming paired-end reads: first is read1, second is read2
input1 = input_files[0]
input2 = input_files[1]

# Derive output file names
output_dir = 'trimmed'  # Create a trimmed directory for outputs
os.makedirs(output_dir, exist_ok=True)

base1 = os.path.basename(input1)
base2 = os.path.basename(input2)

output1 = os.path.join(output_dir, base1.replace('.fastq', '_trimmed.fastq'))
output2 = os.path.join(output_dir, base2.replace('.fastq', '_trimmed.fastq'))

# Fastp command with specified parameters
command = [
    'fastp',
    '-i', input1,
    '-I', input2,
    '-o', output1,
    '-O', output2,
    '-W', '4',
    '-M', '20',
    '-q', '15',
    '-u', '40',
    '-n', '5',
    '-Y', '0'
]

# Run the command
try:
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    print("Fastp trimming completed successfully.")
    print("Output files:", output1, output2)
    print("Standard output:", result.stdout)
    print("Standard error:", result.stderr)
except subprocess.CalledProcessError as e:
    print("Error running fastp:", e)
    print("Standard output:", e.stdout)
    print("Standard error:", e.stderr)