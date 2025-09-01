import pandas as pd
from Bio import Entrez
from io import StringIO
import sys

# --- Configuration ---

# ALWAYS tell NCBI who you are. This is a requirement.
# Replace with your actual email address.
Entrez.email = "jiaxianyan@mail.ustc.edu.cn"

# The search query provided in the request.
search_term = "16S V4[All Fields] AND illumina[Platform] AND fecal[All Fields]"

# Number of records to retrieve details for.
# NCBI might limit large requests, so start with a reasonable number.
num_results_to_fetch = 100

# --- Script ---

print(f"Searching SRA for: '{search_term}'")

# Step 1: Use esearch to get the unique IDs (UIDs) of the search results.
try:
    search_handle = Entrez.esearch(
        db="sra",
        term=search_term,
        retmax=num_results_to_fetch  # Get UIDs for the top N results
    )
    search_record = Entrez.read(search_handle)
    search_handle.close()

except Exception as e:
    print(f"An error occurred during the NCBI search: {e}")
    sys.exit(1)


id_list = search_record["IdList"]
total_count = int(search_record["Count"])

print(f"\nFound {total_count} total results.")

if not id_list:
    print("No matching records were found.")
    sys.exit(0)

print(f"Fetching metadata for the top {len(id_list)} results...")

# Step 2: Use efetch to retrieve the full records for the given UIDs.
# We will fetch 'runinfo' which is a convenient CSV format.
try:
    fetch_handle = Entrez.efetch(
        db="sra",
        id=id_list,
        rettype="runinfo",
        retmode="text"
    )
    runinfo_csv = fetch_handle.read()
    fetch_handle.close()

except Exception as e:
    print(f"An error occurred while fetching data from NCBI: {e}")
    sys.exit(1)


# Step 3: Process and display the results using the pandas library.
# Pandas can easily parse the CSV-formatted string into a structured DataFrame.
if runinfo_csv:
    try:
        # Use StringIO to treat the CSV string as a file
        sra_df = pd.read_csv(StringIO(runinfo_csv))

        print("\n--- SRA Run Information (Top Results) ---")
        print(sra_df.head()) # Display the first 5 rows

        # --- Optional: Save the results to a CSV file ---
        output_filename = "sra_search_results.csv"
        sra_df.to_csv(output_filename, index=False)
        print(f"\nSuccessfully saved the full table of {len(sra_df)} results to '{output_filename}'")

    except pd.errors.EmptyDataError:
        print("\nFetched data was empty. No run information to display.")
    except Exception as e:
        print(f"An error occurred while processing the data with pandas: {e}")

else:
    print("\nCould not retrieve run information for the found IDs.")