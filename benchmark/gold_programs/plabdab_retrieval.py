import time
from pathlib import Path
from Bio import Entrez, SeqIO

Entrez.email = "your_email@example.com"
# Entrez.api_key = "your_ncbi_api_key" 
# ['antibody', 'immunoglobulin', 'IgG', 'Fab', 'scFv']
# nanobody_keywords = [
#     'nanobody', 'single domain', 'VHH', 'vhh',
#     'camelid', 'llama', 'alpaca', 'shark', 'camel',
#     'sdAb', 'heavy-chain only', 'single-domain'
# ]

def PLAbDab_fetch_and_filter(
    out_fasta: str = "stage1_filtered.fasta", 
    max_rets: int = 10000,
    batch_size: int = 1000,
):
    query = ('(antibody OR antibodies OR immunoglobulin OR scfv OR bcr) '
                'NOT (nanobody OR nanobodies)')
    print("Start NCBI search...")
    search_handle = Entrez.esearch(db="protein", term=query, usehistory="y")
    search_results = Entrez.read(search_handle)
    count = min(max_rets, int(search_results["Count"]))
    
    Path(out_fasta).parent.mkdir(exist_ok=True, parents=True)
    with open(out_fasta, "w") as out_f:
        for start in range(0, count, batch_size):
            fetch_handle = Entrez.efetch(
                db="protein", rettype="fasta", retmode="text",
                retstart=start, retmax=batch_size,
                webenv=search_results["WebEnv"], 
                query_key=search_results["QueryKey"],
            )
            
            # Streaming the current batch
            for record in SeqIO.parse(fetch_handle, "fasta"):
                seq_len = len(record.seq)
                # sequence length filtering
                if 70 <= seq_len <= 1000:
                    SeqIO.write(record, out_f, "fasta")
            
            print(f"Progress: {start}/{count}", end="\r")
            time.sleep(0.1)

if __name__ == "__main__":
    PLAbDab_fetch_and_filter(out_fasta = "benchmark/gold_results/antibody_seq_retrieval.fasta")