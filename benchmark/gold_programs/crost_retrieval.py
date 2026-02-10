import time
import json
from Bio import Entrez

Entrez.email = "your_email@example.com"

SPECIES = '"Homo sapiens"[Mesh] OR "human"[Title/Abstract]'
TOPIC = '"spatial transcriptomics"[Title/Abstract] OR "spatially resolved transcriptomics"[Title/Abstract]'
DATE_RANGE = '"2021/01/01"[PDAT] : "2021/02/01"[PDAT]'
EXCLUDE_TYPES = 'NOT ("review"[PT] OR "editorial"[PT] OR "comment"[PT] OR "letter"[PT])'
OUTPUT_FILE = "benchmark/gold_results/crost_retrieval_tmp.json"

def build_query():
    query = f"({SPECIES}) AND ({TOPIC}) AND ({DATE_RANGE}) {EXCLUDE_TYPES}"
    return query

def search_ncbi(query):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=100000)
    record = Entrez.read(handle)
    handle.close()
    id_list = record["IdList"]
    return id_list

def save_results(pmid_list, file_path):
    with open(file_path, 'w', encoding='utf-8') as f:
        json.dump(pmid_list, f, indent=4)

if __name__ == "__main__":
    final_query = build_query()
    results_pmids = search_ncbi(final_query)
    if results_pmids:
        save_results(results_pmids, OUTPUT_FILE)