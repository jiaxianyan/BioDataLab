import json
from Bio import Entrez

Entrez.email = "your_email@example.com"

SPECIES_QUERY = '"Homo sapiens"[Mesh] OR "human"[Title/Abstract]'
METHOD_QUERY = (
    '"single-cell RNA sequencing"[Title/Abstract] OR "scRNA-seq"[Title/Abstract] OR '
    '"single-nucleus RNA sequencing"[Title/Abstract] OR "snRNA-seq"[Title/Abstract]'
)
TISSUE_QUERY = (
    '"brain"[Title/Abstract] OR "spinal cord"[Title/Abstract] OR '
    '"retina"[Title/Abstract] OR "embryo"[Title/Abstract]'
)
DATA_AVAILABILITY_QUERY = (
    '"dataset"[Title/Abstract] OR "datasets"[Title/Abstract] OR '
    '"GSE"[All Fields] OR "accession number"[All Fields] OR "publicly available"[Title/Abstract]'
)
DATE_QUERY = '"2021/01/01"[PDAT] : "2021/04/01"[PDAT]'
OUTPUT_FILE = "benchmark/gold_results/scan_retrieval_tmp.json"

def build_query():
    query = (
        f"({SPECIES_QUERY}) AND "
        f"({METHOD_QUERY}) AND "
        f"({TISSUE_QUERY}) AND "
        f"({DATA_AVAILABILITY_QUERY}) AND "
        f"({DATE_QUERY})"
    )
    return query

def search_pubmed(query):
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=10000)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        return []

if __name__ == "__main__":
    final_query = build_query()
    pmid_list = search_pubmed(final_query)
    if pmid_list:
        with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
            json.dump(pmid_list, f, indent=4)