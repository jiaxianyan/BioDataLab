from Bio import Entrez
import json

Entrez.email = "your_email@example.com"
Entrez.tool = "DDI_Search_Script"

def search_ciprofloxacin_ddi():
    query = (
        "("
        '"Ciprofloxacin"[MeSH Terms] OR "Ciprofloxacin"[Title/Abstract]' 
        ") AND ("
        '"Drug Interactions"[MeSH Terms] OR "drug-drug interaction"[Title/Abstract] OR "DDI"[Title/Abstract] OR "drug interaction"[Title/Abstract]'
        ") AND ("
        '"2021/01/01"[Date - Publication] : "2021/03/01"[Date - Publication]'
        ") NOT ("
        '"Food-Drug Interactions"[MeSH Terms] OR "food"[Title] OR '
        '"Pharmacogenetics"[MeSH Terms]  OR "gene"[Title/Abstract]'
        ")"
    )
    
    try:
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=1000,
            retmode="xml"
        )
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]
        return id_list
    except Exception:
        return []

if __name__ == "__main__":
    pmid_list = search_ciprofloxacin_ddi()
    with open("benchmark/gold_results/ddinter_retrieval_tmp.json", "w") as f:
        json.dump(pmid_list, f, indent=4)