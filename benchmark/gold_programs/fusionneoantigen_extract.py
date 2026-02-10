from Bio import Entrez
import json

Entrez.email = "your_email@example.com"
Entrez.tool = "Neoantigen_Search_Bot"

def search_bcr_abl_neoantigens():
    keywords = (
        "("
        '"BCR-ABL1"[All Fields] OR "BCR-ABL"[All Fields] OR "BCR/ABL"[All Fields]'
        ") AND ("
        '"neoantigen"[All Fields] OR "neo-antigen"[All Fields] OR '
        '"neoepitope"[All Fields] OR "tumor specific antigen"[All Fields]'
        ")"
    )
    
    date_range = '"2010/01/01"[Date - Publication] : "2025/01/01"[Date - Publication]'
    final_query = f"{keywords} AND {date_range}"
    
    try:
        database = "pmc"
        handle = Entrez.esearch(
            db=database, 
            term=final_query, 
            retmax=1000,
            retmode="xml"
        )
        record = Entrez.read(handle)
        handle.close()

        pmid_list = record["IdList"]
        
        if pmid_list:
            with open("benchmark/gold_results/fusionneoantigen_extract_tmp.json", "w") as f:
                json.dump(pmid_list, f, indent=4)

    except Exception as e:
        return []


if __name__ == "__main__":
    search_bcr_abl_neoantigens()

