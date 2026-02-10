from Bio import Entrez
import json

def search_pubmed_pmids():
    Entrez.email = "your_email@example.com"
    keywords = ('("Curcumin"[Mesh] OR Curcumin[Title/Abstract]) AND '
                '("Drug Therapy, Combination"[Mesh] OR "Drug Synergism"[Mesh] OR '
                'combin*[Title/Abstract] OR synerg*[Title/Abstract])')
    
    date_range = '"2025/11/01"[PDAT] : "2025/11/05"[PDAT]'
    article_type = '"Review"[Publication Type]'
    full_query = f"{keywords} AND {date_range} AND {article_type}"

    try:
        handle = Entrez.esearch(db="pubmed", term=full_query, retmax=10000)
        record = Entrez.read(handle)
        handle.close()
        
        pmid_list = record["IdList"]
        
        if pmid_list:
            with open("benchmark/gold_results/npcdr_retrieval_tmp.json", "w") as f:
                json.dump(pmid_list, f, indent=4)
            
    except Exception as e:
        pass

if __name__ == "__main__":
    search_pubmed_pmids()
