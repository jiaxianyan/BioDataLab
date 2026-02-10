from Bio import Entrez
import json

def search_pubmed_pmids():
    Entrez.email = "your_email@example.com"
    query = '(transcription factor[Title/Abstract] OR transcription factors[Title/Abstract]) AND (("2025/01/01"[Date - Publication] : "2025/01/05"[Date - Publication])) AND (meta-analysis[ptyp]) AND (hasabstract[text])'

    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=10000)
        record = Entrez.read(handle)
        handle.close()
        
        pmid_list = record["IdList"]
        
        if pmid_list:
            with open("benchmark/gold_results/tf_marker_retrieval_tmp.json", "w") as f:
                json.dump(pmid_list, f, indent=4)
            
    except Exception as e:
        pass

if __name__ == "__main__":
    search_pubmed_pmids()
