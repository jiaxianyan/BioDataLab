import json
from Bio import Entrez

def save_hsc_pmids_to_json():
    Entrez.email = "your.email@example.com"
    
    query = (
        '("Hematopoietic Stem Cells"[Title/Abstract] OR "HSC"[Title/Abstract]) '
        'AND ("single-cell RNA sequencing"[Title/Abstract] OR "scRNA-seq"[Title/Abstract]) '
        'NOT ("Hepatic Stellate Cells")'
    )
    
    try:
        handle = Entrez.esearch(
            db="pubmed", 
            term=query, 
            mindate="2023/02/10", 
            maxdate="2023/03/10", 
            datetype="pdat", 
            retmax=1000
        )
        
        results = Entrez.read(handle)
        handle.close()
        
        pmid_list = results["IdList"]
        
        output_filename = "benchmark/gold_results/stemdriver_retrieval_tmp.json"
        with open(output_filename, 'w', encoding='utf-8') as f:
            json.dump(pmid_list, f, indent=4)
            
    except Exception as e:
        pass

if __name__ == "__main__":
    save_hsc_pmids_to_json()