from Bio import Entrez
import json

def search_geo_datasets():
    Entrez.email = "your_email@example.com"  
    
    query = (
    # 癌症部分 - 使用 MeSH 术语提高查全率
    '("Neoplasms"[MeSH Terms] OR "cancer"[Title/Abstract] OR "carcinoma"[Title/Abstract] OR '
    '"tumor"[Title/Abstract] OR "neoplasm"[Title/Abstract] OR "malignancy"[Title/Abstract]) AND '
    
    # 药物/治疗部分
    '("Drug Therapy"[MeSH Terms] OR "treatment"[Title/Abstract] OR "therapy"[Title/Abstract] OR '
    '"drug"[Title/Abstract] OR "response"[Title/Abstract] OR "resistance"[Title/Abstract] OR '
    '"sensitivity"[Title/Abstract]) AND '
    
    # 患者/临床部分
    '("Patients"[MeSH Terms] OR "clinical"[Title/Abstract] OR "patient"[Title/Abstract] OR '
    '"man"[Title/Abstract] OR "woman"[Title/Abstract] OR "sample"[Title/Abstract]) AND '
    
    # 物种限制
    '"Homo sapiens"[Organism]'
)
    date_range = '"2019/10/01"[PDAT] : "2019/10/15"[PDAT]'
    entry_type = '"GSE"[Filter]' 
    
    full_query = f"{query} AND {date_range} AND {entry_type}"
    
    try:
        handle = Entrez.esearch(
            db="gds",           
            term=full_query,
            retmax=1000,          
            retmode="xml"
        )
        record = Entrez.read(handle)
        handle.close()
        
        pmid_list = record["IdList"]
        
        if pmid_list:
            with open("benchmark/gold_results/cds_db_retrieval_tmp.json", "w") as f:
                json.dump(pmid_list, f, indent=4)
            
    except Exception as e:
        pass

if __name__ == "__main__":
    search_geo_datasets()