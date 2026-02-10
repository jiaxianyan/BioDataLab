import json
import requests
import time
from Bio import Entrez
from datetime import datetime

Entrez.email = "your_email@example.com"

DRUG = "Docetaxel"
SPECIES = '"Homo sapiens"[Organism]'
CANCER = "cancer OR tumor OR carcinoma OR neoplasm"
START_DATE = "2021/01/01"
END_DATE = "2021/06/01"

def search_geo():
    query = f"({DRUG}) AND ({SPECIES}) AND ({CANCER}) AND \"gse\"[Filter] AND (\"{START_DATE}\"[PDAT] : \"{END_DATE}\"[PDAT])"
    handle = Entrez.esearch(db="gds", term=query, retmax=100)
    record = Entrez.read(handle)
    handle.close()
    id_list = record["IdList"]
    gse_list = []
    if id_list:
        handle = Entrez.esummary(db="gds", id=",".join(id_list))
        summary_records = Entrez.read(handle)
        handle.close()
        for record in summary_records:
            if 'Accession' in record:
                accession = record['Accession']
                if accession.startswith('GSE'):
                    gse_list.append(accession)
    return gse_list

if __name__ == "__main__":
    gse_ids = search_geo()
    output_file = "benchmark/gold_results/ctr_db_retrieval_tmp.json"
    with open(output_file, 'w') as f:
        json.dump(gse_ids, f, indent=4)