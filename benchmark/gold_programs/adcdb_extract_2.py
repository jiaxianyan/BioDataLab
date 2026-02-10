import time
import json
import pandas as pd
from Bio import Entrez

Entrez.email = "your_email@example.com"

SEARCH_KEYWORDS = ['Zilovertamab vedotin', 'MK-2140', 'VLS-101']
BIOLOGICAL_ACTIVITY_KEYWORDS = [
    'biological activity', 'pharmacology', 'mechanism of action',
    'efficacy', 'cytotoxicity', 'pharmacokinetics', 'pharmacodynamics',
    'biodistribution', 'internalization', 'potency', 'ADC'
]

SPECIES_QUERY = ''
DATE_RANGE = ''

OUTPUT_FILE = "benchmark/gold_results/adcdb_extract_2_tmp.json"

def build_search_query():
    drug_or = " OR ".join([f'"{kw}"[Title/Abstract]' for kw in SEARCH_KEYWORDS])
    activity_or = " OR ".join([f'"{kw}"[Title/Abstract]' for kw in BIOLOGICAL_ACTIVITY_KEYWORDS])
    
    final_query =  f"({drug_or}) AND ({activity_or})"
    
    if DATE_RANGE:
        final_query += f" AND {DATE_RANGE}"
    if SPECIES_QUERY:
        final_query += f" AND {SPECIES_QUERY}"
        
    return final_query

def search_pubmed(query):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=100000)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def verify_and_collect(id_list):
    batch_size = 200
    final_pmids = []

    for i in range(0, len(id_list), batch_size):
        batch_ids = id_list[i:i+batch_size]
        try:
            handle = Entrez.efetch(db="pubmed", id=",".join(batch_ids), rettype="xml", retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            for article in records['PubmedArticle']:
                medline = article['MedlineCitation']
                article_data = medline['Article']
                pmid = str(medline['PMID'])

                languages = [l.lower() for l in article_data.get('Language', [])]
                if 'eng' not in languages:
                    continue

                final_pmids.append(pmid)

            time.sleep(0.3)
        except Exception as e:
            continue

    return final_pmids

if __name__ == "__main__":
    pubmed_query = build_search_query()
    
    ids = search_pubmed(pubmed_query)
    
    if ids:
        results = verify_and_collect(ids)
        
        with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=4)