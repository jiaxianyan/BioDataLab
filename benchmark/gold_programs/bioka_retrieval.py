import time
import re
import nltk
import pandas as pd
import json
from Bio import Entrez
from tqdm import tqdm

Entrez.email = "your_email@example.com"
SEARCH_KEYWORDS = ['biomarker', 'marker', 'indicator', 'predictor']
SPECIES_QUERY = '("cats"[MeSH Terms] OR "feline"[Title/Abstract] OR "cat"[Title/Abstract])'
DATE_RANGE = '"2022/01/01"[PDAT] : "2022/01/05"[PDAT]'
ENTITY_TERMS = r"biomarker|marker|indicator|target"
QUALIFYING_WORDS = r"diagnostic|prognostic|therapeutic|valuable"
EXCLUDE_TYPES = ['review', 'comment', 'letter', 'editorial', 'meta-analysis']

def download_nltk_data():
    try:
        nltk.data.find('tokenizers/punkt')
    except LookupError:
        nltk.download('punkt')
        nltk.download('punkt_tab')

def build_search_query():
    keywords_or = " OR ".join([f'"{kw}"[Title/Abstract]' for kw in SEARCH_KEYWORDS])
    final_query = f"({keywords_or}) AND {SPECIES_QUERY} AND {DATE_RANGE}"
    return final_query

def search_pubmed(query):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=100000)
    record = Entrez.read(handle)
    handle.close()
    id_list = record["IdList"]
    return id_list

def fetch_details_and_filter(id_list):
    batch_size = 200
    filtered_articles = []
    re_entity = re.compile(ENTITY_TERMS, re.IGNORECASE)
    re_qualifier = re.compile(QUALIFYING_WORDS, re.IGNORECASE)

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
                title = article_data.get('ArticleTitle', '')
                journal = article_data.get('Journal', {}).get('Title', '')
                languages = article_data.get('Language', [])
                if 'eng' not in languages:
                    continue
                pub_types = article_data.get('PublicationTypeList', [])
                pub_types_str = [str(pt).lower() for pt in pub_types]
                if any(exclude in pt for pt in pub_types_str for exclude in EXCLUDE_TYPES):
                    continue
                abstract_text = ""
                if 'Abstract' in article_data and 'AbstractText' in article_data['Abstract']:
                    abs_list = article_data['Abstract']['AbstractText']
                    if isinstance(abs_list, list):
                        abstract_text = " ".join([str(t) for t in abs_list])
                    else:
                        abstract_text = str(abs_list)
                if not abstract_text:
                    continue
                sentences = nltk.sent_tokenize(abstract_text)
                match_found = False
                matched_sentences = []
                for sentence in sentences:
                    has_entity = re_entity.search(sentence)
                    has_qualifier = re_qualifier.search(sentence)
                    if has_entity and has_qualifier:
                        match_found = True
                        matched_sentences.append(sentence)
                if match_found:
                    filtered_articles.append({
                        'PMID': pmid,
                        'Title': title,
                        'Journal': journal,
                        'Type': ", ".join(pub_types),
                        'Abstract': abstract_text,
                        'Matched_Sentences': " || ".join(matched_sentences)
                    })
            time.sleep(0.5)
        except Exception as e:
            continue
    return pd.DataFrame(filtered_articles)

if __name__ == "__main__":
    download_nltk_data()
    query = build_search_query()
    ids = search_pubmed(query)
    if len(ids) > 0:
        df_results = fetch_details_and_filter(ids)
        pmid_list = df_results['PMID'].tolist()
        with open('benchmark/gold_results/bioka_retrieval_tmp.json', 'w') as f:
            json.dump(pmid_list, f, indent=4)