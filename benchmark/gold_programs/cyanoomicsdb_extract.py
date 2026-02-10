import pandas as pd
from Bio import Entrez

def fetch_strict_cyanobacteria_papers():
    Entrez.email = "your_email@example.com"
    search_query = 'cyanobacteria AND ("1900/01/01"[Date - Publication] : "2025/12/31"[Date - Publication])'
    try:
        search_handle = Entrez.esearch(
            db="pubmed",
            term=search_query,
            sort="pub_date",
            retmax=5
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()
        id_list = search_results["IdList"]
        if not id_list:
            return
        fetch_handle = Entrez.efetch(
            db="pubmed",
            id=id_list,
            retmode="xml"
        )
        papers = Entrez.read(fetch_handle)
        fetch_handle.close()
        data = []
        for article in papers['PubmedArticle']:
            try:
                citation = article['MedlineCitation']
                article_data = citation['Article']
                title = article_data.get('ArticleTitle', 'N/A')
                journal = article_data.get('Journal', {}).get('Title', 'N/A')
                doi = "N/A"
                if 'ELocationID' in article_data:
                    for eid in article_data['ELocationID']:
                        if getattr(eid, 'attributes', {}).get('EIdType') == 'doi':
                            doi = str(eid)
                            break
                if doi == "N/A" and 'PubmedData' in article:
                    article_ids = article['PubmedData'].get('ArticleIdList', [])
                    for aid in article_ids:
                        if getattr(aid, 'attributes', {}).get('IdType') == 'doi':
                            doi = str(aid)
                            break
                pub_date_data = article_data.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
                pub_year = pub_date_data.get('Year', '')
                pub_month = pub_date_data.get('Month', '')
                pub_date_str = f"{pub_year} {pub_month}".strip()
                data.append({
                    "Title": title,
                    "Journal": journal,
                    "DOI": doi,
                    "PubDate": pub_date_str
                })
            except Exception as e:
                continue
        df = pd.DataFrame(data)
        output_df = df[["Title", "Journal", "DOI"]]
        filename = "benchmark/gold_results/cyanoomicsdb_extract.csv"
        output_df.to_csv(filename, index=False, encoding='utf-8-sig')
    except Exception as e:
        pass

if __name__ == "__main__":
    fetch_strict_cyanobacteria_papers()