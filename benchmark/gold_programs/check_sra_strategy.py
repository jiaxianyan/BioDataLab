from Bio import Entrez
import pandas as pd
import sys
import io

# Entrez.email = "your_email@example.com"

def fetch_runinfo(run_id):
    """Fetch SRA runinfo table for a given SRR/ERR accession"""
    handle = Entrez.esearch(db="sra", term=run_id)
    record = Entrez.read(handle)
    handle.close()

    if len(record["IdList"]) == 0:
        print("No record found!")
        return None

    sra_id = record["IdList"][0]

    handle = Entrez.efetch(db="sra", id=sra_id, rettype="runinfo", retmode="text")
    runinfo = handle.read()
    handle.close()
    # print(runinfo)
    df = pd.read_csv(io.StringIO(runinfo.decode('utf-8')))
    print(df)
    return df


def classify_library(strategy):
    """Classify sequencing type based on LibraryStrategy"""
    if strategy == "AMPLICON":
        return "16S amplicon"
    elif strategy == "METAGENOMIC":
        return "Shotgun metagenome"
    elif strategy == "METATRANSCRIPTOMIC":
        return "Metatranscriptome (metaRNA)"
    elif strategy == "RNA-Seq":
        return "Host RNA-seq"
    else:
        return f"Else strategy: {strategy}"


if __name__ == "__main__":
    run_id = sys.argv[1]

    df = fetch_runinfo(run_id)
    if df is None:
        sys.exit()

    strategy = df.loc[0, "LibraryStrategy"]
    layout = df.loc[0, "LibraryLayout"]
    save_data_path = "~/biodatalab/benchmark/dataset/mBodyMap/metadata_info/"
    df.to_csv(save_data_path + run_id +"_info.csv")
    print("===================================")
    print("Run accession:", run_id)
    print("LibraryStrategy:", strategy)
    print("LibraryLayout:", layout)
    print("Classification:", classify_library(strategy))
    print("===================================")
