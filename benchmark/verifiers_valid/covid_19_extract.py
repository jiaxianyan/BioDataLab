import os

def covid_19_extract(d):
    return os.path.exists(f'{d}/covid_19_extract.json')