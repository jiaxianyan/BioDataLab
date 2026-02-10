import os

def cancermirnome_annotate(d):
    return os.path.exists(f'{d}/cancermirnome_annotate.txt')