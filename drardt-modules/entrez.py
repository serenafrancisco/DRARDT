import requests
from Bio import Entrez
import json

def get_publication_count(gene_name):
    search_term = f"{gene_name} [Title/Abstract] AND 2000:2024 [PDat]"
    handle = Entrez.esearch(db="pubmed", term=search_term, retmax=1)
    record = Entrez.read(handle)
    return int(record["Count"])

def assign_publication_index(pub_count):
    if pub_count > 250:
        return "High"
    elif 100 <= pub_count <= 250:
        return "Medium"
    else:
        return "Low"



