import xml.etree.ElementTree as ET
import requests
import json

def get_human_uniprot_id(gene_name):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_name}+AND+organism_id:9606&fields=accession"
    response = requests.get(url)
    data = response.json()
    return data['results'][0]['primaryAccession'] if data['results'] else None

def get_uniprot_length(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}&fields=length"
    response = requests.get(url)
    response.raise_for_status()
    data = json.loads(response.text)
    if data['results']:
        return data['results'][0]['sequence']['length']
    return "Length not found"

def get_alphafold_prediction(uniprot_id):
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    response = requests.get(url)
    if response.status_code == 200:
        return f"AlphaFold prediction available: https://alphafold.ebi.ac.uk/entry/{uniprot_id}"
    else:
        return "No AlphaFold prediction found"

def get_uniprot_3d(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}&fields=xref_pdb"
    response = requests.get(url)
    response.raise_for_status()
    data = json.loads(response.text)
    structures = []
    
    pdb_count = 0
    
    if data['results'] and 'uniProtKBCrossReferences' in data['results'][0]:
        pdb_structures = data['results'][0]['uniProtKBCrossReferences']
        for pdb in pdb_structures:
            if pdb['database'] == 'PDB':
                pdb_count += 1
                structure_info = (
                    f"PDB ID: {pdb['id']}, Method: {pdb['properties'][0]['value']}, "
                    f"Resolution: {pdb['properties'][1]['value']}, Chains: {pdb['properties'][2]['value']}"
                )
                structures.append(structure_info)
    
    if not structures:
        structures.append("No experimental 3D structures found")
    
    # Check AlphaFold prediction
    alphafold_available = get_alphafold_prediction(uniprot_id)
    
    if alphafold_available:
        structures.append(f"AlphaFold prediction available: https://alphafold.ebi.ac.uk/entry/{uniprot_id}")
    
    # Assign index based on criteria
    if alphafold_available and pdb_count > 0:
        index = "High"
    else:
        index = "Low"  # You can adjust this or add more conditions if needed
    
    return structures, index

def get_uniprot_disease(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}&fields=cc_disease"
    response = requests.get(url)
    response.raise_for_status()
    data = json.loads(response.text)
    
    diseases = []
    if data['results'] and 'comments' in data['results'][0]:
        disease_comments = [comment for comment in data['results'][0]['comments'] if comment['commentType'] == 'DISEASE']
        for comment in disease_comments:
            disease = comment['disease']
            diseases.append(f"Disease: {disease['diseaseId']}. {disease['description']}")
    
    if diseases:
        return "\n".join(diseases)
    return "No disease information found"