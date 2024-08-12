import requests

def get_kegg_pathways(gene_name):
    url = f"https://rest.kegg.jp/find/pathway/{gene_name}"
    response = requests.get(url)
    response.raise_for_status()
    interactors = response.text.split('\n')
    return len(interactors) - 1  # subtracting 1 to account for empty string at end

def get_string_interactors(gene_name):
    url = f"https://string-db.org/api/tsv-no-header/interaction_partners?identifiers={gene_name}&species=9606&network_type=physical"
    response = requests.get(url)
    response.raise_for_status()
    
    pathways = response.text.strip().split('\n')
    filtered_interactors = []
    
    for pathway in pathways:
        if pathway:  # Ensure we don't process empty lines
            data = pathway.split('\t')
            experimental_score = float(data[10])  # 11th column (index 10) is the experimental score
            if experimental_score > 0.700:
                filtered_interactors.append(data)
    
    return len(filtered_interactors)

def assign_interactome_index(pathways_count, interactors_count):
    if interactors_count > 3:
        return "High"
    elif pathways_count >= 1 and interactors_count < 3:
        return "Medium"
    else:
        return "Low"