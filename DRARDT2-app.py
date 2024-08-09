import streamlit as st
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

def get_uniprot_id(gene_name):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_name}&fields=accession"
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
            diseases.append(f"<b style='font-size: 18px;'Disease: </b> {disease['diseaseId']}.\n </b> {disease['description']}")
    
    if diseases:
        return "\n".join(diseases)
    return "No disease information found"


###################### Streamlit app code ######################

def main():
    st.title("DRARDT: Drug Repurposing Assessment for Rare Disease Targets ")

    st.markdown(f"<i style='font-size: 24px;'>A method to evaluate the druggability potential of rare disease target proteins</i>", unsafe_allow_html=True)

    # Display the logo centered using st.image
    logo = "./cassmedchem_logo.png"  # Replace with your logo URL or local file path
    st.image(logo, use_column_width='auto', caption="DRARDT was developed at CASSMedChem (University of Turin, Italy).\n Website: https://www.cassmedchem.unito.it/")
    
    st.markdown(f"<div style='font-size: 20px;'>Understanding the targetability potential of rare disease targets is a requirement to  pursue drug discovery campaigns aiming at treating these disorders, which are frequently neglected by pharmaceutical companies. </div>", unsafe_allow_html=True)
    st.markdown(f"<div style='font-size: 20px;'>DRARDT takes into account the following aspects to assign targetability indexes to RD targets: </div>", unsafe_allow_html=True)
    
    # Display the DRARDT scheme centered using st.image
    scheme = "./drardt.jpg"
    st.image(scheme, use_column_width='auto', caption="Add ref ProtSci")

    # Input fields for gene name and email
    gene_name = st.text_input("Enter gene name:")
    # Set your email here (required by Entrez)
    email = st.text_input("Enter a valid email address:")

    # Button to trigger data processing
    if st.button("Submit"):
        if not email or not gene_name:
            st.error("Please provide both gene name and email address.")
        else:
            # Set email for Entrez
            Entrez.email = email 
            
            # Process the data
            uniprot_id = get_uniprot_id(gene_name)
            st.markdown(f"<b style='font-size: 20px;'>UniProt ID for {gene_name}:</b>", unsafe_allow_html=True)
            st.markdown(f"<div style='font-size: 18px;'>{uniprot_id}</div>", unsafe_allow_html=True)
            st.markdown("<br>", unsafe_allow_html=True)

            if uniprot_id:
                pub_count = get_publication_count(gene_name)
                st.markdown(f"<b style='font-size: 20px;'>Number of publications about {gene_name} from 2000 to present:</b>", unsafe_allow_html=True)
                st.markdown(f"<div style='font-size: 18px;'>{pub_count}</div>", unsafe_allow_html=True)
                st.markdown("<br>", unsafe_allow_html=True)
               
                pub_index = assign_publication_index(pub_count)
                st.markdown(f"<b style='font-size: 20px;'>Publication index for {gene_name}:</b>", unsafe_allow_html=True)
                st.markdown(f"<div style='font-size: 18px;'>{pub_index}</div>", unsafe_allow_html=True)
                st.markdown("<br>", unsafe_allow_html=True)

                pathways_count = get_kegg_pathways(gene_name)
                st.markdown(f"<b style='font-size: 20px;'>Number of pathways {gene_name} is involved in:</b>", unsafe_allow_html=True)
                st.markdown(f"<div style='font-size: 18px;'>{pathways_count}</div>", unsafe_allow_html=True)
                st.markdown("<br>", unsafe_allow_html=True)

                interactors_count = get_string_interactors(gene_name)
                st.markdown(f"<b style='font-size: 20px;'>Number of {gene_name} interactors from the STRING-DB:</b>", unsafe_allow_html=True)
                st.markdown(f"<div style='font-size: 18px;'>{interactors_count}</div>", unsafe_allow_html=True)
                st.markdown("<br>", unsafe_allow_html=True)

                interactome_index = assign_interactome_index(pathways_count, interactors_count)
                st.markdown(f"<b style='font-size: 20px;'>Interactome index for {gene_name}:</b>", unsafe_allow_html=True)
                st.markdown(f"<div style='font-size: 18px;'>{interactome_index}</div>", unsafe_allow_html=True)
                st.markdown("<br>", unsafe_allow_html=True)

                uniprot_prot_length = get_uniprot_length(uniprot_id)    
                st.markdown(f"<b style='font-size: 20px;'>Length of {gene_name}:</b>", unsafe_allow_html=True)
                st.markdown(f"<div style='font-size: 18px;'>{uniprot_prot_length}</div>", unsafe_allow_html=True)
                st.markdown("<br>", unsafe_allow_html=True)

                uniprot_3d_struct, structure_index = get_uniprot_3d(uniprot_id)
                st.markdown(f"<b style='font-size: 20px;'>3D structure information for {gene_name}:</b>", unsafe_allow_html=True)
                for structure in uniprot_3d_struct:
                    st.markdown(f"<div style='font-size: 18px;'>{structure}</div>", unsafe_allow_html=True)
                    st.markdown("<br>", unsafe_allow_html=True)
                
                st.markdown(f"<b style='font-size: 20px;'>Structure index for {gene_name}:</b>", unsafe_allow_html=True)
                st.markdown(f"<div style='font-size: 18px;'>{structure_index}</div>", unsafe_allow_html=True)
                st.markdown("<br>", unsafe_allow_html=True)

                uniprot_pathology = get_uniprot_disease(uniprot_id)
                st.markdown(f"<b style='font-size: 20px;'>Involvement of {gene_name} in diseases:</b>", unsafe_allow_html=True)
                diseases = uniprot_pathology.split("\n")
                for disease in diseases:
                    st.markdown(f"<div style='font-size: 18px;'>{disease}</div>", unsafe_allow_html=True)
                    st.markdown("<br>", unsafe_allow_html=True)  # Add space between diseases

if __name__ == "__main__":
    main()