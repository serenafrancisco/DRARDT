import streamlit as st
from Bio import Entrez
from time import sleep

from entrez import get_publication_count, assign_publication_index
from interactome import get_kegg_pathways, get_string_interactors, assign_interactome_index
from uniprot import get_human_uniprot_id, get_uniprot_length, get_alphafold_prediction, get_uniprot_3d, get_uniprot_disease
from clinvar import fetch_variant_ids, fetch_variant_details

def main():
    st.title("Gene Information Dashboard")

    email = st.text_input("Enter your email address:")
    gene_name = st.text_input("Enter gene name:")

    if st.button("Submit"):
        if not email or not gene_name:
            st.error("Please provide both gene name and email address.")
        else:
            Entrez.email = email

            # Fetch UniProt ID
            uniprot_id = get_human_uniprot_id(gene_name)
            if uniprot_id:
                st.write(f"**UniProt ID for {gene_name}:** {uniprot_id}")

                # Publication Count and Index
                pub_count = get_publication_count(gene_name)
                pub_index = assign_publication_index(pub_count)
                st.write(f"**Number of publications about {gene_name} (2000 to present):** {pub_count}")
                st.write(f"**Publication index for {gene_name}:** {pub_index}")

                # KEGG Pathways
                pathways_count = get_kegg_pathways(gene_name)
                st.write(f"**Number of pathways {gene_name} is involved in:** {pathways_count}")

                # STRING Interactors
                interactors_count = get_string_interactors(gene_name)
                st.write(f"**Number of interactors from STRING-DB:** {interactors_count}")

                # Interactome Index
                interactome_index = assign_interactome_index(pathways_count, interactors_count)
                st.write(f"**Interactome index for {gene_name}:** {interactome_index}")

                # UniProt Length
                uniprot_prot_length = get_uniprot_length(uniprot_id)
                st.write(f"**Length of {gene_name}:** {uniprot_prot_length}")

                # UniProt 3D Structure
                uniprot_3d_struct, structure_index = get_uniprot_3d(uniprot_id)
                st.write(f"**3D structure information for {gene_name}:**")
                for structure in uniprot_3d_struct:
                    st.write(structure)
                st.write(f"**Structure index for {gene_name}:** {structure_index}")

                # Disease Involvement
                uniprot_pathology = get_uniprot_disease(uniprot_id)
                st.write(f"**Involvement of {gene_name} in diseases:**")
                diseases = uniprot_pathology.split("\n")
                for disease in diseases:
                    st.write(disease)

                # ClinVar Variant Information
                search_api_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={gene_name}[gene]+AND+single_gene[prop]+AND+missense_variant[molcons]&retmax=900'
                variant_ids = fetch_variant_ids(search_api_url)
                if variant_ids:
                    st.write(f"**List of missense mutations for {gene_name} from ClinVar:**")
                    st.write("_This operation might take a while. We suggest you to grab a coffee :)_")
                    sleep(5)  # simulate waiting time
                    df_all_variants = fetch_variant_details(variant_ids)
                    if not df_all_variants.empty:
                        st.write(df_all_variants)
                    else:
                        st.write("No variant details available.")
                else:
                    st.write("No variant IDs found.")
            else:
                st.write("Gene not found or UniProt ID could not be retrieved.")

if __name__ == "__main__":
    main()
