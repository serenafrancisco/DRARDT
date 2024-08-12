import xml.etree.ElementTree as ET
import pandas as pd
import requests

def fetch_variant_ids(api_url):
    response = requests.get(api_url)
    if response.status_code == 200:
        try:
            root = ET.fromstring(response.text)
            id_list = root.find('IdList')
            return [id_elem.text for id_elem in id_list.findall('Id')] if id_list is not None else []
        except ET.ParseError as e:
            print(f"Failed to parse XML: {e}")
    else:
        print(f"Request failed with status code {response.status_code}")
    return []

def split_variation_name(variation_name):
    # Assuming the format: NM_004722.4(AP4M1):c.1259A>G (p.Gln420Arg)
    gene_info = variation_name.split(':')[0].strip()
    rest = variation_name.split(':')[1].strip()
    nucleotide_change = rest.split(' ')[0].strip()
    amino_acid_change = rest.split(' ')[1].strip() if len(rest.split(' ')) > 1 else ''
    # Remove parentheses from amino acid change
    if '(' in amino_acid_change and ')' in amino_acid_change:
        amino_acid_change = amino_acid_change.split('(')[-1].split(')')[0]
    return gene_info, nucleotide_change, amino_acid_change

def fetch_variant_details(variant_ids):
    df_list = []
    for variant_id in variant_ids:
        api_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={variant_id}&retmode=json'
        response = requests.get(api_url)
        if response.status_code == 200:
            summary_data = response.json()
            result_key = str(variant_id)
            if result_key in summary_data['result']:
                variant_detail = summary_data['result'][result_key]
                variation_name = variant_detail.get('title', '')
                gene_info, nucleotide_change, amino_acid_change = split_variation_name(variation_name)
                data = {
                    'UID': [variant_detail.get('uid', '')],
                    'Transcript(Gene)': [gene_info],
                    'Nucleotide Change': [nucleotide_change],
                    'Amino Acid Change': [amino_acid_change]
                }
                df_list.append(pd.DataFrame(data))
            else:
                print(f"Variant ID {variant_id} not found in the response.")
        else:
            print(f"Request for Variant_ID {variant_id} failed with status code {response.status_code}")
    return pd.concat(df_list, ignore_index=True) if df_list else pd.DataFrame()