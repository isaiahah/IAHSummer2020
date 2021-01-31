import pandas as pd

def extract_prot(filepath):
    # Given a set of all COVID protein - other protein interactions, create a file for each COVID protein
    # detailing only its interactions
    ppis = pd.read_csv(filepath)

    for protein in ["NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6", "NSP7", "NSP8", "NSP8", "NSP9", 
                "NSP10", "NSP11", "NSP12", "NSP13", "NSP14", "NSP15", "NSP16", "S", "ORF3A", "ORF3B",
                "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF9A", "ORF9B", "ORF10"]:
        protein_dataframe = ppis[ppis.BaitSymbol == protein]
        first_deg_temp = ppis[ppis.PreySymbol.isin(protein_dataframe.PreySymbol)]
        first_deg = first_deg_temp[first_deg_temp.BaitSymbol.isin(protein_dataframe.PreySymbol)]
        protein_dataframe = protein_dataframe.append(first_deg)
        protein_dataframe = protein_dataframe.drop_duplicates()
        protein_dataframe.to_csv('processed/merged_ppi_lists/ppi-rna/CoV_protein_maps/' + protein + '.csv', index = False, na_rep = "NaN")    
