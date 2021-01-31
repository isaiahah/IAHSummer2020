import mygene
import pandas as pd

# This file uses MyGene's api to standardize gene symbols in the HuRI dataset.
mg = mygene.MyGeneInfo()

# Convert the Ensembl IDs in HuRI's Lit-BM document to Gene Symbols.
# Uncomment to run.
"""
Lit_BM = pd.read_csv('Human_PPIs/HuRI/Supplementary_Tables_new/Supplementary Table 14.txt', 
                     usecols = ["ensembl_gene_id_a", "ensembl_gene_id_b"], 
                     delim_whitespace = True)

Lit_BM = Lit_BM.rename(columns = {"ensembl_gene_id_a" : "BaitID", 
                                  "ensembl_gene_id_b" : "PreyID"})
Lit_BM.insert(1, "BaitIsCoV", [0] * len(Lit_BM.BaitID))
Lit_BM.insert(2, "BaitSymbol", mg.querymany(Lit_BM.BaitID, scopes = 'ensemblgene', 
                                            species = 'human', fields = "symbol", 
                                            as_dataframe = True).symbol.reset_index(drop = True))
Lit_BM.insert(4, "PreySymbol", mg.querymany(Lit_BM.PreyID, scopes = 'ensemblgene', 
                                            species = "human", fields = "symbol", 
                                            as_dataframe = True).symbol.reset_index(drop = True))

Lit_BM.to_csv('Human_PPIs/HuRI/Supplementary_Tables_new/Lit-BM.csv', index = False, na_rep = "NaN")
"""

# Convert the ensembl IDs in HuRI to Gene Symbols
# Uncomment to run.
"""
HuRI = pd.read_csv('Human_PPIs/HuRI/Supplementary_Tables_new/Supplementary Table 9.txt', 
                     usecols = ["Ensembl_gene_id_a", "Ensembl_gene_id_b"], 
                     delim_whitespace = True)

HuRI = HuRI.rename(columns = {"Ensembl_gene_id_a" : "BaitID", 
                                  "Ensembl_gene_id_b" : "PreyID"})
HuRI.insert(1, "BaitIsCoV", [0] * len(HuRI.BaitID))
HuRI.insert(2, "BaitSymbol", mg.querymany(HuRI.BaitID, scopes = 'ensemblgene', 
                                            species = 'human', fields = "symbol", 
                                            as_dataframe = True).symbol.reset_index(drop = True))
HuRI.insert(4, "PreySymbol", mg.querymany(HuRI.PreyID, scopes = 'ensemblgene', 
                                            species = "human", fields = "symbol", 
                                            as_dataframe = True).symbol.reset_index(drop = True))

HuRI.to_csv('Human_PPIs/HuRI/Supplementary_Tables_new/HuRI.csv', index = False, na_rep = "NaN")
"""

# Generate Entrez Gene IDs from the Gene Symbols in the CoV_Human_PPI data
# Uncomment to run.
'''
input_file = 'processed/merged_ppi_lists/Roth1-Human-PPI.csv'
ppis = pd.read_csv(input_file)
# 7079 is the last line of CoV in Roth1-Past-Human-PPIs.csv. 7077 in DF
ppis.BaitID.loc[323:] = mg.querymany(ppis.BaitSymbol.loc[323:], scopes = 'symbol', 
                               species = 'human', fields = 'entrezgene',
                               as_dataframe = True).entrezgene.reset_index(drop = True)
ppis.PreyID = mg.querymany(ppis.PreySymbol, scopes = 'symbol', 
                               species = 'human', fields = 'entrezgene',
                               as_dataframe = True).entrezgene.reset_index(drop = True)
output_file = 'processed/merged_ppi_lists/Roth1-Human-PPI2.csv'
ppis.to_csv(output_file, index = False, na_rep = "NaN")
'''