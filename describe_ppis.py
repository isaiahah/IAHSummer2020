import math
import numpy as np
import pandas as pd

# This file contains a set of operations to combine, clean, and analyze 
# datasets. Each operation is its own block with a short explanation. 
# Uncomment the desired operations to run.
# Some methods are OUTDATED but included for reference. DO NOT USE them.

# If you see an Exception raised, check you aren't calling an outdated method
# before reporting the error.

# Adding the RNA data to the PPI data
"""
ppi_data = pd.read_csv("AP-MS/PPI-Output/Merged-PPI-List-filtered.csv")
rs_data = pd.read_csv("RNASeq/rna-seq-combined.csv")

rs_data = rs_data.rename(columns = {'GeneName' : 'PreyGene', 'log2ChangeA549' : 'l2CA549', 
'padjA549' : 'padjA549', 'log2ChangeNHBE' : 'l2CNHBE', 'padjNHBE': 'padjNHBE',
'log2ChangeAvg' : 'l2CAvg', 'padjAvg' : 'padjAvg'})

ppi_data = ppi_data.join(rs_data.set_index('PreyGene'), on = 'PreyGene')

ppi_data.to_csv("AP-MS/PPI-Output/ppi-rna/ppi-rna-filtered.csv", index = False, na_rep = 'NaN')
"""


# An earlier attempt at counting proteins by PPIs. DO NOT USE, OUTDATED
"""
raise Exception
def df_dict(df, rowname, colname):
    ''' (pandas DataFrame, str row name in DataFrame, str col name in DataFrame) -> Value in DataFrame
    
    Return an item in the input DataFrame based on row and column, but will not throw an error. 
    Why did I make this?
    '''
    if rowname in df.index:
        return df.at[rowname, colname]
    else:
        return None

def count(list_of_col, file):
    ''' (list of pandas Series, filepath) -> NoneType
    
    Counts the number of occurrences of each item in a list of Series, and writes it to the file given.
    The File is written as <items> : <count> for all items, each on a different line.    
    '''
    
    total_proteins = pd.Series()  # Why is this a Series? Should have used a list.
    for col in list_of_col:
        total_proteins = total_proteins.append(col)
    protein_count = {}
    for protein in total_proteins:
        if protein in protein_count:
            protein_count[protein] += 1
        else:
            protein_count[protein] = 1
    protein_count_list = []
    for protein in protein_count:
        protein_count_list.append([protein_count[protein], protein])
    protein_count_list.sort(reverse = True)

    outputfile = open(file, 'wt')
    outputfile.write("ALL PROTEINS\n---\n")
    for pair in protein_count_list:
        write_str = str(pair[1]) + " : " + str(pair[0]) + "\n"
        outputfile.write(write_str)
    outputfile.close()

count([ppi_data.Bait, ppi_data.PreyGene], "AP-MS/PPI-Output/protein-rank.txt")
"""

# Adding Ash's combined PPI Data
"""
CoV_PPIs = pd.read_csv('processed/merged_ppi_lists/Merged-CoV-PPI-List.csv', skiprows= 1, usecols = ["BaitGene", "BaitIsCoV", "Prey", "PreyGene"])
ash_data = pd.read_csv('AP-MS/Ash_Combined/SARS_MERS_SARS2 Known PPIs.csv', usecols = ["Virus", "Host gene ID", "Host gene symbol", "Virus gene"])[["Virus", "Virus gene", "Host gene ID", "Host gene symbol"]]
ash_data = ash_data[ash_data.Virus == "SARS-CoV-2"]
ash_data = ash_data.drop(columns = 'Virus')
ash_data = ash_data.rename(columns = {"Virus gene" : "BaitGene", "Host gene ID" : "Prey", 
                                      "Host gene symbol" : "PreyGene"})
ash_data.insert(1, "BaitIsCoV", [1] * len(ash_data.BaitGene))
CoV_PPIs = CoV_PPIs.append(ash_data)
CoV_PPIs.to_csv('processed/merged_ppi_lists/Merged-CoV-PPI-List.csv', index = False, na_rep = "NaN")
"""

# Adding Human PPIs to a sourced file
"""
input_file = 'processed/merged_ppi_lists/Roth1-Past-PPI-Sourced.csv'
CoV_PPIs = pd.read_csv(input_file)

BioPlex_PPIs = pd.read_csv('Human_PPIs/Bioplex3.tsv', skiprows = 1, delim_whitespace = True, 
                           usecols = ["SymbolA", "SymbolB"])[["SymbolA", "SymbolB"]]
BioPlex_PPIs = BioPlex_PPIs.rename(columns = {"SymbolA" : "BaitSymbol", "SymbolB" : "PreySymbol"})
BioPlex_PPIs.insert(1, "BaitIsCoV", 0)
BioPlex_PPIs.insert(3, "PPISymbol", BioPlex_PPIs.BaitSymbol + "_" + BioPlex_PPIs.PreySymbol)
BioPlex_PPIs.insert(4, "Source", "BioPlex")
BioPlex_Rel_PPIs = BioPlex_PPIs[BioPlex_PPIs.BaitSymbol.isin(CoV_PPIs.PreySymbol)]
BioPlex_Rel_PPIs = BioPlex_Rel_PPIs[BioPlex_Rel_PPIs.PreySymbol.isin(CoV_PPIs.PreySymbol)]
BioPlex_Rel_PPIs = BioPlex_Rel_PPIs.drop_duplicates(subset = ["BaitSymbol", "PreySymbol"])
CoV_PPIs = CoV_PPIs.append(BioPlex_Rel_PPIs)

QUBIC = pd.read_csv('Human_PPIs/QUBIC.csv', usecols = ["bait.Gene.name", "prey.Gene.name"])[["bait.Gene.name", "prey.Gene.name"]]
QUBIC = QUBIC.rename(columns = {"bait.Gene.name" : "BaitSymbol", "prey.Gene.name" : "PreySymbol"})
QUBIC.insert(1, "BaitIsCoV", 0)
QUBIC.insert(3, "PPISymbol", QUBIC.BaitSymbol + "_" + QUBIC.PreySymbol)
QUBIC.insert(4, "Source", "QUBIC")
QUBIC_Rel_PPIs = QUBIC[QUBIC.BaitSymbol.isin(CoV_PPIs.PreySymbol)]
QUBIC_Rel_PPIs = QUBIC_Rel_PPIs[QUBIC_Rel_PPIs.PreySymbol.isin(CoV_PPIs.PreySymbol)]
QUBIC_Rel_PPIs = QUBIC_Rel_PPIs.drop_duplicates(subset = ["BaitSymbol", "PreySymbol"])
CoV_PPIs = CoV_PPIs.append(QUBIC_Rel_PPIs)

CoFrac = pd.read_csv('Human_PPIs/CoFrac/nature14871-s2/CoFrac_PPIs.csv', 
                     usecols= ["Name1", "Name2"])[["Name1", "Name2"]]
CoFrac = CoFrac.rename(columns = {"Name1" : "BaitSymbol", "Name2" : "PreySymbol"})
CoFrac_Rel_PPIs = CoFrac[CoFrac.BaitSymbol.isin(CoV_PPIs.PreySymbol)]
CoFrac_Rel_PPIs = CoFrac_Rel_PPIs[CoFrac_Rel_PPIs.PreySymbol.isin(CoV_PPIs.BaitSymbol)]
CoFrac_Rel_PPIs.insert(1, "BaitIsCoV", 0)
CoFrac_Rel_PPIs.insert(3, "PPISymbol", CoFrac.BaitSymbol + "_" + CoFrac.PreySymbol)
CoFrac_Rel_PPIs.insert(4, "Source", "CoFrac")
CoFrac_Rel_PPIs = CoFrac_Rel_PPIs.drop_duplicates(subset = ["BaitSymbol", "PreySymbol"])
CoV_PPIs = CoV_PPIs.append(CoFrac_Rel_PPIs)

Lit_BM = pd.read_csv('Human_PPIs/HuRI/Supplementary_Tables_new/Lit-BM-17.csv',
                     usecols = ["BaitGene", "PreyGene"])
Lit_BM = Lit_BM.rename(columns = {"BaitGene" : "BaitSymbol", "PreyGene" : "PreySymbol"})
Lit_BM_Rel_PPIs = Lit_BM[Lit_BM.BaitSymbol.isin(CoV_PPIs.PreySymbol)]
Lit_BM_Rel_PPIs = Lit_BM_Rel_PPIs[Lit_BM_Rel_PPIs.PreySymbol.isin(CoV_PPIs.PreySymbol)]
Lit_BM_Rel_PPIs.insert(1, "BaitIsCoV", 0)
Lit_BM_Rel_PPIs.insert(3, "PPISymbol", Lit_BM_Rel_PPIs.BaitSymbol + "_" + Lit_BM_Rel_PPIs.PreySymbol)
Lit_BM_Rel_PPIs.insert(4, "Source", "Lit-BM-17")
Lit_BM_Rel_PPIs = Lit_BM_Rel_PPIs.drop_duplicates(subset = ["BaitSymbol", "PreySymbol"])
CoV_PPIs = CoV_PPIs.append(Lit_BM_Rel_PPIs)

HuRI = pd.read_csv('Human_PPIs/HuRI/Supplementary_Tables_new/HuRI.csv',
                   usecols = ["BaitGene", "PreyGene"])
HuRI = HuRI.rename(columns = {"BaitGene" : "BaitSymbol", "PreyGene" : "PreySymbol"})
HuRI_Rel_PPIs = HuRI[HuRI.BaitSymbol.isin(CoV_PPIs.PreySymbol)]
HuRI_Rel_PPIs = HuRI_Rel_PPIs[HuRI_Rel_PPIs.PreySymbol.isin(CoV_PPIs.PreySymbol)]
HuRI_Rel_PPIs.insert(1, "BaitIsCoV", 0)
HuRI_Rel_PPIs.insert(3, "PPISymbol", HuRI_Rel_PPIs.BaitSymbol + "_" + HuRI_Rel_PPIs.PreySymbol)
HuRI_Rel_PPIs.insert(4, "Source", "HuRI")
HuRI_Rel_PPIs = HuRI_Rel_PPIs.drop_duplicates(subset = ["BaitSymbol", "PreySymbol"])
CoV_PPIs = CoV_PPIs.append(HuRI_Rel_PPIs)

output_file = 'processed/merged_ppi_lists/Roth1-Past-Human-PPI-Sourced.csv'
CoV_PPIs.to_csv(output_file, na_rep = "NaN", index = False)
"""

# An OUTDATED method to add RNAseq from 
# https://www.biorxiv.org/content/10.1101/2020.03.24.004655v1.full.pdf to Human + CoV PPIs
# Re-written to use published data, which is in better formatting.
# Listed here for reference.
# DO NOT USE
"""
raise Exception
ppi_data = pd.read_csv('processed/merged_ppi_lists/CoV_Human_PPIs.csv')
rna_data = pd.read_csv('RNAseq/rna-seq-combined.csv', 
                       usecols= ["GeneName", "log2ChangeA549", "log2ChangeNHBE", "log2ChangeAvg"])
rna_data = rna_data.rename(columns = {'GeneName' : 'PreySymbol', 'log2ChangeA549' : 'l2CA549', 
                                      'log2ChangeNHBE' : 'l2CNHBE', 'log2ChangeAvg' : 'l2CAvg'})
ppi_data = ppi_data.join(rna_data.set_index('PreySymbol'), on = 'PreySymbol')

ppi_data.to_csv('processed/merged_ppi_lists/ppi-rna/ppi-rna-CoV-Human.csv', index = False, 
                na_rep = 'NaN')
"""

# Removing duplicate PPIs with a switched order
"""
ppis = pd.read_csv('processed/merged_ppi_lists/Roth1-Past-Human-PPI.csv')
ppi_dict = {}
drop_rows = []
for row in ppis.index:
    bait = ppis.loc[row, "BaitSymbol"]
    prey = ppis.loc[row, "PreySymbol"]
    if bait in ppi_dict:
        if prey in ppi_dict[bait]:
            drop_rows.append(row)
        else:
            ppi_dict[bait].append(prey)
    else:
        if prey in ppi_dict:
            if bait in ppi_dict[prey]:
                drop_rows.append(row)
            else:
                ppi_dict[prey].append(bait)
        else:
            ppi_dict[bait] = [prey]

ppis = ppis.drop(drop_rows)
ppis.to_csv('processed/merged_ppi_lists/Roth1-Past-Human-PPI2.csv', index = False, na_rep = "NaN")
"""

# Dropping NA baits and preys, which likely entered from the original dataset
# being formatted in Excel. 
"""
ppis = pd.read_csv('processed/merged_ppi_lists/Roth1-Past-PPI.csv', na_values = "NaN")
ppis = ppis[ppis.BaitSymbol.notna()]
ppis = ppis[ppis.PreySymbol.notna()]
ppis.to_csv('processed/merged_ppi_lists/Roth1-Past-PPI.csv', index = False, na_rep = "NaN")
"""


# Adding the new RNA-seq data from https://doi.org/10.1016/j.cell.2020.04.026 and 
# https://doi.org/10.1038/s41586-020-2332-7
"""
mb_mmc1 = pd.read_csv('RNASeq/blanco-mello/bm-mmc1.csv', 
                      usecols = ["GeneName", "SARS-CoV-2(Calu-3)_L2FC", "SARS-CoV-2(A549)_L2FC",
                                 "SARS-CoV-2(A549-ACE2)LowMOI_L2FC", "SARS-CoV-2(A549-ACE2)HiMOI_L2FC",
                                 "SARS-CoV-2(A549-ACE2)-Ruxolitinib_L2FC"],
                      na_values = ["#NUM!", "#DIV/0!"])
mb_mmc1 = mb_mmc1.rename(columns = {"GeneName" : "PreySymbol", 
                                    "SARS-CoV-2(Calu-3)_L2FC": "L2F(Calu-3)", 
                                    "SARS-CoV-2(A549)_L2FC" : "L2F(A549)",
                                    "SARS-CoV-2(A549-ACE2)LowMOI_L2FC" : "L2F(A549-ACE2-LowMOI)",
                                    "SARS-CoV-2(A549-ACE2)HiMOI_L2FC" : "L2F(A549-ACE2-HighMOI)",
                                    "SARS-CoV-2(A549-ACE2)-Ruxolitinib_L2FC" : "L2F(A549-ACE2-Ruxolitinib"})

mb_mmc2 = pd.read_csv('RNASeq/blanco-mello/bm-mmc2.csv',
                      usecols = ["GeneName", "SARS-CoV-2_L2FC"],
                      na_values = ["#NUM!", "#DIV/0!"])
mb_mmc2 = mb_mmc2.rename(columns = {"GeneName" : "PreySymbol", "SARS-CoV-2_L2FC" : "L2F(NHBE)"})

mb_mmc3_ferret = pd.read_csv('RNASeq/blanco-mello/bm-mmc3-ferret.csv', 
                             usecols = ["GeneSymbol", "SARS-CoV-2 Day1_L2FC", 
                                        "SARS-CoV-2 Day3_L2FC", "SARS-CoV-2 Day7_L2FC", 
                                        "SARS-CoV-2 Day14_L2FC"],
                             na_values = ["#NUM!", "#DIV/0!"])
mb_mmc3_ferret = mb_mmc3_ferret.rename(columns = {"GeneSymbol" : "PreySymbol",
                                                  "SARS-CoV-2 Day1_L2FC" : "L2F(FerretD1)",
                                                  "SARS-CoV-2 Day3_L2FC" : "L2F(FerretD3)",
                                                  "SARS-CoV-2 Day7_L2FC" : "L2F(FerretD7)",
                                                  "SARS-CoV-2 Day14_L2FC" : "L2F(FerretD14)"})

mb_mmc3_trachea = pd.read_csv('RNASeq/blanco-mello/bm-mmc3-trachea.csv',
                              usecols = ["GeneSymbol", "SARS-CoV-2 Day3_L2FC"],
                              na_values = ["#NUM!", "#DIV/0!"])
mb_mmc3_trachea = mb_mmc3_trachea.rename(columns = {"GeneSymbol" : "PreySymbol", 
                                                    "SARS-CoV-2 Day3_L2FC" : "L2F(Trachea)"})

mb_mmc4 = pd.read_csv('RNASeq/blanco-mello/bm-mmc4.csv',
                      usecols = ["Gene_name", "log2FoldChange"],
                      na_values = ["#NUM!", "#DIV/0!"])
mb_mmc4 = mb_mmc4.rename(columns = {"Gene_name" : "PreySymbol", "log2FoldChange" : "L2F(Patient)"})

bojkova_proteome = pd.read_csv('RNASeq/bojkova/proteome-time-data.csv',
                               usecols = ["Gene Symbol", "Ratio 2h", "Ratio 6h", "Ratio 10h ", "Ratio 24h"],
                               na_values = ["#NUM!", "#DIV/0!"])
bojkova_proteome = bojkova_proteome.rename(columns = {"Gene Symbol" : "PreySymbol", "Ratio 2h" : "protL2F_2h",
                                            "Ratio 6h" : "protL2F_6h", "Ratio 10h " : "protL2F_10h", 
                                            "Ratio 24h" : "protL2F_24h"})

bojkova_translateome = pd.read_csv('RNASeq/bojkova/translateome-time-data.csv',
                               usecols = ["Gene Symbol01", "Ratio 2h", "Ratio 6h", "Ratio 10h ", "Ratio 24h"],
                               na_values = ["#NUM!", "#DIV/0!"])
bojkova_translateome = bojkova_translateome.rename(columns = {"Gene Symbol01" : "PreySymbol", 
                                                    "Ratio 2h" : "transL2F_2h", "Ratio 6h" : "transL2F_6h", 
                                                    "Ratio 10h " : "transL2F_10h", "Ratio 24h" : "transL2F_24h"})

for dataset in ["Past-CoV-Human-PPI", "Past-CoV-PPI", "Roth1-Human-PPI", "Roth1-Past-Human-PPI",
                "Roth1-Past-PPI", "Roth1-PPI"]:
    input_file = 'processed/merged_ppi_lists/' + dataset + '.csv'
    ppis = pd.read_csv(input_file)
    
    ppis = ppis.join(mb_mmc1.set_index('PreySymbol'), on = 'PreySymbol')
    ppis = ppis.join(mb_mmc2.set_index('PreySymbol'), on = 'PreySymbol')
    ppis = ppis.join(mb_mmc3_ferret.set_index('PreySymbol'), on = "PreySymbol")
    ppis = ppis.join(mb_mmc3_trachea.set_index('PreySymbol'), on = "PreySymbol")
    ppis = ppis.join(mb_mmc4.set_index('PreySymbol'), on = 'PreySymbol')
    ppis = ppis.join(bojkova_translateome.set_index('PreySymbol'), on = 'PreySymbol')
    ppis = ppis.join(bojkova_proteome.set_index('PreySymbol'), on = 'PreySymbol')
    ppis = ppis.drop_duplicates(subset = ["BaitSymbol", "PreySymbol"])

    output_file = 'processed/merged_ppi_lists/ppi-rna/RNA-' + dataset + '.csv'
    ppis.to_csv(output_file, index = False, na_rep = "NaN")
"""

# Adds phosphoprotemomics data from https://doi.org/10.1101/2020.06.17.156455
"""
phospho = pd.read_csv('MMann/phosphoproteomics-1.csv', skiprows = 1, 
                      usecols = ["site_ID", "residue", "position", "gene name", 
                                 "significant_hit-6", "fold_change_vs_mock_(log2)-6",
                                 "significant_hit-24", "fold_change_vs_mock_(log2)-24", "known site"])
phospho = phospho[["gene name", "site_ID", "residue", "position", "significant_hit-6", 
                   "fold_change_vs_mock_(log2)-6", "significant_hit-24", 
                   "fold_change_vs_mock_(log2)-24", "known site"]]
phospho = phospho.rename(columns = {"site_ID" : "site_ID", "residue" : "residue", 
                                    "position" : "position", "gene name" : "PreySymbol", 
                                    "significant_hit-6" : "P_6_hit", 
                                    "fold_change_vs_mock_(log2)-6" : "P_6_L2F", 
                                    "significant_hit-24" : "P_24_hit", 
                                    "fold_change_vs_mock_(log2)-24" : "P_24_L2F",
                                    "known site" : "known_site"})
phospho.P_6_hit = phospho.P_6_hit == "+"
phospho.P_24_hit = phospho.P_24_hit == "+"
phospho.known_site = phospho.known_site == "+"

for dataset in ["Past-CoV-Human-PPI", "Past-CoV-PPI", "Roth1-Human-PPI", "Roth1-Past-Human-PPI",
                "Roth1-Past-PPI", "Roth1-PPI"]:
    dataset_input = 'processed/merged_ppi_lists/' + dataset + '.csv'
    ppis = pd.read_csv(dataset_input)
    ppis = ppis.merge(phospho.set_index("PreySymbol"), on = "PreySymbol")
    dataset_output = 'processed/merged_ppi_lists/phospho_prot/' + dataset + '.csv'
    ppis.to_csv(dataset_output, index = False, na_rep = "NaN")
"""

# Setting NAs to be consistent in the data. Again, thanks Excel.
"""
for dataset in ["Past-CoV-Human-PPI-List", "Past-CoV-PPI", "Roth1-Human-PPI", "Roth1-Past-Human-PPI",
                "Roth1-Past-PPI", "Roth1-PPI"]:
    dataset_filepath = 'processed/merged_ppi_lists/ppi-rna/RNA-' + dataset + '.csv'
    ppis = pd.read_csv(dataset_filepath, na_values = ["#NUM!", "#DIV/0!"])
    ppis.to_csv(dataset_filepath, index = False, na_rep = "NaN")
"""

# Adding RNA-seq Rank Data
"""
for dataset in ['processed/merged_ppi_lists/ppi-rna/RNA-Roth1-PPI',
                'processed/merged_ppi_lists/ppi-rna/RNA-Roth1-Past-PPI',
                'processed/merged_ppi_lists/ppi-rna/RNA-Roth1-Past-Human-PPI',
                'processed/merged_ppi_lists/ppi-rna/RNA-Roth1-Human-PPI']:
    ppis = pd.read_csv(dataset + '.csv', na_values = ["#DIV/0!", "#NUM!"])
    ppis.insert(len(ppis.columns), "RNA_Rank", 
                ((ppis.L2F_A549 - ppis.L2F_A549.mean()) / (ppis.L2F_A549.std(ddof = 0))))
    ppis.to_csv(dataset + '.csv', index = False, na_rep = "NaN")
"""

# Counting frequency of 'important' proteins in significant sets
"""
for dataset in ['RNA-Roth1-PPI', 'RNA-Roth1-Past-PPI', 
                'RNA-Roth1-Past-Human-PPI', 'RNA-Roth1-Human-PPI']:
    ppis = pd.read_csv('processed/merged_ppi_lists/subsets/network_annot/' + dataset + '.csv')
    ppis = ppis.drop_duplicates(subset = ["PreySymbol"])
    freq = ppis[ppis.Important > 0].Important.count() / ppis.PreySymbol.count()
    print(dataset + ":", freq)
"""

# Creating a non-sourced file from a sourced file
"""
for dataset in ["Past-CoV-Human-PPI", "Past-CoV-PPI", "Roth1-Human-PPI", "Roth1-Past-Human-PPI",
                "Roth1-Past-PPI"]:
    input_file = 'processed/merged_ppi_lists/' + dataset + '-Sourced.csv'
    ppis = pd.read_csv(input_file)
    ppis = ppis.drop(columns = ["Source"])
    ppis = ppis.drop_duplicates(subset = ["BaitSymbol", "PreySymbol"])
    output_file = 'processed/merged_ppi_lists/' + dataset + '.csv'
    ppis.to_csv(output_file, index = False, na_rep = "NaN")
"""

# Generate file of only Phosphoproteomic Proteins
"""
phospho_ppi = pd.read_csv('MMann/phosphoproteomics.csv', skiprows = 1,
                          usecols = ["gene name"])
phospho_ppi = phospho_ppi.drop_duplicates()
phospho_ppi.to_csv('MMann/phosphoproteomics_prot.csv', index = False)
"""

# Convert the Falter-Braun list to a PPI list
'''
fb_prot = pd.read_table('AP-MS/Falter-Braun/FB_Gene_List.tab', 
                        usecols = ["Gene_names_primary"])
fb_prot = fb_prot.rename(columns = {"Gene_names_primary" : "PreySymbol"})
fb_prot.to_csv('AP-MS/Falter-Braun/Falter-Braun_PPIs.csv', 
               index = False, na_rep = "NaN")
'''
