import pandas as pd
import numpy as np

# This file contains a set of operations to rank proteins in datasets in importance.
# This replaces interaction_subsetting.py.
# Each operation is its own block with a short explanation. 
# Uncomment the desired operations to run.


# Add Rank columns to the data
"""
ppi_rna = pd.read_csv("AP-MS/PPI-Output/ppi-rna/ppi-rna-with-etienne.csv", na_values= "NAN")

def count(col):
    ''' (pandas Series) -> pandas Series
    
    Given a Series, counts each item's occurrences and returns them as a series paired with the item
    
    >>> count(pd.Series([1, 1, 2, 3, 4, 4]))
    pd.Series([2, 2, 1, 1, 2, 2])
    '''
    
    count_dict = {}
    for protein in col:
        if protein in count_dict:
            count_dict[protein] += 1
        else:
            count_dict[protein] = 1
    count_col = []
    for protein in col:
        count_col.append(count_dict[protein])
    return pd.Series(count_col)

ppi_rna.insert(11, "l2CAvgRank", ppi_rna.l2CAvg.rank(ascending = False, pct = True))
ppi_rna.insert(12, "padjAvgRank", ppi_rna.padjAvg.rank(ascending = False, pct = True))
ppi_rna.insert(13, "PPICount", count(ppi_rna.Prey))
ppi_rna.insert(14, "PPICountRank", ppi_rna.PPICount.rank(ascending = False, pct = True))
ppi_rna.insert(15, "rankSum", ppi_rna.l2CAvgRank + ppi_rna.padjAvgRank + ppi_rna.PPICountRank)
ppi_rna.insert(16, "squaredRankSum", (ppi_rna.l2CAvgRank ** 2 + ppi_rna.PPICountRank ** 2) / 2)
ppi_rna.to_csv('AP-MS/PPI-Output/ppi-rna/ppi-rna-with-etienne-ranked.csv', index = False, na_rep = "NaN")
"""

# Rank whole protein set
"""
ppis = pd.read_csv('AP-MS/PPI-Output/ppi-rna/ppi-rna-with-etienne-ranked.csv')
significant_list_file = open('AP-MS/PPI-Output/ppi-rna/significant_ppi.txt', 'wt')
significant_list = pd.Series()

def write_top_nper(rank_col, nper):
    ''' (Column name in PPI DataFrame, num) -> pandas Series
    
    Returns the top nper % of proteins based on the rankings in rank_col.
    Writes these proteins to significant_list_file, one per line.
    '''
    ppis_rank_col_sort = ppis.copy().set_index(rank_col).sort_index()['PreyGene']  # copy() to avoid aliasing.
    ppis_rank_col_sort_top_nper = ppis_rank_col_sort[ppis_rank_col_sort.index <= nper / 100].drop_duplicates()
    for protein in ppis_rank_col_sort_top_nper:
        significant_list_file.write(str(protein) + "\n")
    return ppis_rank_col_sort_top_nper


significant_list_file.write("PPI Count Top 1%\n")
significant_list = significant_list.append(write_top_nper('PPICountRank', 1))

significant_list_file.write("\nLog2 Change Top 1%\n")
significant_list = significant_list.append(write_top_nper('l2CAvgRank', 1))

significant_list_file.write("\np-value Top 1%\n")
significant_list = significant_list.append(write_top_nper('padjAvgRank', 1))

significant_list_file.write("\nSquare Avg PPI Log2 Count & Change Top 1%\n")
significant_list = significant_list.append(write_top_nper('squaredRankSum', 1))

significant_list_file.write("\nCombined List\n")
for protein in significant_list:
    significant_list_file.write(str(protein) + "\n")

significant_list_file.close()
"""

# Rank by Human PPIs (Awaiting Human PPI Data)
"""
ppis = pd.read_csv('AP-MS/PPI-Output/ppi-rna/ppi-rna-with-etienne-ranked.csv')
def count(list_of_col):
    ''' (List of pandas Series) -> dict
    
    Returns a dictionary where the key is items in series given in the list 
    and the value is the count of the item occurring across all lists.
    '''
    dict_count = {}
    for col in list_of_col:
        for protein in col:
            if protein in dict_count:
                dict_count[protein] += 1
            else:
                dict_count[protein] = 1
    return dict_count


def df_dict_apply(df_col, dict_apply):
    ''' (pandas Series, dict) -> pandas Series
    
    Returns a Series of the values given by using df_col as the keys to dict_apply. 
    Adds None of the key does not exits.
    '''
    return_list = []
    for item in df_col:
        if item in dict_apply:
            return_list.append(dict_apply[item])
        else:
            return_list.append(None)
    return pd.Series(return_list)
    

human_ppis = pd.read_csv('path_to_human_ppis')  # Replace this
human_ppi_proteins = human_ppis.PreyGene.append(human_ppis.Bait).drop_duplicates()
h_ppi_prot_count = pd.DataFrame({'PreyGene' : human_ppi_proteins, 
                                 'hPPICount' : df_dict_apply(human_ppis.Prey, 
                                                             count([human_ppis.PreyGene, 
                                                                    human_ppis.Bait]))})
ppis = ppis.join(h_ppi_prot_count.set_index('PreyGene'), on = 'PreyGene')
# This ends by adding human PPI counts as a row to the PPI file. More processing afterward needed.
"""

# Rank each bait's significant Preys. 
"""
significant_list_file = open('AP-MS/PPI-Output/ppi-rna/significant_ppi_per_bait.txt', 'wt')
ppis = pd.read_csv('AP-MS/PPI-Output/ppi-rna/ppi-rna-with-etienne-ranked.csv', index_col = 'Bait')

percentile = 1

for cov_protein in ppis.index.drop_duplicates():  # Finally using two loops to reduce memory usage!
    for i in range(0, 2):
        if i == 0:
            significant_list_file.write(str(cov_protein) + " PPI Count Top 1%\n")
            ppis_sorted = ppis.copy().loc[cov_protein].sort_values('PPICount')[['PreyGene', 'PPICount']]
            ppis_sorted.insert(2, 'PPIRank', ppis_sorted.PPICount.rank(ascending = False, pct = True))
            ppis_sorted = ppis_sorted.set_index('PPIRank')
            ppis_sorted_top_pct = ppis_sorted[ppis_sorted.index <= percentile / 100].drop_duplicates()
            for protein in ppis_sorted_top_pct.PreyGene:
                significant_list_file.write(str(protein) + "\n")
        if i == 1:
            significant_list_file.write(str(cov_protein) + " Log2 Change Top 1%\n")
            ppis_sorted = ppis.copy().loc[cov_protein].sort_values('l2CAvg')[['PreyGene', 'l2CAvg']]
            ppis_sorted.insert(2, 'l2CAvgRank', ppis_sorted.l2CAvg.rank(ascending = False, pct = True))
            ppis_sorted = ppis_sorted.set_index('l2CAvgRank')
            ppis_sorted_top_pct = ppis_sorted[ppis_sorted.index <= percentile / 100].drop_duplicates()
            for protein in ppis_sorted_top_pct.PreyGene:
                significant_list_file.write(str(protein) + "\n")
        significant_list_file.write("\n")
    significant_list_file.write("\n")
"""

# Finding protein repeat counts across past data.
"""
def series_list_index(l, index):
    ''' (List, index) -> Item in list
    
    Returns the list item at the given index. Only exists to use on a pandas series of lists.
    There's definitely a better way to do this, but I don't know what it is.
    '''
    return l[index]


dg = pd.read_csv('processed/merged_ppi_lists/Individual_Datasets/DGordon-High-Confidence-Pairs.csv', 
                 usecols = ["Bait", "PreyGene"])
dg = dg.rename(columns = {"Bait" : "BaitSymbol", "PreyGene" : "PreySymbol"})
dg.BaitSymbol = dg.BaitSymbol.str.split(" ").apply(series_list_index, args = (1,))
dg.BaitSymbol[dg.BaitSymbol == "Spike"] = "S"
dg.BaitSymbol = dg.BaitSymbol.str.upper()

li = pd.read_csv('processed/merged_ppi_lists/Individual_Datasets/Li-SARS-2-PPI-networks.csv',
                 usecols = ["protein", "genename"])
li = li.rename(columns = {"protein" : "BaitSymbol", "genename" : "PreySymbol"})
li.BaitSymbol[li.BaitSymbol == "spike"] = "S"
li.BaitSymbol[li.BaitSymbol == "envelope"] = "E"
li.BaitSymbol[li.BaitSymbol == "membrane"] = "M"
li.BaitSymbol[li.BaitSymbol == "nucleocapsid"] = "N"
li.BaitSymbol = li.BaitSymbol.str.upper()

mm = pd.read_csv('processed/merged_ppi_lists/Individual_Datasets/MMann-PPI-M2si.csv',
                 usecols = ["bait_full_name", "gene_name"])
mm = mm.rename(columns = {"bait_full_name" : "BaitSymbol", "gene_name" : "PreySymbol"})
mm = mm[mm.BaitSymbol.str.startswith("SARS_CoV2")]
mm.BaitSymbol = mm.BaitSymbol.str.split("_").apply(series_list_index, args = (2,))
mm.BaitSymbol = mm.BaitSymbol.str.upper()

et = pd.read_csv('AP-MS/etienne/Etienne-High-Conf.csv', skiprows = 1, 
                 usecols = ["Bait", "Gene_name"])
et = et.rename(columns = {"Bait" : "BaitSymbol", "Gene_name" : "PreySymbol"})
et.BaitSymbol = et.BaitSymbol.str.upper()

fr = pd.read_csv('processed/merged_ppi_lists/Individual_Datasets/Roth1-PPI-GoFiltered.csv',
                 usecols = ["BaitSymbol", "PreySymbol"])

ppi_count_dict = {}

for ppi_set in [dg, li, mm, et, fr]:
    for ppi_row in ppi_set.index:
        ppi_bait = ppi_set.loc[ppi_row, "BaitSymbol"]
        ppi_prey = ppi_set.loc[ppi_row, "PreySymbol"]
        if not (ppi_bait, ppi_prey) in ppi_count_dict:
            ppi_count_dict[(ppi_bait, ppi_prey)] = 1
        elif (ppi_bait, ppi_prey) in ppi_count_dict:
            ppi_count_dict[(ppi_bait, ppi_prey)] += 1

ppis_combined = dg.append(li).append(mm).append(et).append(fr).drop_duplicates(ignore_index = True)


count_col = []
for ppi_row in ppis_combined.index:
    ppi_bait = ppis_combined.loc[ppi_row, "BaitSymbol"]
    ppi_prey = ppis_combined.loc[ppi_row, "PreySymbol"]
    count_col.append(ppi_count_dict[(ppi_bait, ppi_prey)])
ppis_combined.insert(2, "Count", count_col)

output_file = 'processed/merged_ppi_lists/Past_Counts.csv'
ppis_combined.to_csv(output_file, index = False, na_rep = "NaN")
"""

# An attempt to rank based on significant prot and trans change over time
"""
for dataset in ['processed/merged_ppi_lists/ppi-rna/RNA-Roth1-PPI',
                'processed/merged_ppi_lists/ppi-rna/RNA-Roth1-Past-PPI',
                'processed/merged_ppi_lists/ppi-rna/RNA-Roth1-Past-Human-PPI',
                'processed/merged_ppi_lists/ppi-rna/RNA-Roth1-Human-PPI']:
    ppis = pd.read_csv(dataset + '.csv')
    trans_change = (ppis[["transL2F_2h", "transL2F_6h", "transL2F_10h", "transL2F_24h"]].max(axis = 1) - ppis[["transL2F_2h", "transL2F_6h", "transL2F_10h", "transL2F_24h"]].min(axis = 1))
    trans_change_rank = (trans_change - trans_change.mean()) / trans_change.std(ddof = 0)
    ppis.insert(13, "trans_change_rank", trans_change_rank)
    prot_change = (ppis[["protL2F_2h", "protL2F_6h", "protL2F_10h", "protL2F_24h"]].max(axis = 1) - 
                   ppis[["protL2F_2h", "protL2F_6h", "protL2F_10h", "protL2F_24h"]].min(axis = 1))
    
    prot_change_rank = (prot_change - prot_change.mean()) / prot_change.std(ddof = 0)
    ppis.insert(18, "prot_change_rank", prot_change_rank)

    ppis.to_csv(dataset + '.csv', index = False, na_rep = "NaN")
"""

# Generating Phosphoproteomics annotations (L2F sums and reative position)
"""
phospho_file = 'MMann/phosphoproteomics.csv'
phospho_data = pd.read_csv(phospho_file, skiprows = 1, 
                           usecols = ["residue", "position", "gene name", 
                                 "significant_hit-6", "fold_change_vs_mock_(log2)-6",
                                 "significant_hit-24", "fold_change_vs_mock_(log2)-24", "known site"])

phospho_data = phospho_data.rename(columns = {"residue" : "residue", 
                                              "position" : "position", "gene name" : "PreySymbol", 
                                              "significant_hit-6" : "P_6_hit", 
                                              "fold_change_vs_mock_(log2)-6" : "P_6_L2F", 
                                              "significant_hit-24" : "P_24_hit", 
                                              "fold_change_vs_mock_(log2)-24" : "P_24_L2F",
                                              "known site" : "known_site"})

P_6_L2F_dict = {}
P_24_L2F_dict = {}

for protein in phospho_data.PreySymbol.drop_duplicates():
    prot_phospho = phospho_data[phospho_data.PreySymbol == protein]
    prot_phospho.P_6_L2F = prot_phospho.P_6_L2F.abs()
    prot_phospho.P_24_L2F = prot_phospho.P_24_L2F.abs()
    P_6_L2F_dict[protein] = prot_phospho.P_6_L2F.sum()
    P_24_L2F_dict[protein] = prot_phospho.P_24_L2F.sum()

def dict_apply(protein, source_dict):
    if protein in source_dict:
        return source_dict[protein]
    return np.NaN

phospho_data.insert(5, "P_6_L2F_Sum", phospho_data.PreySymbol.apply(dict_apply, args = (P_6_L2F_dict,)))
phospho_data.insert(8, "P_24_L2F_Sum", phospho_data.PreySymbol.apply(dict_apply, args = (P_24_L2F_dict,)))

prot_length_data = pd.read_table('MMann/phosphoproteomics(1).tab', usecols = ["Gene_names_primary", "Length"])
prot_length_data = prot_length_data.drop_duplicates(subset = ["Gene_names_primary"])

prot_length_dict = {}
for index_id in prot_length_data.index:
    prot_length_dict[prot_length_data.loc[index_id, "Gene_names_primary"]] = prot_length_data.loc[index_id, "Length"]

prot_length = phospho_data.PreySymbol.apply(dict_apply, args = (prot_length_dict,))
phospho_data.insert(2, "rel_position", phospho_data.position / prot_length)

phospho_data.to_csv('MMann/phosphoproteomics_annotated.csv', index = False, na_rep = "NaN")
"""

# Can generate notable protein lists for: Interaction count, PPI multiplicity, RNA-seq, 
# proteome time-dependent change, translateome time-dependent change, phosphorylation
# Uncomment the metrics you want inside the section. 
# Uncomment the full section after to run.
"""
ppi_file = 'processed/merged_ppi_lists/Roth1-Past-PPI.csv'
ppi_data = pd.read_csv(ppi_file)

# Interaction Count
'''
ppi_rep_data = pd.read_csv('processed/merged_ppi_lists/Roth1-Past-PPI-Sourced.csv')
ppi_rep_data = ppi_rep_data.drop(columns = ["Index"])
ppi_rep_data = ppi_rep_data.drop_duplicates()

rep_ppis = ppi_rep_data.drop_duplicates(subset = ["PPISymbol"])
ppis_reps = []
for ppi in rep_ppis.PPISymbol:
    ppis_reps.append(ppi_rep_data.PPISymbol[ppi_rep_data.PPISymbol == ppi].count())
ppis_reps = pd.DataFrame({"BaitSymbol" : rep_ppis.BaitSymbol,
                           "PreySymbol" : rep_ppis.PreySymbol,
                           "PPISymbol" : rep_ppis.PPISymbol,
                           "PPICount" : ppis_reps})
ppis_reps = ppis_reps.sort_values(by = "PPICount", ascending = False)
ppis_reps.to_csv('processed/merged_ppi_lists/rankings/replicates.csv', 
                index = False, na_rep = "NaN")
'''

# PPI Multiplicity
'''
count_ppis = ppi_data.PreySymbol.drop_duplicates()

def count_ppi(prey):
    return ppi_data[ppi_data.PreySymbol == prey].PreySymbol.count()

prot_counts = count_ppis.apply(count_ppi)
prot_counts = pd.DataFrame({"PreySymbol" : count_ppis,
                            "PreyCount" : prot_counts})
prot_counts = prot_counts.sort_values(by = "PreyCount", ascending = False)
prot_counts.to_csv('processed/merged_ppi_lists/rankings/counts.csv', 
                   index = False, na_rep = "NaN")
'''

# RNA-seq
'''
rna_data1 = pd.read_csv('RNASeq/blanco-mello/bm-mmc1.csv', 
                       usecols = ["GeneName", "SARS-CoV-2(Calu-3)_L2FC", "SARS-CoV-2(A549)_L2FC",
                                  "SARS-CoV-2(A549-ACE2)LowMOI_L2FC", "SARS-CoV-2(A549-ACE2)HiMOI_L2FC"])
rna_data1 = rna_data1.rename(columns = {"GeneName" : "PreySymbol", 
                                        "SARS-CoV-2(Calu-3)_L2FC" : "L2F_Calu_3", 
                                        "SARS-CoV-2(A549)_L2FC" : "L2F_A549", 
                                        "SARS-CoV-2(A549-ACE2)LowMOI_L2FC" : "L2F_A549_ACE2_LowMOI", 
                                        "SARS-CoV-2(A549-ACE2)HiMOI_L2FC" : "L2F_A549_ACE2_HiMOI"})
rna_data1 = rna_data1[rna_data1.PreySymbol.isin(ppi_data.PreySymbol)]

rna_calu3 = pd.DataFrame({"PreySymbol" : rna_data1.PreySymbol,
                          "L2F_Calu3" : rna_data1.L2F_Calu_3})
rna_calu3.insert(2, "L2F_abs", rna_calu3.L2F_Calu3.abs())
rna_calu3 = rna_calu3.sort_values(by = "L2F_abs", ascending = False)
rna_calu3.to_csv('processed/merged_ppi_lists/rankings/rna_calu3.csv',
                 index = False, na_rep = "NaN")

rna_A549 = pd.DataFrame({"PreySymbol" : rna_data1.PreySymbol,
                          "L2F_A549" : rna_data1.L2F_A549})
rna_A549.insert(2, "L2F_abs", rna_A549.L2F_A549.abs())
rna_A549 = rna_A549.sort_values(by = "L2F_abs", ascending = False)
rna_A549.to_csv('processed/merged_ppi_lists/rankings/rna_A549.csv',
                 index = False, na_rep = "NaN")

rna_A549_ACE2_LowMOI = pd.DataFrame({"PreySymbol" : rna_data1.PreySymbol,
                          "L2F_A549_ACE2_LowMOI" : rna_data1.L2F_A549_ACE2_LowMOI})
rna_A549_ACE2_LowMOI.insert(2, "L2F_abs", rna_A549_ACE2_LowMOI.L2F_A549_ACE2_LowMOI.abs())
rna_A549_ACE2_LowMOI = rna_A549_ACE2_LowMOI.sort_values(by = "L2F_abs", ascending = False)
rna_A549_ACE2_LowMOI.to_csv('processed/merged_ppi_lists/rankings/rna_A549_ACE2_LowMOI.csv',
                 index = False, na_rep = "NaN")

rna_A549_ACE2_HiMOI = pd.DataFrame({"PreySymbol" : rna_data1.PreySymbol,
                          "L2F_A549_ACE2_HiMOI" : rna_data1.L2F_A549_ACE2_HiMOI})
rna_A549_ACE2_HiMOI.insert(2, "L2F_abs", rna_A549_ACE2_HiMOI.L2F_A549_ACE2_HiMOI.abs())
rna_A549_ACE2_HiMOI = rna_A549_ACE2_HiMOI.sort_values(by = "L2F_abs", ascending = False)
rna_A549_ACE2_HiMOI.to_csv('processed/merged_ppi_lists/rankings/rna_A549_ACE2_HiMOI.csv',
                 index = False, na_rep = "NaN")

rna_data2 = pd.read_csv('RNASeq/blanco-mello/bm-mmc2.csv', 
                        usecols = ["GeneName", "SARS-CoV-2_L2FC"])
rna_data2 = rna_data2.rename(columns = {"GeneName" : "PreySymbol", 
                                        "SARS-CoV-2_L2FC" : "L2F_NHBE"})
rna_data2 = rna_data2[rna_data2.PreySymbol.isin(ppi_data.PreySymbol)]

rna_NHBE = pd.DataFrame({"PreySymbol" : rna_data2.PreySymbol,
                          "L2F_NHBE" : rna_data2.L2F_NHBE})
rna_NHBE.insert(2, "L2F_abs", rna_NHBE.L2F_NHBE.abs())
rna_NHBE = rna_NHBE.sort_values(by = "L2F_abs", ascending = False)
rna_NHBE.to_csv('processed/merged_ppi_lists/rankings/rna_NHBE.csv',
                 index = False, na_rep = "NaN")


rna_data3 = pd.read_csv('RNASeq/blanco-mello/bm-mmc3-ferret.csv', 
                        usecols = ["GeneSymbol", "SARS-CoV-2 Day7_L2FC"])
rna_data3 = rna_data3.rename(columns = {"GeneSymbol" : "PreySymbol", 
                                        "SARS-CoV-2 Day7_L2FC" : "L2F_D7_Ferret"})
rna_data3 = rna_data3[rna_data3.PreySymbol.isin(ppi_data.PreySymbol)]

rna_Ferret = pd.DataFrame({"PreySymbol" : rna_data3.PreySymbol,
                          "L2F_Ferret" : rna_data3.L2F_D7_Ferret})
rna_Ferret.insert(2, "L2F_abs", rna_Ferret.L2F_Ferret.abs())
rna_Ferret = rna_Ferret.sort_values(by = "L2F_abs", ascending = False)
rna_Ferret.to_csv('processed/merged_ppi_lists/rankings/rna_Ferret.csv',
                 index = False, na_rep = "NaN")

rna_data4 = pd.read_csv('RNASeq/blanco-mello/bm-mmc3-trachea.csv', 
                        usecols = ["GeneSymbol", "SARS-CoV-2 Day3_L2FC"])
rna_data4 = rna_data4.rename(columns = {"GeneSymbol" : "PreySymbol", 
                                        "SARS-CoV-2 Day3_L2FC" : "L2F_Trachea"})
rna_data4 = rna_data4[rna_data4.PreySymbol.isin(ppi_data.PreySymbol)]

rna_Trachea = pd.DataFrame({"PreySymbol" : rna_data4.PreySymbol,
                          "L2F_Trachea" : rna_data4.L2F_Trachea})
rna_Trachea.insert(2, "L2F_abs", rna_Trachea.L2F_Trachea.abs())
rna_Trachea = rna_Trachea.sort_values(by = "L2F_abs", ascending = False)
rna_Trachea.to_csv('processed/merged_ppi_lists/rankings/rna_Trachea.csv',
                 index = False, na_rep = "NaN")

rna_data5 = pd.read_csv('RNASeq/blanco-mello/bm-mmc4.csv', 
                        usecols = ["Gene_name", "log2FoldChange"])
rna_data5 = rna_data5.rename(columns = {"Gene_name" : "PreySymbol", 
                                        "log2FoldChange" : "L2F_Patient"})
rna_data5 = rna_data5[rna_data5.PreySymbol.isin(ppi_data.PreySymbol)]

rna_Patient = pd.DataFrame({"PreySymbol" : rna_data5.PreySymbol,
                          "L2F_Patient" : rna_data5.L2F_Patient})
rna_Patient.insert(2, "L2F_abs", rna_Patient.L2F_Patient.abs())
rna_Patient = rna_Patient.sort_values(by = "L2F_abs", ascending = False)
rna_Patient.to_csv('processed/merged_ppi_lists/rankings/rna_Patient.csv',
                 index = False, na_rep = "NaN")
'''

# Translateome time-dependent change
'''
trans_data = pd.read_csv('RNASeq/bojkova/translateome-time-data.csv', 
                         usecols = ["Gene Symbol01", "Ratio 2h", "Ratio 6h", "Ratio 10h", "Ratio 24h"],
                         na_values = ["#NUM!", "#DIV/0!"])
trans_data = trans_data.rename(columns = {"Gene Symbol01" : "PreySymbol", 
                                          "Ratio 2h" : "L2F_2h", 
                                          "Ratio 6h" : "L2F_6h", 
                                          "Ratio 10h" : "L2F_10h", 
                                          "Ratio 24h" : "L2F_24h"})
trans_data = trans_data[trans_data.PreySymbol.isin(ppi_data.PreySymbol)]

trans_range = (trans_data[["L2F_2h", "L2F_6h", "L2F_10h", "L2F_24h"]].max(axis = 1) - 
                trans_data[["L2F_2h", "L2F_6h", "L2F_10h", "L2F_24h"]].min(axis = 1))
trans_data.insert(5, "translateome_range", trans_range)
trans_data = trans_data.sort_values(by = "translateome_range", ascending = False)
trans_data.to_csv('processed/merged_ppi_lists/rankings/translateome.csv',
                  index = False, na_rep = "NaN")
'''

# Proteome time-dependent change
'''
prot_data = pd.read_csv('RNASeq/bojkova/proteome-time-data.csv', 
                         usecols = ["Gene Symbol", "Ratio 2h", "Ratio 6h", "Ratio 10h", "Ratio 24h"],
                         na_values = ["#NUM!", "#DIV/0!"])
prot_data = prot_data.rename(columns = {"Gene Symbol" : "PreySymbol", 
                                          "Ratio 2h" : "L2F_2h", 
                                          "Ratio 6h" : "L2F_6h", 
                                          "Ratio 10h" : "L2F_10h", 
                                          "Ratio 24h" : "L2F_24h"})
prot_data = prot_data[prot_data.PreySymbol.isin(ppi_data.PreySymbol)]

prot_range = (prot_data[["L2F_2h", "L2F_6h", "L2F_10h", "L2F_24h"]].max(axis = 1) - 
              prot_data[["L2F_2h", "L2F_6h", "L2F_10h", "L2F_24h"]].min(axis = 1))
prot_data.insert(5, "proteome_range", prot_range)
prot_data = prot_data.sort_values(by = "proteome_range", ascending = False)
prot_data.to_csv('processed/merged_ppi_lists/rankings/proteome.csv',
                  index = False, na_rep = "NaN")
'''

# Phosphorylation
'''
phospho_data = pd.read_csv('MMann/phosphoproteomics.csv',
                        skiprows = 1,
                        usecols = ["gene name", "fold_change_vs_mock_(log2)-6",
                                   "fold_change_vs_mock_(log2)-24"])
phospho_data = phospho_data.rename(columns = {"gene name" : "PreySymbol", 
                                              "fold_change_vs_mock_(log2)-6" : "L2F_6h", 
                                              "fold_change_vs_mock_(log2)-24" : "L2F_24h"})
phospho_data = phospho_data[phospho_data.PreySymbol.isin(ppi_data.PreySymbol)]

phospho_6h = pd.DataFrame({"PreySymbol" : phospho_data.PreySymbol,
                           "L2F_6h" : phospho_data.L2F_6h})
phospho_6h.insert(2, "L2F_abs", phospho_6h.L2F_6h.abs())
phospho_6h = phospho_6h.sort_values(by = "L2F_abs", ascending = False)
phospho_6h.to_csv('processed/merged_ppi_lists/rankings/phospho_6h.csv',
                 index = False, na_rep = "NaN")

phospho_24h = pd.DataFrame({"PreySymbol" : phospho_data.PreySymbol,
                           "L2F_24h" : phospho_data.L2F_24h})
phospho_24h.insert(2, "L2F_abs", phospho_24h.L2F_24h.abs())
phospho_24h = phospho_24h.sort_values(by = "L2F_abs", ascending = False)
phospho_24h.to_csv('processed/merged_ppi_lists/rankings/phospho_24h.csv',
                 index = False, na_rep = "NaN")
'''

# Write data to a file
for rank_file in ['counts', 'phospho_6h', 'phospho_24h', 'proteome', 'replicates', 'rna_A549',
                  'rna_A549_ACE2_HiMOI', 'rna_A549_ACE2_LowMOI', 'rna_calu3', 'rna_Ferret',
                  'rna_NHBE', 'rna_Patient', 'rna_Trachea', 'translateome']:
    file_path = 'processed/merged_ppi_lists/rankings/' + rank_file + '.csv'
    dataset = pd.read_csv(file_path)
    dataset = dataset[dataset.PreySymbol.notna()]
    dataset = dataset.drop_duplicates(subset = ["PreySymbol"])
    
    dataset_length = dataset.PreySymbol.count()
    dataset_cutoff = int(dataset_length / 20) + 1
    pct5_dataset = dataset[0:dataset_cutoff]
    
    output_path = 'processed/merged_ppi_lists/rankings/rankings_5pct/' + rank_file + '.csv'
    pct5_dataset.to_csv(output_path, index = False, na_rep = "NaN")

prot_list = ppi_data.PreySymbol.drop_duplicates()
prot_sig_dict = {}

# Create a total ranking row
for rank_file in ['counts', 'phospho_6h', 'phospho_24h', 'proteome', 'replicates', 'rna_A549',
                  'rna_A549_ACE2_HiMOI', 'rna_A549_ACE2_LowMOI', 'rna_calu3', 'rna_Ferret',
                  'rna_NHBE', 'rna_Patient', 'rna_Trachea', 'translateome']:
    input_file = 'processed/merged_ppi_lists/rankings/rankings_5pct/' + rank_file + '.csv'
    dataset = pd.read_csv(input_file)
    for protein in dataset.index:
        if dataset.loc[protein, "PreySymbol"] in prot_sig_dict:
            prot_sig_dict[dataset.loc[protein, "PreySymbol"]] += 1
        else:
            prot_sig_dict[dataset.loc[protein, "PreySymbol"]] = 1

def dict_apply(protein, source_dict):
    if protein in source_dict:
        return source_dict[protein]
    return 0

prot_sig_list = prot_list.apply(dict_apply, args = (prot_sig_dict,))
full_ranks = pd.DataFrame({"PreySymbol" : prot_list,
                           "SigCount" : prot_sig_list})
full_ranks = full_ranks.sort_values(by = "SigCount", ascending = False)
full_ranks.to_csv('processed/merged_ppi_lists/rankings/rankings_5pct/combined.csv',
                  index = False, na_rep = "NaN")
'''

# Uncomment to add an in_roth column
'''
ppis = pd.read_csv('processed/merged_ppi_lists/rankings/rankings_5pct/combined.csv')
roth = pd.read_csv('processed/merged_ppi_lists/Roth1-PPI.csv')
ppis = ppis[ppis.PreySymbol.isin(roth.PreySymbol)]
ppis.to_csv('processed/merged_ppi_lists/rankings/rankings_5pct/combined_roth.csv',
            index = False, na_rep = "NaN")
'''
"""


# Finding overlap between top 5% ranked datasets
"""
for rank_file1 in ['counts', 'phospho_6h', 'phospho_24h', 'proteome', 'replicates', 'rna_A549',
                  'rna_A549_ACE2_HiMOI', 'rna_A549_ACE2_LowMOI', 'rna_calu3', 'rna_Ferret',
                  'rna_NHBE', 'rna_Patient', 'rna_Trachea', 'translateome']:
    input_file1 = 'processed/merged_ppi_lists/rankings/rankings_5pct/' + rank_file1 + '.csv'
    dataset1 = pd.read_csv(input_file1)
    for rank_file2 in ['counts', 'phospho_6h', 'phospho_24h', 'proteome', 'replicates', 'rna_A549',
                       'rna_A549_ACE2_HiMOI', 'rna_A549_ACE2_LowMOI', 'rna_calu3', 'rna_Ferret',
                       'rna_NHBE', 'rna_Patient', 'rna_Trachea', 'translateome']:
        input_file2 = 'processed/merged_ppi_lists/rankings/rankings_5pct/' + rank_file2 + '.csv'
        dataset2 = pd.read_csv(input_file2)
        shared_set = dataset1[dataset1.PreySymbol.isin(dataset2.PreySymbol)]
        shared_pct = shared_set.PreySymbol.count() / dataset1.PreySymbol.count()
        print(rank_file1 + " has " + str(shared_pct) + " overlap with " + rank_file2)
"""

# Finding proteins in a significant GO subset which are also ranked
"""
ppis = pd.read_csv('processed/merged_ppi_lists/rankings/rankings_5pct/combined.csv')
tgf_beta = pd.read_table('Human_PPIs/Sig_Sets/tgf_beta_signalling.txt')
tgf_preys = tgf_beta.GeneSymbol.drop_duplicates()
ppis = ppis[ppis.SigCount > 0]
rel_ppis = ppis[ppis.PreySymbol.isin(tgf_preys)]
print(rel_ppis.PreySymbol.to_list())
print(rel_ppis.PreySymbol.count() / ppis.PreySymbol.count())
"""

# Adding and ranking GO tags for significantly ranked proteins
"""
ppis = pd.read_csv('processed/merged_ppi_lists/rankings/rankings_5pct/combined.csv')
roth_ppis = pd.read_csv('processed/merged_ppi_lists/rankings/rankings_5pct/combined_roth.csv')
ppis[ppis.SigCount > 0][["PreySymbol"]].to_csv('processed/merged_ppi_lists/rankings/rankings_5pct/combined_GO.csv',
                            index = False)
roth_ppis[roth_ppis.SigCount > 0][["PreySymbol"]].to_csv('processed/merged_ppi_lists/rankings/rankings_5pct/combined_roth_GO.csv',
                            index = False)
"""

"""
combined = pd.read_table('processed/merged_ppi_lists/rankings/rankings_5pct/combined_GO.tab',
                         usecols = ["Gene_names_primary", "GO_IDs"])
combined_prot = pd.read_csv('processed/merged_ppi_lists/rankings/rankings_5pct/combined.csv')

combined = combined[combined.Gene_names_primary.isin(combined_prot.PreySymbol)]
combined = combined.drop_duplicates()

combined_GO_dict = {}

def split_GO(GOs):
    GO_list = str(GOs).split("; ")
    for GO in GO_list:
        if GO in combined_GO_dict:
            combined_GO_dict[GO] += 1
        else:
            combined_GO_dict[GO] = 1

combined.GO_IDs = combined.GO_IDs.apply(split_GO)
combined_GO_df = pd.DataFrame()
combined_GO_df = combined_GO_df.from_dict(combined_GO_dict, orient = 'index')
combined_GO_df = combined_GO_df.sort_values(by = 0, ascending = False)
combined_GO_df.to_csv('processed/merged_ppi_lists/rankings/rankings_5pct/combined_GO.csv')


combined_roth = pd.read_table('processed/merged_ppi_lists/rankings/rankings_5pct/combined_roth_GO.tab',
                         usecols = ["Gene_names_primary", "GO_IDs"])
combined_roth_prot = pd.read_csv('processed/merged_ppi_lists/rankings/rankings_5pct/combined_roth.csv')

combined_roth = combined_roth[combined_roth.Gene_names_primary.isin(combined_roth_prot.PreySymbol)]
combined_roth = combined_roth.drop_duplicates()

combined_roth_GO_dict = {}

def split_GO(GOs):
    GO_list = str(GOs).split("; ")
    for GO in GO_list:
        if GO in combined_roth_GO_dict:
            combined_roth_GO_dict[GO] += 1
        else:
            combined_roth_GO_dict[GO] = 1

combined_roth.GO_IDs = combined_roth.GO_IDs.apply(split_GO)
combined_roth_GO_df = pd.DataFrame()
combined_roth_GO_df = combined_roth_GO_df.from_dict(combined_roth_GO_dict, orient = 'index')
combined_roth_GO_df = combined_roth_GO_df.sort_values(by = 0, ascending = False)
combined_roth_GO_df.to_csv('processed/merged_ppi_lists/rankings/rankings_5pct/combined_roth_GO.csv')
"""


# Creating a combined ranking file of all RNA data
"""
proteins = pd.read_csv('processed/merged_ppi_lists/rankings/rankings_5pct/combined.csv')

for sig_set in ['counts', 'rna_A549', 
                'rna_A549_ACE2_HiMOI', 'rna_A549_ACE2_LowMOI', 'rna_calu3', 'rna_Ferret',
                'rna_NHBE', 'rna_Patient', 'rna_Trachea']:
    sig_file = 'processed/merged_ppi_lists/rankings/' + sig_set + '.csv'
    sig_data = pd.read_csv(sig_file)
    proteins = proteins.join(sig_data.set_index("PreySymbol"), on = "PreySymbol")

proteins.to_csv('processed/merged_ppi_lists/rankings/all_rankings.csv',
                index = False, na_rep = "NaN")    
"""
"""
all_rankings = pd.read_csv('processed/merged_ppi_lists/rankings/all_rankings.csv',
                           usecols = ["PreySymbol", "SigCount", "PPICount",
                                      "L2F_A549", "L2F_A549_ACE2_HiMOI", "L2F_A549_ACE2_LowMOI",
                                      "L2F_Calu3", "L2F_Ferret", "L2F_NHBE", "L2F_Patient",
                                      "L2F_Trachea"])
all_rna = all_rankings[["L2F_A549", "L2F_A549_ACE2_HiMOI", "L2F_A549_ACE2_LowMOI", "L2F_Calu3", 
                        "L2F_Ferret", "L2F_NHBE", "L2F_Patient", "L2F_Trachea"]]

concensus = all_rna.apply(lambda x : (sum(x > 0) - (x.count() / 2)) / (x.count() / 2), axis = 1)
all_rankings.insert(3, "L2F_concensus", concensus)
all_rankings.to_csv('processed/merged_ppi_lists/rankings/all_rankings.csv',
                    index = False, na_rep = "NaN")
"""
