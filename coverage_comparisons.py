import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Uncomment to run a coverage analysis of various datasets relative 
# to the hORFEOME dataset. 
# Print total Prey proteins involved, total PPIs listed, and coverage of Prey / hORFEOME
"""
hORFEOME = 'Human_PPIs/HuRI/Supplementary_Tables_new/Supplementary Table 2.txt'
hORF = pd.read_csv(hORFEOME, delim_whitespace = True, usecols = ["symbol"])
hORF = hORF.drop_duplicates()
hORF = hORF.rename(columns = {"symbol" : "PreySymbol"})

def print_coverage(input_file, prey_col, dataset_name):
    dataset = pd.read_csv(input_file, usecols = [prey_col])
    dataset = dataset.rename(columns = {prey_col : "PreySymbol"})
    ppi_count = dataset.PreySymbol.count()
    dataset = dataset.drop_duplicates()
    overlap_set = dataset[dataset.PreySymbol.isin(hORF.PreySymbol)]
    prey_count = overlap_set.PreySymbol.count()
    coverage = overlap_set.PreySymbol.count() / hORF.PreySymbol.count() * 100
    print(dataset_name + ": " + str(prey_count) + " (" + str(coverage) + "%) from " + str(ppi_count))

print_coverage('processed/merged_ppi_lists/Individual_Datasets/DGordon-High-Confidence-Pairs.csv', 
               "PreyGene", "DGordon")
print_coverage('processed/merged_ppi_lists/Individual_Datasets/Li-SARS-2-PPI-networks.csv', 
               "genename", "Li")
print_coverage('processed/merged_ppi_lists/Individual_Datasets/MMann-PPI-M2si.csv', 
               "gene_name", "MMann")
print_coverage('AP-MS/etienne/Etienne-High-Conf.csv', "Gene_name", "Etienne")
print_coverage('processed/merged_ppi_lists/Individual_Datasets/Roth1-PPI-GoFiltered.csv', 
               "PreySymbol", "Roth")
"""


# ALl code below generates a histogram of PPI distributions for each dataset
def series_list_index(l, index):
    ''' (List, index) -> List
    
    Returns the list item at the given index. Only exists to use on a pandas series of lists.
    There's probably a better way to do this, but this'll do.
    '''
    return l[index]

# Process the DGordon dataset to a better form
dg = pd.read_csv('processed/merged_ppi_lists/Individual_Datasets/DGordon-High-Confidence-Pairs.csv',
                 usecols = ["Bait", "PreyGene"])
dg = dg.rename(columns = {"Bait" : "BaitSymbol", "PreyGene" : "PreySymbol"})
dg.BaitSymbol = dg.BaitSymbol.str.split(" ").apply(series_list_index, args = (1,))
dg.BaitSymbol[dg.BaitSymbol == "Spike"] = "S"
dg.BaitSymbol = dg.BaitSymbol.str.upper()

# Process the Li dataset to a better form
li = pd.read_csv('processed/merged_ppi_lists/Individual_Datasets/Li-SARS-2-PPI-networks.csv',
                 usecols = ["protein", "genename"])
li = li.rename(columns = {"protein" : "BaitSymbol", "genename" : "PreySymbol"})
li = li.drop(li[li.BaitSymbol.str.startswith("nsp3")].index)
li.BaitSymbol[li.BaitSymbol == "spike"] = "S"
li.BaitSymbol[li.BaitSymbol == "envelope"] = "E"
li.BaitSymbol[li.BaitSymbol == "membrane"] = "M"
li.BaitSymbol[li.BaitSymbol == "nucleocapsid"] = "N"
li.BaitSymbol = li.BaitSymbol.str.upper()

# Process MMann dataset to a better form
MM = pd.read_csv('processed/merged_ppi_lists/Individual_Datasets/MMann-PPI-M2si.csv',
                 usecols = ["bait_full_name", "gene_name"])
MM = MM.rename(columns = {"bait_full_name" : "BaitSymbol", "gene_name" : "PreySymbol"})
MM = MM[MM.BaitSymbol.str.startswith("SARS_CoV2")]
MM.BaitSymbol = MM.BaitSymbol.str.split("_").apply(series_list_index, args = (2,))
MM.PreySymbol = MM.PreySymbol.str.split(";").apply(series_list_index, args = (0,))
MM.BaitSymbol = MM.BaitSymbol.str.upper()
MM = MM.drop(MM[MM.BaitSymbol == 'NSP3_MACROD'].index)

# Process Etienne dataset to a better form
et = pd.read_csv('AP-MS/etienne/Etienne-High-Conf.csv', 
                 usecols = ["Bait", "Gene_name"])
et = et.rename(columns = {"Bait" : "BaitSymbol", "Gene_name" : "PreySymbol"})
et.BaitSymbol = et.BaitSymbol.str.upper()

# Process Roth dataset to a better form
fr = pd.read_csv('processed/merged_ppi_lists/Individual_Datasets/Roth1-PPI-GoFiltered.csv',
                 usecols = ["BaitSymbol", "PreySymbol"])
fr.BaitSymbol = fr.BaitSymbol.str.upper()
fr = fr.drop(fr[fr.BaitSymbol == "ACE2"].index)


# Process Prey data
combined_preys = dg.PreySymbol.append(li.PreySymbol).append(MM.PreySymbol).append(et.PreySymbol).append(fr.PreySymbol).drop_duplicates()
combined_baits = dg.BaitSymbol.append(li.BaitSymbol).append(MM.BaitSymbol).append(et.BaitSymbol).append(fr.BaitSymbol).drop_duplicates()

dg_prey_count = []
li_prey_count = []
MM_prey_count = []
et_prey_count = []
fr_prey_count = []
dg_bait_count = []
li_bait_count = []
MM_bait_count = []
et_bait_count = []
fr_bait_count = []


def count_prey(gene):
    dg_prey_count.append(dg[dg.PreySymbol == gene].PreySymbol.count())
    li_prey_count.append(li[li.PreySymbol == gene].PreySymbol.count())
    MM_prey_count.append(MM[MM.PreySymbol == gene].PreySymbol.count())
    et_prey_count.append(et[et.PreySymbol == gene].PreySymbol.count())
    fr_prey_count.append(fr[fr.PreySymbol == gene].PreySymbol.count())
    return gene

def conut_bait(gene):
    dg_bait_count.append(dg[dg.BaitSymbol == gene].BaitSymbol.count())
    li_bait_count.append(li[li.BaitSymbol == gene].BaitSymbol.count())
    MM_bait_count.append(MM[MM.BaitSymbol == gene].BaitSymbol.count())
    et_bait_count.append(et[et.BaitSymbol == gene].BaitSymbol.count())
    fr_bait_count.append(fr[fr.BaitSymbol == gene].BaitSymbol.count())
    return gene

combined_preys = combined_preys.apply(count_prey)


prey_counts = pd.DataFrame({"PreySymbol" : combined_preys, "dg" : dg_prey_count, "li" : li_prey_count,
                            "MM" : MM_prey_count, "et" : et_prey_count, "fr" : fr_prey_count})
prey_counts.insert(5, "Past", (prey_counts.dg + prey_counts.li + prey_counts.MM + prey_counts.et))


prey_freq = pd.DataFrame({"PreySymbol" : prey_counts.PreySymbol, 
                          "Gordon et al" : (prey_counts.fr/prey_counts.fr.sum()) / (prey_counts.dg / prey_counts.dg.sum()),
                          "Li et al" : (prey_counts.fr/prey_counts.fr.sum()) / (prey_counts.li / prey_counts.li.sum()),
                          "Mann et al" : (prey_counts.fr/prey_counts.fr.sum()) / (prey_counts.MM / prey_counts.MM.sum()),
                          "Etienne et al" : (prey_counts.fr/prey_counts.fr.sum()) / (prey_counts.et / prey_counts.MM.sum()),
                          "Past" : (prey_counts.fr/prey_counts.fr.sum()) / (prey_counts.Past / prey_counts.Past.sum())})
prey_freq = prey_freq.set_index("PreySymbol")
prey_freq[prey_freq.isna()] = 0

for column in prey_freq.columns:
    prey_freq[column][prey_freq[column] == np.inf] = -1
    prey_freq[column][prey_freq[column] == -1] = prey_freq[column].max()

# Uncomment to perform bait analysis
'''
combined_baits = combined_baits.apply(count_bait)
bait_counts = pd.DataFrame({"BaitSymbol" : combined_baits, "dg" : dg_bait_count, "li" : li_bait_count,
                            "MM" : MM_bait_count, "et" : et_bait_count, "fr" : fr_bait_count})
bait_counts.insert(5, "Past", (bait_counts.dg + bait_counts.li + bait_counts.MM + bait_counts.et))
bait_counts = bait_counts.set_index('BaitSymbol')
bait_counts = bait_counts.rename(columns = {"dg" : "Gordon et al", "li" : "Li et al",
                                            "MM" : "Mann et al", "et" : "Etienne et al",
                                            "fr" : "Roth Lab"})


bait_freq = pd.DataFrame({"BaitSymbol" : bait_counts.BaitSymbol, 
                          "Gordon et al" : (bait_counts.fr/bait_counts.fr.sum()) / (bait_counts.dg / bait_counts.dg.sum()),
                          "Li et al" : (bait_counts.fr/bait_counts.fr.sum()) / (bait_counts.li / bait_counts.li.sum()),
                          "Mann et al" : (bait_counts.fr/bait_counts.fr.sum()) / (bait_counts.MM / bait_counts.MM.sum()),
                          "Etienne et al" : (bait_counts.fr/bait_counts.fr.sum()) / (bait_counts.et / bait_counts.MM.sum()),
                          "Past" : (bait_counts.fr/bait_counts.fr.sum()) / (bait_counts.Past / bait_counts.Past.sum())})
bait_freq = bait_freq.set_index("BaitSymbol")
bait_freq[bait_freq.isna()] = 0

for column in bait_freq.columns:
    bait_freq[column][bait_freq[column] == np.inf] = -1
    bait_freq[column][bait_freq[column] == -1] = bait_freq[column].max()

bait_counts.plot.bar(rot = 90, subplots = True, figsize = (6, 5), legend = False)
# prey_counts.plot.bar(rot = 90, subplots = True, figsize = (6, 5), legend = False)
# bait_freq.plot.bar(rot = 90, subplots = True, figsize = (6,5), legend = False)
'''
prey_freq.plot.bar(rot = 0, subplots = True, figsize = (6,5), legend = False)
plt.tight_layout()
plt.show()



# Uncomment to generate GO Frequency Plots
"""
dg = pd.read_table('processed/merged_ppi_lists/Individual_Datasets/DGordon-GO.tab', 
                 usecols = ["PreySymbol", "GO_IDs"])
li = pd.read_table('processed/merged_ppi_lists/Individual_Datasets/Li-GO.tab', 
                 usecols = ["PreySymbol", "GO_IDs"])
MM = pd.read_table('processed/merged_ppi_lists/Individual_Datasets/MMann-Go.tab', 
                 usecols = ["PreySymbol", "GO_IDs"])
et = pd.read_table('AP-MS/etienne/Etienne-GO.tab', 
                 usecols = ["PreySymbol", "GO_IDs"])
fr = pd.read_table('processed/merged_ppi_lists/Individual_Datasets/Roth1-GO.tab', 
                 usecols = ["PreySymbol", "GO_IDs"])

dg_GO = []
li_GO = []
MM_GO = []
et_GO = []
fr_GO = []

def semicolon_parse(item, go_list):
    go_list.extend(str(item).split('; '))
    return item

dg.GO_IDs = dg.GO_IDs.apply(semicolon_parse, args = (dg_GO,))
li.GO_IDs = li.GO_IDs.apply(semicolon_parse, args = (li_GO,))
MM.GO_IDs = MM.GO_IDs.apply(semicolon_parse, args = (MM_GO,))
et.GO_IDs = et.GO_IDs.apply(semicolon_parse, args = (et_GO,))
fr.GO_IDs = fr.GO_IDs.apply(semicolon_parse, args = (fr_GO,))

dg_GO = pd.Series(dg_GO)
li_GO = pd.Series(li_GO)
MM_GO = pd.Series(MM_GO)
et_GO = pd.Series(et_GO)
fr_GO = pd.Series(fr_GO)

all_GO = dg_GO.append(li_GO).append(MM_GO).append(et_GO).append(fr_GO).drop_duplicates()

dg_GO_count = []
li_GO_count = []
MM_GO_count = []
et_GO_count = []
fr_GO_count = []

def count_GO(item):
    dg_GO_count.append(dg_GO[dg_GO == item].count())
    li_GO_count.append(li_GO[li_GO == item].count())
    MM_GO_count.append(MM_GO[MM_GO == item].count())
    et_GO_count.append(et_GO[et_GO == item].count())
    fr_GO_count.append(fr_GO[fr_GO == item].count())
    return item

all_GO = all_GO.apply(count_GO)

GO_counts = pd.DataFrame({"PreySymbol" : all_GO, "dg" : dg_GO_count, "li" : li_GO_count,
                          "MM" : MM_GO_count, "et" : et_GO_count, "fr" : fr_GO_count})
GO_counts.insert(5, "Past", (GO_counts.dg + GO_counts.li + GO_counts.MM + GO_counts.et))

# GO_counts.plot.bar(rot = 0, subplots = True, figsize = (6,5), legend = False, 
#                    ylim = [0, 5])

GO_freq = pd.DataFrame({"PreySymbol" : GO_counts.PreySymbol,
                        "dg" : (GO_counts.fr / GO_counts.fr.sum()) / (GO_counts.dg / GO_counts.dg.sum()),
                        "li" : (GO_counts.fr / GO_counts.fr.sum()) / (GO_counts.li / GO_counts.li.sum()),
                        "MM" : (GO_counts.fr / GO_counts.fr.sum()) / (GO_counts.MM / GO_counts.MM.sum()),
                        "et" : (GO_counts.fr / GO_counts.fr.sum()) / (GO_counts.et / GO_counts.et.sum()),
                        "Past" : (GO_counts.fr / GO_counts.fr.sum()) / (GO_counts.Past / GO_counts.Past.sum())})
GO_freq = GO_freq.set_index("PreySymbol")

# print(GO_freq[GO_freq.Past == np.inf])
# GO_freq = GO_freq.sort_values(by = ["Past"])
# print(GO_freq)

print(GO_freq[GO_freq == np.inf])


# GO_freq.plot.bar(rot = 0, subplots = True, figsize = (6,5), legend = False)
# plt.tight_layout()
# plt.show()
"""