import pandas as pd

# This file contains a set of operations to merge datasets. Each operation is
# its own block with a short explanation. Uncomment the desired operations to
# run.

# Merges the PPIs in the D. Gordon, Li, and M. Mann Papers.
"""
# Should have used a for loop to limit memory usage
dgordon = pd.read_csv("AP-MS/PPI-Output/DGordon-High-Confidence-Pairs.csv") 
li = pd.read_csv("AP-MS/PPI-Output/Li-SARS-2-PPI-networks.csv")
MMann = pd.read_csv("AP-MS/PPI-Output/MMann-PPI-M2si.csv")

dgordon = dgordon[["Bait", "Preys", "PreyGene"]]
dgordon = dgordon.rename(columns = {"Bait" : "Bait", "Preys" : "Prey", "PreyGene" : "PreyGene"})
li = li[["protein", "accession", "genename"]]
li = li.rename(columns = {"protein" : "Bait", "accession" : "Prey", "genename" : "PreyGene"})
MMann = MMann[["bait_full_name", "gene_name", "majority_protein_acs"]]
MMann[["majority_protein_acs", "gene_name"]] = MMann[["gene_name", "majority_protein_acs"]]
MMann = MMann.rename(columns = {"bait_full_name" : "Bait", "gene_name" : "Prey", "majority_protein_acs" : "PreyGene"})

MMann = MMann[MMann.Bait.str.startswith("SARS_CoV2")]  # MMann went above and beyond, but we don't need it

def series_list_index(l, index):
    ''' (List, index) -> Item in list
    
    Returns the list item at the given index. Only exists to use on a pandas series of lists.
    There's definitely a better way to do this, but I don't know what it is.
    '''
    return l[index]


# Standardizing SARS-CoV2 Gene Names
dgordon.Bait = dgordon.Bait.str.split(" ").apply(series_list_index, args = (1,))
dgordon.Bait[dgordon.Bait == "Spike"] = "S"
dgordon.Bait = dgordon.Bait.str.upper()

li.Bait[li.Bait == "spike"] = "S"
li.Bait[li.Bait == "envelope"] = "E"
li.Bait[li.Bait == "membrane"] = "M"
li.Bait[li.Bait == "nucleocapsid"] = "N"
li.Bait = li.Bait.str.upper()

MMann.Bait = MMann.Bait.str.split("_").apply(series_list_index, args = (2,))
MMann.Prey = MMann.Prey.str.split(";").apply(series_list_index, args = (0,))
MMann.Bait = MMann.Bait.str.upper()

full_list = pd.DataFrame()
full_list = full_list.append(dgordon)
full_list = full_list.append(li)
full_list = full_list.append(MMann)
full_list = full_list.drop_duplicates()
full_list.insert(1, "BaitIsCoV", [1] * len(full_list.Bait))
# The Below does not work. I manually did what it's supposed to do, fix if this ever gets run again.
full_list.insert(3, "PreyIsCoV", full_list.Prey.isin(full_list.Bait).astype(int))

full_list.to_csv("AP-MS/PPI-Output/Merged-List.csv", index = False)
"""

# Merges all Etienne files into one DataFrame. Adds this to a pre-existing PPI file. 
"""
def etienne_file_dir(protein : str):
    ''' (str) -> str
    
    Concatenates the filename of one of Etienne's files to be the file's path.
    
    >>> etienne_file_dir('nsp5')
    'AP-MS/etienne/Etienne-nsp5.csv'
    '''
    
    return "AP-MS/etienne/Etienne-" + protein + ".csv"

def under_value(pd_series, value):
    ''' (pandas series, int) -> list
    
    Given a series and a cutoff, returns the indices of values in the series which are below the cutoff. 
    Use to remove values below a threshold from a DataFrame.
    
    >>> under_value(pd.Series([1, 2, 3, 4, 5]), 3)
    [1, 2]
    '''
    
    below_value = []
    for index, item in enumerate(pd_series):
        if item.isdigit():
            if int(item) < value:
                below_value.append(index)
        else:
            below_value.append(index)
    return below_value

etienne_files = ["e", "m", "n", "nsp1", "nsp2", "nsp4", "nsp5", "nsp6", "nsp7", "nsp8", "nsp9", 
"nsp10", "nsp12", "nsp13", "nsp14", "nsp15", "nsp16", "orf3a", "orf3b", "orf6", "orf7a", "orf7b", 
"orf8", "orf9b", "plpro", "s"]

etienne_ppis = pd.DataFrame()
for etienne_file in etienne_files:
    temp_ppi_data = pd.read_csv(etienne_file_dir(etienne_file), skiprows = 2)
    temp_ppi_data = temp_ppi_data[["Uniprot ID", "Name", "# baits"]]
    temp_ppi_data = temp_ppi_data.rename(columns = {"Uniprot ID" : "Prey", "Name" : "PreyGene", "# baits" : "baitsCount"})
    temp_ppi_data = temp_ppi_data.drop(under_value(temp_ppi_data.baitsCount, 2))  # Remove this line to include all values
    temp_ppi_data = temp_ppi_data.drop(columns = ['baitsCount'])
    temp_ppi_data.insert(0, "Bait", [etienne_file.upper()] * len(temp_ppi_data.Prey))
    temp_ppi_data.insert(1, "BaitIsCoV", [1] * len(temp_ppi_data.Prey))
    temp_ppi_data.insert(3, "PreyIsCoV", [0] * len(temp_ppi_data.Prey))
    etienne_ppis = etienne_ppis.append(temp_ppi_data)

merged_ppis = pd.read_csv("AP-MS/PPI-Output/Merged-PPI-List-wo-etienne.csv")
merged_ppis = merged_ppis.append(etienne_ppis)
merged_ppis = merged_ppis.drop_duplicates()

merged_ppis.to_csv("AP-MS/PPI-Output/Merged-PPI-List-filtered.csv", index = False)
"""

# Merging datasets to include a source column. Duplicates not removed.
"""
def series_list_index(l, index):
    ''' (List, index) -> Item in list
    
    Returns the list item at the given index. Only exists to use on a pandas series of lists.
    There's definitely a better way to do this, but I don't know what it is.
    '''
    return l[index]

all_ppis = pd.DataFrame()
CoV_Baits = ["NSP1", "NSP2", "NSP3", "NSP4", "NSP5", "NSP6", "NSP7", "NSP8", "NSP8", "NSP9", 
             "NSP10", "NSP11", "NSP12", "NSP13", "NSP14", "NSP15", "NSP16", "S", "ORF3A", "ORF3B",
             "E", "M", "ORF6", "ORF7A", "ORF7B", "ORF8", "N", "ORF9A", "ORF9B", "ORF10"]

dg = pd.read_csv('processed/merged_ppi_lists/Individual_Datasets/DGordon-High-Confidence-Pairs.csv',
                 usecols = ["Bait", "PreyGene"])
dg = dg.rename(columns = {"Bait" : "BaitSymbol", "PreyGene" : "PreySymbol"})
dg.BaitSymbol = dg.BaitSymbol.str.split(" ").apply(series_list_index, args = (1,))
dg.BaitSymbol[dg.BaitSymbol == "Spike"] = "S"
dg.BaitSymbol = dg.BaitSymbol.str.upper()
dg.insert(2, "PPISymbol", dg.BaitSymbol + '_' + dg.PreySymbol)
dg.insert(3, "Source", "Gordon et al. (2020)")

li = pd.read_csv('processed/merged_ppi_lists/Individual_Datasets/Li-SARS-2-PPI-networks.csv',
                 usecols = ["protein", "genename"])
li = li.rename(columns = {"protein" : "BaitSymbol", "genename" : "PreySymbol"})
li = li.drop(li[li.BaitSymbol.str.startswith("nsp3")].index)
li.BaitSymbol[li.BaitSymbol == "spike"] = "S"
li.BaitSymbol[li.BaitSymbol == "envelope"] = "E"
li.BaitSymbol[li.BaitSymbol == "membrane"] = "M"
li.BaitSymbol[li.BaitSymbol == "nucleocapsid"] = "N"
li.BaitSymbol = li.BaitSymbol.str.upper()
li.insert(2, "PPISymbol", li.BaitSymbol + '_' + li.PreySymbol)
li.insert(3, "Source", "Li et al. (2020)")

MM = pd.read_csv('processed/merged_ppi_lists/Individual_Datasets/MMann-PPI-M2si.csv',
                 usecols = ["bait_full_name", "gene_name"])
MM = MM.rename(columns = {"bait_full_name" : "BaitSymbol", "gene_name" : "PreySymbol"})
MM = MM[MM.BaitSymbol.str.startswith("SARS_CoV2")]
MM.BaitSymbol = MM.BaitSymbol.str.split("_").apply(series_list_index, args = (2,))
MM.PreySymbol = MM.PreySymbol.str.split(";").apply(series_list_index, args = (0,))
MM.BaitSymbol = MM.BaitSymbol.str.upper()
MM = MM.drop(MM[MM.BaitSymbol == 'NSP3_MACROD'].index)
MM.BaitSymbol[MM.BaitSymbol == 'ORF3'] = 'ORF3A'
MM.insert(2, "PPISymbol", MM.BaitSymbol + '_' + MM.PreySymbol)
MM.insert(3, "Source", "Stukalov et al. (2020)")

et = pd.read_csv('AP-MS/etienne/Etienne-High-Conf.csv', 
                 usecols = ["Bait", "Gene_name"])
et = et.rename(columns = {"Bait" : "BaitSymbol", "Gene_name" : "PreySymbol"})
et.BaitSymbol = et.BaitSymbol.str.upper()
et.BaitSymbol[et.BaitSymbol == 'ORF9B'] = 'ORF9A'
et.BaitSymbol[et.BaitSymbol == 'ORF14'] = 'ORF9B'
et.insert(2, "PPISymbol", et.BaitSymbol + '_' + et.PreySymbol)
et.insert(3, "Source", "Etienne et al. (2020)")

ash = pd.read_csv('AP-MS/Ash_Combined/SARS_MERS_SARS2 Known PPIs Combined_PPIs.csv',
                  usecols = ["Virus", "Host_gene_symbol", "Virus_gene", 
                             "Source"])[["Virus", "Virus_gene", "Host_gene_symbol", "Source"]]
ash = ash[ash.Virus == "SARS-CoV-2"]
ash = ash.drop(columns = ["Virus"])
ash = ash.rename(columns = {"Virus_gene" : "BaitSymbol", "Host_gene_symbol" : "PreySymbol", 
                            "Source" : "Source"})
ash.BaitSymbol = ash.BaitSymbol.str.upper()
ash = ash[ash.BaitSymbol.isin(CoV_Baits)]
ash.insert(2, "PPISymbol", ash.BaitSymbol + '_' + ash.PreySymbol)

all_ppis = all_ppis.append(dg).append(li).append(MM).append(et).append(ash)
all_ppis.insert(1, "BaitIsCoV", 1)

all_ppis.to_csv('processed/merged_ppi_lists/Past-CoV-PPI-Sourced.csv', 
                index = False, na_rep = "NaN")
"""

# Generating the Roth1 + Past sourced dataset
"""
past = pd.read_csv('processed/merged_ppi_lists/Past-CoV-PPI-Sourced.csv')
roth = pd.read_csv('processed/merged_ppi_lists/Roth1-PPI.csv')
roth.insert(4, "Source", "Roth Lab")
ppis = past.append(roth)
ppis.to_csv('processed/merged_ppi_lists/Roth1-Past-PPI-Sourced.csv')
"""