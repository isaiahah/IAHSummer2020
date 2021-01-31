import pandas as pd

# This file contains a set of operations to subset datasets and attempt to find important subsets.
# It is largely outdated by protein_ranking.py. 
# Each operation is its own block with a short explanation. 
# Uncomment the desired operations to run.

# Outdated subsetting with poor file lists using non-GO lists and producing subset files. 
# DO NOT USE.
# Left for reference.
'''
raise Exception
full_gene_list = []

def protein_subset(subset_file, ppis_file, review_filter, output_file):
    subset = pd.read_table(subset_file, usecols = ["Status", "Gene_names", "Gene_names_(primary)"])
    if review_filter:
        subset = subset[subset.Status == "reviewed"]
    subset = subset.drop(columns = "Status")
    subset = subset.rename(columns = {"Gene_names" : "AllGenes",
                              "Gene_names_(primary)" : "MainGene"})

    ppis = pd.read_csv(ppis_file)
    
    subset.AllGenes = subset.AllGenes.apply(series_concat)

    rel_ppis = ppis[ppis.PreySymbol.isin(full_gene_list)]
    rel_ppis = rel_ppis.append(ppis[ppis.BaitSymbol.isin(full_gene_list)])
    rel_ppis = rel_ppis.drop_duplicates(subset = ["BaitSymbol", "PreySymbol"])
    rel_ppis.to_csv(output_file, index = False, na_rep = "NaN")
    
    
def strict_protein_subset(subset_file, ppis_file, review_filter, output_file):
    subset = pd.read_table(subset_file, usecols = ["Status", "Gene_names", "Gene_names_(primary)"])
    if review_filter:
        subset = subset[subset.Status == "reviewed"]
    subset = subset.drop(columns = "Status")
    subset = subset.rename(columns = {"Gene_names" : "AllGenes",
                              "Gene_names_(primary)" : "MainGene"})

    ppis = pd.read_csv(ppis_file)
    
    subset.AllGenes = subset.AllGenes.apply(series_concat)

    rel_ppis1 = ppis[ppis.PreySymbol.isin(full_gene_list)]
    rel_ppis2 = ppis[ppis.BaitSymbol.isin(full_gene_list)]
    rel_ppis2 = rel_ppis2[rel_ppis2.PreySymbol.isin(full_gene_list)]
    rel_ppis = rel_ppis1.append(rel_ppis2)
    rel_ppis = rel_ppis.drop_duplicates(subset = ["BaitSymbol", "PreySymbol"])

    rel_ppis.to_csv(output_file, index = False, na_rep = "NaN")
    
    
def series_concat(items):
    for item in str(items).split():
        full_gene_list.append(item)
    return items

# Generating the Ubiquitation-related subset of proteins.
# Uncomment to run.
"""
protein_subset('Human_PPIs/Sig_Sets/ubiquitation.tab', 
               'processed/merged_ppi_lists/ppi-rna/RNA-Roth1-Past-Human-PPI.csv',
               False,
               'processed/merged_ppi_lists/subsets/ubiquitin_all.csv')

strict_protein_subset('Human_PPIs/Sig_Sets/ubiquitation.tab', 
               'processed/merged_ppi_lists/ppi-rna/RNA-Roth1-Past-Human-PPI.csv',
               False,
               'processed/merged_ppi_lists/subsets/ubiquitin_only.csv')
"""

# Generating the Inflammation-related subset of proteins
# Uncomment to run.
"""
protein_subset('Human_PPIs/Sig_Sets/inflammation.tab', 
               'processed/merged_ppi_lists/ppi-rna/RNA-Roth1-Past-Human-PPI.csv',
               False,
               'processed/merged_ppi_lists/subsets/inflammation_all.csv')

strict_protein_subset('Human_PPIs/Sig_Sets/inflammation.tab', 
               'processed/merged_ppi_lists/ppi-rna/RNA-Roth1-Past-Human-PPI.csv',
               False,
               'processed/merged_ppi_lists/subsets/inflammation_only.csv')
"""

# Generating the Mediator Complex-related subset of proteins.
# Uncomment to run.
"""
ppis = 'processed/merged_ppi_lists/ppi-rna/RNA-Roth1-Past-Human-PPI.csv'
ppis = pd.read_csv(ppis)
mediator = ["MED1", "MED4", "MED6", "MED7", "MED8", "MED9", "MED10", "MED11", "MED12", "MED13", 
            "MED13L", "MED14", "MED15", "MED16", "MED17", "MED18", "MED19", "MED20", "MED21", 
            "MED22", "MED23", "MED24", "MED25", "MED26", "MED27", "MED29", "MED30", "MED31", 
            "CCNC", "CDK8", "CDC2L6", "CDK11"]

rel_ppis = ppis[ppis.PreySymbol.isin(mediator)]
rel_ppis = rel_ppis.append(ppis[ppis.BaitSymbol.isin(mediator)])
rel_ppis = rel_ppis.drop_duplicates(subset = ["BaitSymbol", "PreySymbol"])
output = 'processed/merged_ppi_lists/subsets/mediator_all.csv'
rel_ppis.to_csv(output, index = False, na_rep = "NaN")
"""
"""
ppis = 'processed/merged_ppi_lists/ppi-rna/RNA-Roth1-Past-Human-PPI.csv'
ppis = pd.read_csv(ppis)
mediator = ["MED1", "MED4", "MED6", "MED7", "MED8", "MED9", "MED10", "MED11", "MED12", "MED13", 
            "MED13L", "MED14", "MED15", "MED16", "MED17", "MED18", "MED19", "MED20", "MED21", 
            "MED22", "MED23", "MED24", "MED25", "MED26", "MED27", "MED29", "MED30", "MED31", 
            "CCNC", "CDK8", "CDC2L6", "CDK11"]

rel_ppis1 = ppis[ppis.PreySymbol.isin(mediator)]
rel_ppis2 = ppis[ppis.BaitSymbol.isin(mediator)]
rel_ppis2 = rel_ppis2[rel_ppis2.PreySymbol.isin(mediator)]
rel_ppis = rel_ppis1.append(rel_ppis2)
rel_ppis = rel_ppis.drop_duplicates(subset = ["BaitSymbol", "PreySymbol"])

output = 'processed/merged_ppi_lists/subsets/mediator_only.csv'
rel_ppis.to_csv(output, index = False, na_rep = "NaN")
"""

# Generating the cAMP Binding-related subset of proteins
# Uncomment to run.
"""
protein_subset('Human_PPIs/Sig_Sets/cAMP_binding.tab', 
               'processed/merged_ppi_lists/ppi-rna/RNA-Roth1-Past-Human-PPI.csv',
               True,
               'processed/merged_ppi_lists/subsets/cAMP_all.csv')

strict_protein_subset('Human_PPIs/Sig_Sets/cAMP_binding.tab', 
               'processed/merged_ppi_lists/ppi-rna/RNA-Roth1-Past-Human-PPI.csv',
               True,
               'processed/merged_ppi_lists/subsets/cAMP_only.csv')
"""
'''
# END OUTDATED CODE.

# Adding an index row to data to note whether the PPI is in GO subsets
"""
def series_concat(items, gene_list):
    for item in str(items).split("|"):
        gene_list.append(item)
    return items


def index_to_subset(subset_index):
    if subset_index == 0:
        return "Ubiquitin"
    elif subset_index == 1:
        return "Mediator_Complex"
    elif subset_index == 2:
        return "Inflammation"
    else:
        return "Inflammation_Full"

for dataset in ['RNA-Roth1-PPI', 'RNA-Roth1-Past-PPI', 
                'RNA-Roth1-Past-Human-PPI', 'RNA-Roth1-Human-PPI']:
    ppis = pd.read_csv('processed/merged_ppi_lists/ppi-rna/' + dataset + '.csv')

    ubiq = ['polyubiquitin_binding', 'protein_ubiquitation', 'ubiquitin_enzyme_activity',
            'ubiquitin_enzyme_binding', 'ubiquitin_ligase_binding', 'ubiquitin_transferase']
    med_comp = ['mediator_complex']
    inflam = ['response_to_lipopolysaccharide']
    full_inflam = ['response_to_lipopolysaccharide', 'inflammatory_response']


    for subset_index, subset_group in enumerate([ubiq, med_comp, inflam, full_inflam]):
        subset_genes = []
        for subset in subset_group:
            subset_file = 'Human_PPIs/Sig_Sets/' + subset + '.tab'
            subset_data = pd.read_table(subset_file, skiprows = 1)
            subset_genes.extend(subset_data.Gene_Name)
            subset_data.All_Names = subset_data.All_Names.apply(series_concat, 
                                                                args = (subset_genes,))
        subset_genes = pd.Series(subset_genes)
        subset_genes = subset_genes.drop_duplicates()
        subset_genes = subset_genes.drop(subset_genes[subset_genes == 'nan'].index)
        ppis.insert(len(ppis.columns), index_to_subset(subset_index), 
                    ppis.PreySymbol.isin(subset_genes).astype(int))
    
    output_file = 'processed/merged_ppi_lists/subsets/network_annot/' + dataset + '.csv'
    ppis.to_csv(output_file, index = False, na_rep = "NaN")
"""

# Adding an importance column to the PPIs
# The absolute importance is irrelevant, it's a binary encoding of metrics.
# First Digit: 2 stdev from mean in RNA regulation.
# Second Digit: 2 stdev from mean in Translation.
# Third Digit: 2 stdev from mean in protein interactions.
"""
for dataset in ['RNA-Roth1-PPI', 'RNA-Roth1-Past-PPI', 
                'RNA-Roth1-Past-Human-PPI', 'RNA-Roth1-Human-PPI']:
    input_file = 'processed/merged_ppi_lists/subsets/network_annot/' + dataset + '.csv'
    ppis = pd.read_csv(input_file)
    
    count_data = pd.read_csv('processed/merged_ppi_lists/Past_Counts.csv')
    count_data = count_data[count_data.Count > 2]
    is_imp = (1 * (ppis.RNA_Rank >= 1.64).astype(int) + 
              1 * (ppis.RNA_Rank <= -1.64).astype(int) +
              2 * (ppis.trans_change_rank >= 1.64).astype(int) + 
              2 * (ppis.trans_change_rank <= -1.64).astype(int) + 
              4 * (ppis.prot_change_rank >= 1.64).astype(int) + 
              4 * (ppis.prot_change_rank <= -1.64).astype(int))

    ppis.insert(len(ppis.columns), "Important", is_imp)
    output_file = 'processed/merged_ppi_lists/subsets/network_annot/' + dataset + '.csv'
    ppis.to_csv(output_file, index = False, na_rep = "NaN")
"""