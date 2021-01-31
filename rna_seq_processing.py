import pandas as pd

# This file contains a set of operations to process RNA data.
# Each operation is its own block with a short explanation. 
# Uncomment the desired operations to run.

# Merges the RNA Spreadsheets into one spreadsheet. Adds Averages of l2F and p
"""
A549 = pd.read_csv('RNASeq/rna-seq-A549.csv')
NHBE = pd.read_csv('RNASeq/rna-seq-NHBE.csv')

A549 = A549[['GeneName', 'log2FoldChange', 'padj']]
NHBE = NHBE[['GeneName', 'log2FoldChange', 'padj']]
A549 = A549.rename(columns = {'GeneName' : "GeneName", 'log2FoldChange' : "log2changeA549", 'padj' : "padjA549"})
NHBE = NHBE.rename(columns = {'GeneName' : "GeneName", 'log2FoldChange' : "log2changeNHBE", 'padj' : "padjNHBE"})

rna_seq_combined = pd.DataFrame({'GeneName' : A549.GeneName.append(NHBE.GeneName).drop_duplicates()})
rna_seq_combined = rna_seq_combined.join(A549.set_index('GeneName'), on = 'GeneName')
rna_seq_combined = rna_seq_combined.join(NHBE.set_index('GeneName'), on = 'GeneName')

log2FoldChangeAvg = (rna_seq_combined.log2changeA549 + rna_seq_combined.log2changeNHBE) / 2
padjAvg = (rna_seq_combined.padjA549 + rna_seq_combined.padjNHBE) / 2
rna_seq_combined.insert(5, "log2ChangeAvg", log2FoldChangeAvg)
rna_seq_combined.insert(6, "padjAvg", padjAvg)

rna_seq_combined.to_csv('RNASeq/rna-seq-combined.csv', index = False, na_rep = 'NaN')
"""

# Adds Gene Symbols to the Ferret and Trachea BM RNA data
"""
ferret = pd.read_csv('RNASeq/blanco-mello/bm-mmc3-ferret.csv')
ferret_genes = pd.read_table('RNASeq/blanco-mello/bm-mmc3-ferret-ids.tab', 
                             usecols = ["GeneName", "GeneSymbol"])
ferret = ferret.join(ferret_genes.set_index("GeneName"), on = "GeneName")
ferret = ferret.drop(columns = ["GeneName"])
ferret = ferret[ferret.GeneSymbol.notna()]
ferret = ferret[["GeneSymbol", "SARS-CoV-2 Day1_L2FC", "SARS-CoV-2 Day3_L2FC", 
                 "SARS-CoV-2 Day7_L2FC", "SARS-CoV-2 Day14_L2FC", "IAV Day7_L2FC", 
                 "padj_SARS-CoV-2 Day1", "padj_SARS-CoV-2 Day3", "padj_SARS-CoV-2 Day7", 
                 "padj_SARS-CoV-2 Day14", "padj_IAV Day7"]]
ferret.to_csv('RNASeq/blanco-mello/bm-mmc3-ferret.csv', index = False, na_rep = "NaN")

trachea = pd.read_csv('RNASeq/blanco-mello/bm-mmc3-trachea.csv')
trachea_genes = pd.read_table('RNASeq/blanco-mello/bm-mmc3-trachea-ids.tab', 
                             usecols = ["GeneName", "GeneSymbol"])
trachea = trachea.join(ferret_genes.set_index("GeneName"), on = "GeneName")
trachea = trachea.drop(columns = ["GeneName"])
trachea = trachea[trachea.GeneSymbol.notna()]
trachea = trachea[["GeneSymbol", "SARS-CoV-2 Day3_L2FC", "IAV Day3_L2FC", 
                   "padj_SARS-CoV-2 Day3", "padj_IAV Day3"]]
trachea.to_csv('RNASeq/blanco-mello/bm-mmc3-trachea.csv', index = False, na_rep = "NaN")
"""
