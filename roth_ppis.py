import numpy as np
import pandas as pd

# This file contains almost all operations to process the Roth data. 
# Each operation is its own block with a short explanation. 
# Uncomment the desired operations to run.


# Process the Roth pre-filtered file into the consistent format
"""
roth_file = 'Roth_Results/Filtered_PPI_List_1.csv'
roth_ppis = pd.read_csv(roth_file, usecols = ["DBi", "DBORF", "Entrez", "EntrezORF", "replicate"])
roth_ppis = roth_ppis[roth_ppis.replicate > 1]
roth_ppis = roth_ppis.drop(columns = "replicate")
roth_ppis.DBi[roth_ppis.DBi != 1] = np.NaN
roth_ppis.DBi[roth_ppis.DBi == 1] = 59272
roth_ppis.DBORF[roth_ppis.DBORF == "ORF9Bwu"] = "ORF9C"
roth_ppis.DBORF[roth_ppis.DBORF == "ORF10wu"] = "ORF10"
roth_ppis.to_csv('processed/merged_ppi_lists/Individual_Datasets/Roth1-PPI-GoFiltered.csv', index = False, na_rep = "NaN")
"""

# Adding the Roth filtered file into the CoV-only list
"""
past_filepath = 'processed/merged_ppi_lists/Past-CoV-PPI-List-Sourced.csv'
past_ppis = pd.read_csv(past_filepath)
roth_filepath = 'processed/merged_ppi_lists/Individual_Datasets/Roth1-PPI-GoFiltered.csv'
roth_ppis = pd.read_csv(roth_filepath, usecols = ["BaitSymbol", "PreySymbol"])
roth_ppis.insert(1, "BaitIsCoV", 1)
roth_ppis.insert(3, "PPISymbol", roth_ppis.BaitSymbol + "_" + roth_ppis.PreySymbol)
roth_ppis.insert(4, "Source", "Roth Lab")

sum_ppis = past_ppis.append(roth_ppis, ignore_index = True)

sum_ppis.to_csv('processed/merged_ppi_lists/Roth1-Past-PPI-Sourced.csv', index = False, na_rep = "NaN")
"""


# Processing the Roth filtered formatted file to add IsRoth and IsCoV columns
"""
input_file = 'processed/merged_ppi_lists/Individual_Datasets/Roth1-PPI-GoFiltered.csv'
CoV_PPIs = pd.read_csv(input_file)
CoV_PPIs.insert(2, "BaitIsCoV", [1] * len(CoV_PPIs.BaitSymbol))
CoV_PPIs.insert(5, "IsRoth", [1] * len(CoV_PPIs.BaitSymbol))
CoV_PPIs.to_csv('processed/merged_ppi_lists/Individual_Datasets/Roth1-PPI-GoFiltered2.csv', index = False, na_rep = "NaN")
"""

# Processing Roth individual data to the merged PPI dataset format
"""
input_file = 'processed/merged_ppi_lists/Individual_Datasets/Roth1-PPI-GoFiltered.csv'
CoV_PPIs = pd.read_csv(input_file, usecols = ["BaitSymbol", "PreySymbol"])
CoV_PPIs.insert(1, "BaitIsCoV", 1)
CoV_PPIs.insert(3, "PPISymbol", CoV_PPIs.BaitSymbol + "_" + CoV_PPIs.PreySymbol)
CoV_PPIs.to_csv('processed/merged_ppi_lists/Roth1-PPI.csv', index = False, na_rep = "NaN")
"""
