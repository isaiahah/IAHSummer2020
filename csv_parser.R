# Title     : csv_extractor
# Objective : Given a csv file, extract relevant rows/columns
# Created by: isaia
# Created on: 7/1/2020

comment1 <- "Instructions for Use:
file: Path to Input File
to_extract: Columns to extract, as vector of strings
skip_rows: Number of rows at the top to skip
write_to: Path Output File
"
file <- "AP-MS/MMann/PPI-M2si.csv"
to_extract <- c("bait_full_name", "gene_name", "majority_protein_acs", "median_log2", "p_value")
skip_rows <- 1
write_to <- "AP-MS/PPI-Output/MMann-PPI-M2si.csv"

# Functional Code
full_table <- read.csv(file, skip = skip_rows, stringsAsFactors = FALSE,
                       fileEncoding = "UTF-8-BOM", header = TRUE)
cmd_base <- "output_table <- data.frame("
column_vectors <- paste("full_table", to_extract, sep = "$")
columns <- paste(to_extract, column_vectors, sep = " = ", collapse = ", ")
cmd_full <- paste(cmd_base, columns, ", stringsAsFactors = FALSE)")
eval(parse(text = cmd_full))
write.table(output_table, file = write_to, row.names = FALSE, quote = FALSE, sep = ",")