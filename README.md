# IAHSummer2020
Code I wrote for my time at Dr. Frederick Roth's lab in the summer of 2020, the year before starting first year undergraduate.

If you want a one-file example of my work, I would recommend protein_ranking.py. If you want a second file example, I would recommend dataset_merge.py.

Due to confidentiality reasons, I cannot include the datasets used. 

The code is slightly disorganized due to the lab changing dataset formats several times. 
This code was originally written for personal use, so comments are fewer than I'd normally use and little UI exists.
Most code was structured into blocks where you uncomment desired blocks, this was performed so common operations could be run in sequence; functions were not used as the operations would change slightly per input and use case and modifying blocks was better than modifying functions.
As this was before first year undergraduate, code methods were very sloppy - though functional! I would add more functions, classes, and type contracts now.

Notable Terms: 
- PPI means Protein Protein Interaction 
- As Y2H (yeast 2 hybrid) was used to detect PPIs, one protein was termed the Bait and the other was termed the Prey.
- Some information on RNA regulation changes, Translateome changes, Proteome changes, and Phosporylation changes were available. These were generally described under the names "RNA", "trans", "prot", and "phos". l2f means Log_2 Fold change. p indicated p-values.

