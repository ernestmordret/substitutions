# substitutions

This collection of scripts allows the user to detect and quantify amino acid substitutions from mass spectrometry data. It relies on MaxQuant to perform the initial analysis, and applies a series of filters to separate amino acid substitutions from other PTMs detected by MaxQuant's dependent peptides algorithm.

The scripts are written in Python 2, and require the following modules:  
numpy  
scipy  
biopython  
pandas  
re  
seaborn