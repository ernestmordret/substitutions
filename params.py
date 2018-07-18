#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 10:14:41 2018

@author: ernestmordret
"""

output_dir = 'output'

"""
PARAMETERS FOR DETECT
"""

#specify a path to the fasta file (DNA sequence of protein coding sequences)
path_to_fasta = '/Users/ernestmordret/Library/Mobile Documents/com~apple~CloudDocs/Documents/MS_files/sequences/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa'

#path to MaxQuant's table, allPeptides.txt
path_to_allPeptides = 'allPeptides.txt'

#m/z tolerance, used to filter DP–BP couples that resemble substitutions
#and exclude couples that resemble known PTM
tol = 0.005

#these parameters control the filtering of n-term and c-term modifications
n_term_prob_cutoff = 0.05
c_term_prob_cutoff = 0.05

#defines the minimal threshold on MaxQuant's localisation probability. Only
#DP–BP peptides for which the localisation of the modification has be been
#determined with higher confidence will be retained. 
positional_probability_cutoff = 0.95

#FDR for the final target-decoy FDR procedure
fdr = 0.01

#OPTIONAL PARAMETER : provide a regex to extract the gene name from the
#description of the fasta records.
#ex: regex = 'gene_symbol\:(.*)$'
#if the regex parameter is not defined, the program will use the fasta
#id in place of the gene name.
regex = 'gene_symbol\:(.*)(\s|$)'

"""
PARAMETERS FOR QUANTIFY
"""

excluded_samples = []

#path to MaxQuant's table, allPeptides.txt
path_to_evidence = 'evidence.txt'

#path to MaxQuant's table, allPeptides.txt
path_to_matched_features = 'matchedFeatures.txt'

#path to MaxQuant's table, allPeptides.txt
path_to_peptides = 'peptides.txt'

#m/z tolerance for the fetching unidentified features
mz_tol = 10 * 10**-6 # 10 ppm

#retention time tolerance for the fetching unidentified features
rt_tol = 0.3 # min