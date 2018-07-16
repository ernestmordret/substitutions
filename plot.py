# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 14:45:46 2016

@author: ernestmordret
"""
import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
import seaborn as sns
from params import output_dir

def hamming(s1,s2): return sum(a!=b for a,b in zip(s1,s2))    

def get_codon_table():
	return dict(zip(codons, amino_acids))
	
def get_inverted_codon_table():
	ct = get_codon_table()
	inv_codon_table = {}
	for k, v in ct.iteritems():
		inv_codon_table[v] = inv_codon_table.get(v, [])
		inv_codon_table[v].append(k)
	return inv_codon_table

def prepare_count_matrix(df):
    matrix = pd.DataFrame(data = 0, index = codons, columns=list('ACDEFGHKLMNPQRSTVWY'),dtype=float)
    df = df[df['DP Decoy']!='+']
    df = df[pd.notnull(df['codon'])]
    df['codon'] = df['codon'].map(lambda x: x.replace('T','U'))
    for label in matrix.index:
        if codon_table[label] == '*':
            matrix.loc[label] = float('NaN')
        for col in matrix.columns:
            if (label in inverted_codon_table[col]) or (codon_table[label] +' to '+col in exact_PTM_spec_list):
                matrix.loc[label, col] = float('NaN')
    subs_agg = pd.DataFrame(np.array(zip(*df.groupby(['protein','position','origin','destination','codon']).groups.keys())).T, columns=['protein','position','origin','destination','codon'])												
    for x, l in subs_agg.groupby(['codon', 'destination']).groups.items():
        codon, destination = x
        if (codon in matrix.index) and pd.notnull(matrix.loc[codon,destination]):
            matrix.loc[codon,destination] = len(l)
    matrix.rename(columns={"L": "I/L"},inplace=True)
    return matrix
	
def probe_mismatch(codon1, codon2, pos, spec):
    origin, destination = spec
    for i in range(3):
        if i == pos:
            if codon1[i] != origin or codon2[i] != destination:
                return False
        else:
            if codon1[i] != codon2[i]:
                return False
    return True

bases = 'UCAG'
codons = [a+b+c for a in bases for b in bases for c in bases]

amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
RC = {'A':'U', 'C':'G', 'G':'C', 'U':'A'}

codon_table = get_codon_table()
inverted_codon_table = get_inverted_codon_table()
inverted_codon_table['L'] = inverted_codon_table['L'] + inverted_codon_table['I']
tol = 0.005
MW_dict = {"G": 57.02147, 
            "A" : 71.03712, 
            "S" : 87.03203, 
            "P" : 97.05277, 
            "V" : 99.06842, 
            "T" : 101.04768, 
            "I" : 113.08407, 
            "L" : 113.08407, 
            "N" : 114.04293, 
            "D" : 115.02695, 
            "Q" : 128.05858, 
            "K" : 128.09497, 
            "E" : 129.0426, 
            "M" : 131.04049, 
            "H" : 137.05891,
            "F" : 147.06842, 
            "R" : 156.10112, 
            "C" : 160.030654, #CamCys
            "Y" : 163.0633,
            "W" : 186.07932,
            }

aas_sorted_by_mass = [i[0] for i in sorted(MW_dict.items(),key=lambda x:x[1])]
danger_mods = pd.read_pickle('danger_mods')
exact_PTM_spec = pd.DataFrame(index = aas_sorted_by_mass,
                              columns = aas_sorted_by_mass,
                              dtype = int)

for aa1 in MW_dict.keys():
	for aa2 in MW_dict.keys():
		delta_m = MW_dict[aa2] - MW_dict[aa1]
		exact_PTM_spec.loc[aa1,aa2]=len(danger_mods[(danger_mods['delta_m']<delta_m + 0.0005) & (danger_mods['delta_m']>delta_m - 0.0005) & (danger_mods['site']==aa1)]) > 0

exact_PTM_spec_list = [str(i) + ' to ' + str(j) for i in aas_sorted_by_mass for j in  aas_sorted_by_mass if exact_PTM_spec.loc[i,j]] 

mask = pd.DataFrame(data = False,
                    index = codons,
                    columns = list('ACDEFGHKLMNPQRSTVWY'),
                    dtype = float)	

for label in codons:
	near_cognates = np.array([hamming(i,label)==1 for i in codons])
	reachable_aa = set(np.array(list(amino_acids))[near_cognates])
	mask.loc[label] =[i in reachable_aa for i in 'ACDEFGHKLMNPQRSTVWY']
	
for label in mask.index:
	if codon_table[label] == '*':
		mask.loc[label]=float('NaN')
	for col in mask.columns:
		if (label in inverted_codon_table[col]) or (codon_table[label] +' to '+col in exact_PTM_spec_list):
			mask.loc[label, col] = float('NaN')

subs = pd.read_pickle(os.path.join(output_dir,'subs'))
subs = subs[~subs.decoy]

data = prepare_count_matrix(subs)
log_data = np.log2(data+1)
m = log_data.max().max()

sns.set_style('darkgrid')

font = {
        'size': 16,
        }
sns.set_style({'font.family':'sans-serif', 'font.serif':['Arial']})
fig, ax = plt.subplots(figsize= (16,16))
ax = sns.heatmap(log_data,square = True, cmap = 'Reds', ax = ax, cbar_kws = {"shrink": .3},linewidths=3)
ax.set_ylabel('original codon', fontdict=font)
ax.set_xlabel('destination amino acid', fontdict=font)
ax.tick_params(labelsize=12)
fig.savefig('heatmap.pdf',bbox_inches='tight')

data.to_pickle('unique_substitutions_count_matrix')

