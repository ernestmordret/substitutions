# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 13:24:37 2016

@author: ernestmordret

"""

import pandas as pd
import numpy as np
import re
from Bio import SeqIO
from itertools import groupby
from operator import itemgetter

#specify the path to a nucleotide fasta reference transcriptome
path_to_fasta = '/Users/ernestmordret/Library/Mobile Documents/com~apple~CloudDocs/Documents/MS_files/sequences/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa'
path_to_allPeptides = 'allPeptides.txt'

#%%
def codonify(seq):
    """
    input: a nucleotide sequence (not necessarily a string)
    output: a list of codons
    """
    seq = str(seq)
    l = len(seq)
    return [seq[i:i+3] for i in xrange(0,l,3)]


def find_proteins(base_seq):
    """
    input: a peptide sequence (string)
    output: the names of proteins containing that sequence 
    """
    tbr = " ".join([names_list[i] for i in np.searchsorted(boundaries_aa-1, SA_search(base_seq, W_aa, sa))])
    if tbr.strip(" ") == '':
        return ''
    else:
        return tbr
                
def fetch_codon(base_seq, modified_pos):
    """
    input: the original aa sequence of a peptide (base_seq),
            and the relative position of the modification.
    output: returns the list of all codons possibly associated 
            with the substitution, presented as a string separated
            by white spaces.
    """
    possible_codons = []
    proteins = find_proteins(base_seq)
    if proteins:
        proteins = proteins.split(" ")
    else:
        return '_'
    for p in proteins:
        if p in record_dict:
            s = record_dict[p].seq
            seq_i = s.translate().find(base_seq)
            i = seq_i + modified_pos
            possible_codons.append(codonify(s)[i])
        else:
            possible_codons.append('_')
    return " ".join(possible_codons)

def fetch_best_codons(modified_seq):
    """
    input: a modified sequence, e.g. LQV(0.91)A(0.09)EK
    output: the list of codons associated with the most likely
            position
    """
    possible_sites = re.findall('\(([^\)]+)\)', modified_seq)
    best_site = np.argmax([float(i) for i in possible_sites])
    modified_pos_prime = [m.start()-1 for m in re.finditer('\(',modified_seq) ][best_site]
    modified_pos = len(re.sub('\(([^\)]+)\)', '',modified_seq[:modified_pos_prime]))
    base_seq = re.sub('\(([^\)]+)\)', '', modified_seq)
    return fetch_codon(base_seq, modified_pos)


def find_substitution_position_local(modified_seq, protein):
    """
    returns the position of a substitutions relative to the start
    of the protein sequence
    """
    possible_sites = re.findall('\(([^\)]+)\)', modified_seq)
    best_site = np.argmax([float(i) for i in possible_sites])
    modified_pos_prime = [m.start()-1 for m in re.finditer('\(',modified_seq) ][best_site]
    modified_pos = len(re.sub('\(([^\)]+)\)', '', modified_seq[:modified_pos_prime]))
    base_seq = re.sub('\(([^\)]+)\)', '', modified_seq)
    s = record_dict[protein].seq
    seq_i = s.translate().find(base_seq)
    i = seq_i + modified_pos
    return i

def find_positions_local(modified_seq, proteins):
    """
    returns the position of a substitutions relative to the start
    of the protein sequence, across all the codons
    """
    positions = []
    for prot in proteins.split(" "):
        positions.append(str(find_substitution_position_local(modified_seq, prot)))
    return " ".join(positions)

def is_gene(record):
    if len(record.seq)%3 != 0:
        return False
    if not record.seq[:3] in {'ATG','GTG','TTG','ATT','CTG'}:
        return False
    if record.seq[-3:].translate()!='*':
        return False
    return True


def refine_localization_probabilities(modified_seq, threshold = 0.05):
	"""
	returns the AAs that were possibly modified (with p > threshold).
	Input: modified sequence (a string of AA with p of each to contain modification: APKV(.7)ML(.3)L means that V was modified with p = .7 and L with p = .3)
	Output: A string with all candidate AAs separated by ';' (V;L).
	"""
    modified_sites = [modified_seq[m.start()-1] for m in re.finditer('\(',modified_seq) ]
    weights = [float(i) for i in re.findall('\(([^\)]+)\)',modified_seq)]
    site_probabilities = {}    
    for aa, weight in zip(modified_sites, weights):
        if aa in site_probabilities:
            site_probabilities[aa] += weight
        else:
            site_probabilities[aa] = weight
    return ";".join([k for k,v in site_probabilities.items() if v>threshold])


def c_term_probability(modified_sequence):
	"""
	Return the probability that C term AA was modified.
	"""
    if modified_sequence[-1] == ')':
        return float(modified_sequence[:-1].split('(')[-1])
    else:
        return 0.0


def n_term_probability(modified_sequence):
    if modified_sequence[1] == '(':
        return float(modified_sequence[2:].split(')')[0])
    else:
        return 0.0


def is_prot_nterm(sequence):
	"""
	Checks whether the observed peptide lies at the N-term of the protein.
	"""
    for start in SA_search(sequence, W_aa, sa):
        if W_aa[start-1] == '*':
            return True
        if W_aa[start-2] == '*':
            return True
    return False


def is_prot_cterm(sequence):
    l=len(sequence)
    for start in SA_search(sequence, W_aa, sa):
        end = start+l
        if W_aa[end] == '*':
            return True
    return False


def get_codon_table():
    return dict(zip(codons, amino_acids))


def get_inverted_codon_table():
    ct = get_codon_table()
    inv_codon_table = {}
    for k, v in ct.iteritems():
        inv_codon_table[v] = inv_codon_table.get(v, [])
        inv_codon_table[v].append(k)
    return inv_codon_table


def hamming(s1,s2): return sum(a!=b for a,b in zip(s1,s2))    


def is_mispairing(row):
    codon = row['codon']
    destination = row['destination']
    if pd.notnull(codon) and pd.notnull(destination):
        if (codon in mask.index) and destination:
            return mask.loc[codon,destination]
        else:
            return 0
    else:
        return float('NaN')


def suffix_array(text, _step=16):
    """Analyze all common strings in the text.
    
    Short substrings of the length _step a are first pre-sorted. The are the 
    results repeatedly merged so that the garanteed number of compared
    characters bytes is doubled in every iteration until all substrings are
    sorted exactly.
    
    Arguments:
        text:  The text to be analyzed.
        _step: Is only for optimization and testing. It is the optimal length
               of substrings used for initial pre-sorting. The bigger value is
               faster if there is enough memory. Memory requirements are
               approximately (estimate for 32 bit Python 3.3):
                   len(text) * (29 + (_size + 20 if _size > 2 else 0)) + 1MB
    
    Return value:      (tuple)
      (sa, rsa, lcp)
        sa:  Suffix array                  for i in range(1, size):
               assert text[sa[i-1]:] < text[sa[i]:]
        rsa: Reverse suffix array          for i in range(size):
               assert rsa[sa[i]] == i
        lcp: Longest common prefix         for i in range(1, size):
               assert text[sa[i-1]:sa[i-1]+lcp[i]] == text[sa[i]:sa[i]+lcp[i]]
               if sa[i-1] + lcp[i] < len(text):
                   assert text[sa[i-1] + lcp[i]] < text[sa[i] + lcp[i]]
    >>> suffix_array(text='banana')
    ([5, 3, 1, 0, 4, 2], [3, 2, 5, 1, 4, 0], [0, 1, 3, 0, 0, 2])
    
    Explanation: 'a' < 'ana' < 'anana' < 'banana' < 'na' < 'nana'
    The Longest Common String is 'ana': lcp[2] == 3 == len('ana')
    It is between  tx[sa[1]:] == 'ana' < 'anana' == tx[sa[2]:]
    """
    tx = text
    size = len(tx)
    step = min(max(_step, 1), len(tx))
    sa = list(range(len(tx)))
    sa.sort(key=lambda i: tx[i:i + step])
    grpstart = size * [False] + [True]  # a boolean map for iteration speedup.
    # It helps to skip yet resolved values. The last value True is a sentinel.
    rsa = size * [None]
    stgrp, igrp = '', 0
    for i, pos in enumerate(sa):
        st = tx[pos:pos + step]
        if st != stgrp:
            grpstart[igrp] = (igrp < i - 1)
            stgrp = st
            igrp = i
        rsa[pos] = igrp
        sa[i] = pos
    grpstart[igrp] = (igrp < size - 1 or size == 0)
    while grpstart.index(True) < size:
        # assert step <= size
        nextgr = grpstart.index(True)
        while nextgr < size:
            igrp = nextgr
            nextgr = grpstart.index(True, igrp + 1)
            glist = []
            for ig in range(igrp, nextgr):
                pos = sa[ig]
                if rsa[pos] != igrp:
                    break
                newgr = rsa[pos + step] if pos + step < size else -1
                glist.append((newgr, pos))
            glist.sort()
            for ig, g in groupby(glist, key=itemgetter(0)):
                g = [x[1] for x in g]
                sa[igrp:igrp + len(g)] = g
                grpstart[igrp] = (len(g) > 1)
                for pos in g:
                    rsa[pos] = igrp
                igrp += len(g)
        step *= 2
    del grpstart
    del rsa
    return sa

def SA_search(P, W, sa):
    lp = len(P)
    n = len(sa)
    l = 0; r = n
    while l < r:
        mid = (l+r) / 2
        a = sa[mid]
        if P > W[a : a + lp]:
            l = mid + 1
        else:
            r = mid
    s = l; r = n
    while l < r:
        mid = (l+r) / 2
        a = sa[mid]
        if P < W[a : a + lp]:
            r = mid
        else:
            l = mid + 1
    return [sa[i] for i in range(s, r)]


def find_homologous_peptide(P):
    if len(SA_search('K' + P, W_aa_ambiguous, sa_ambiguous)) > 0:
        return False
    if len(SA_search('R' + P, W_aa_ambiguous, sa_ambiguous)) > 0:
        return False
    if len(SA_search('*' + P, W_aa_ambiguous, sa_ambiguous)) > 0:
        return False
    return True


def create_modified_seq(modified_seq, destination):
    possible_sites = re.findall('\(([^\)]+)\)', modified_seq)
    best_site = np.argmax([float(i) for i in possible_sites])
    modified_pos_prime = [m.start()-1 for m in re.finditer('\(',modified_seq)][best_site]
    modified_pos = len(re.sub('\(([^\)]+)\)', '', modified_seq[ : modified_pos_prime]))
    base_seq = re.sub('\(([^\)]+)\)', '', modified_seq)
    return base_seq[: modified_pos] + destination + base_seq[modified_pos + 1 :]


#%% 

bases = 'TCAG'
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
RC = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

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

subs_dict = { i+' to '+j : MW_dict[j] - MW_dict[i] for i in MW_dict for j in MW_dict if i!=j}
del subs_dict['L to I']
del subs_dict['I to L']

for k,v in subs_dict.items():
    if k[-1]=='I':
        subs_dict[k+'/L'] = v
        del subs_dict[k]
        del subs_dict[k[:-1]+'L']

sorted_subs, sorted_sub_masses = zip(*sorted(subs_dict.items(), key= lambda x: x[1]))

record_list = []
translated_record_list = []
names_list = []
record_dict = {}
boundaries_aa = [0]
W_codons = []
for record in SeqIO.parse(open(path_to_fasta,'rU'),'fasta'):
    record.seq = record.seq.upper()    
    if is_gene(record):
        translation = unicode(record.seq.translate())
        bits = record.description.split(' ')
        for i in bits:
            if 'gene_symbol' in i:
                record.name = i.split(':')[-1]
        record_list.append(record)
        translated_record_list.append(translation)
        names_list.append(record.name)
        record_dict[record.name] = record
        boundaries_aa.append(boundaries_aa[-1]+len(translation))
        W_codons.extend(list(codonify(record.seq)))

boundaries_aa = np.array(boundaries_aa[1:])
W_aa = ''.join(translated_record_list)
sa = suffix_array(W_aa)
W_aa_ambiguous = W_aa.replace('I','L')
sa_ambiguous = suffix_array(W_aa_ambiguous)

AAs = 'ACDEFGHIKLMNPQRSTVWY'
sites = list(AAs)+['nterm','cterm']

df_iter = pd.read_csv(path_to_allPeptides, sep='\t', chunksize = 10000, iterator=True)
dp = pd.concat( chunk[pd.notnull(chunk['DP Mass Difference'])] for chunk in df_iter)
dp.reset_index(drop=True, inplace=True)

dp['DPMD'] = dp['DP Mass Difference']
dp['DPTD'] = dp['DP Time Difference']
dp['DPAA_noterm'] = dp['DP Probabilities'].map(refine_localization_probabilities)
dp['logFC'] = np.log10(dp['DP Ratio mod/base'])
dp['nterm'] = dp['DP Probabilities'].map(n_term_probability)
dp['cterm'] = dp['DP Probabilities'].map(c_term_probability)
dp['prot_nterm'] = dp['DP Base Sequence'].map(is_prot_nterm)
dp['prot_cterm'] = dp['DP Base Sequence'].map(is_prot_cterm)

danger_mods = pd.read_pickle('danger_mods') #defined in randomize_substutions.py
dp['danger'] = False
for mod in danger_mods.iterrows():
    mod = mod[1]
    position = mod['position']
    site = mod['site']
    delta_m = mod['delta_m']
    
    mass_filter = (delta_m - (2*tol) < dp.DPMD) & (dp.DPMD < delta_m + (2*tol))
    
    term_filter = True
    if position == 'Protein N-term':
        term_filter = (dp.nterm > 0.05) & dp.prot_nterm
    elif position == 'Protein C-term':
        term_filter = (dp.cterm > 0.05) & dp.prot_cterm
    elif position == 'Any N-term':
        term_filter = dp.nterm > 0.05
    elif position == 'Any C-term':
        term_filter = dp.cterm > 0.05
    
    site_filter = True
    if site in amino_acids:    
        site_filter = dp.DPAA_noterm.str.contains(site)
    
    dp.loc[site_filter & term_filter & mass_filter, 'danger'] = True


dp['substitution'] = False
for i in sorted(subs_dict.keys()):
    delta_m = subs_dict[i]
    original_aa = i[0]
    dp.loc[(dp.DPMD > delta_m - tol) & (dp.DPMD < delta_m + tol) & (dp['DPAA_noterm'] == original_aa) & (dp['DP Positional Probability']>0.95) & ~dp['danger'], 'substitution'] = i


#%%
"""
Create mask for mispairing
"""            
mask = pd.DataFrame(data = False, index = codons, columns=list('ACDEFGHKLMNPQRSTVWY'),dtype=float)    
for label in codons:
    near_cognates = np.array([hamming(i,label)==1 for i in codons])
    reachable_aa = set(np.array(list(amino_acids))[near_cognates])
    mask.loc[label] =[i in reachable_aa for i in 'ACDEFGHKLMNPQRSTVWY']


for label in mask.index:
    for col in mask.columns:
        if label in inverted_codon_table[col]:
            mask.loc[label, col] = float('NaN')

#%%
subs = dp[dp['substitution']!=False].copy()
subs['proteins'] = subs['DP Base Sequence'].map(find_proteins)
subs['protein'] = subs['proteins'].map(lambda x: x.split(' ')[0] if len(x)>0 else float('NaN'))
subs = subs[pd.notnull(subs['protein'])]

subs['codons'] = float('NaN')
subs.loc[ subs['DP Positional Probability'] > 0.95, 'codons' ] = subs[ subs['DP Positional Probability'] > 0.95 ]['DP Probabilities'].map(fetch_best_codons)
subs['codon'] = subs['codons'].map(lambda x: x.split(' ')[0] if len(set(x.split(' ')))==1 else float('NaN'))
subs['destination'] = subs['substitution'].map(lambda x: x[-1] if x else False)
subs['origin'] = subs['substitution'].map(lambda x: x[0] if x else False)
subs['mispairing'] = subs.apply(is_mispairing,1)    

subs['positions'] = subs.apply(lambda row : 
        find_positions_local(row['DP Probabilities'],row['proteins']),
        axis=1)
subs['position'] = subs['positions'].map(lambda x: int(x.split(' ')[0]) if len(set(x.split(' ')))==1 else float('NaN'))
subs['modified_sequence'] = subs.apply(lambda row : create_modified_seq(row['DP Probabilities'], row['destination']), axis=1)
subs['modified_sequence'] = subs['modified_sequence'].map(lambda x: x.replace('I','L'))
subs = subs[subs['modified_sequence'].map(lambda x: find_homologous_peptide(x))]

subs.sort_values('DP PEP', inplace = True)
subs['decoy'] = pd.notnull(subs['DP Decoy'])
cut_off = np.max(np.where(np.array([i/float(j) for i,j in zip(subs['decoy'].cumsum(), range(1,len(subs)+1))])<0.01))
subs = subs.iloc[:cut_off+1]

subs['Solubility'] = subs['Raw file'].map(lambda x: x[-7])
subs['Fraction'] = subs['Raw file'].map(lambda x: int(x[-1]))
subs['Experiment'] = subs['Raw file'].map(lambda x: x[-6:-2])
#%%
subs = subs[~subs.decoy]

subs.to_pickle('subs')
subs.to_csv('subs.csv')
dp.to_pickle('dp')
