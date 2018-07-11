#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 13:56:13 2017

@author: ernestmordret
"""

import pandas as pd
from scipy.interpolate import interp1d
import re
import numpy as np

def create_modified_seq(modified_seq, destination):
    possible_sites = re.findall('\(([^\)]+)\)', modified_seq)
    best_site = np.argmax([float(i) for i in possible_sites])
    modified_pos_prime = [m.start()-1 for m in re.finditer('\(',modified_seq) ][best_site]
    modified_pos = len(re.sub('\(([^\)]+)\)', '', modified_seq[:modified_pos_prime]))
    base_seq = re.sub('\(([^\)]+)\)', '', modified_seq)
    return base_seq[:modified_pos]+destination+base_seq[modified_pos+1:]

columns_evidence = [u'Sequence', u'Proteins', u'Raw file', u'Fraction',
       u'Experiment', u'MS/MS m/z', u'Charge', u'm/z', u'Mass',
       u'Retention time', u'Retention length', u'Calibrated retention time',
       u'Calibrated retention time start', u'Calibrated retention time finish',
       u'Retention time calibration', u'Match time difference',
       u'Match m/z difference', u'Match q-value', u'Match score',
       u'PEP',
       u'Score', u'Delta score', u'Intensity']

evidence = pd.read_csv('evidence.txt',
                       sep = '\t', usecols = columns_evidence)

calibrate = {}
for i,j in evidence.groupby('Raw file'):
    calibrate[i] = interp1d(j['Retention time'], j['Calibrated retention time'],fill_value="extrapolate")

subs = pd.read_pickle('subs')
subs['modified_sequence'] = subs.apply(lambda row : create_modified_seq(row['DP Probabilities'], row['destination']), axis=1) 
subs['base RT'] =  subs['Retention time'] - subs['DP Time Difference']

for i,j in subs.groupby('Raw file'):
    subs.loc[j.index, 'Calibrated retention time'] = subs.loc[j.index, 'Retention time'].map(lambda x: calibrate[i](x))
    subs.loc[j.index, 'Calibrated base RT'] = subs.loc[j.index, 'base RT'].map(lambda x: calibrate[i](x))
subs['Calibrated DPTD'] = subs['Calibrated retention time'] - subs['Calibrated base RT']
subs['Base m/z'] = subs['m/z'] - (subs['DPMD']/subs['Charge'])

mf = pd.read_csv('matchedFeatures.txt',
                 sep = '\t')
mf.replace(0, np.nan, inplace = True)

pep_head = pd.read_csv('peptides.txt',
                  sep='\t', 
                  nrows=10,
                  index_col = 'Sequence')

intensities = [unicode(i) for i in pep_head.columns if 'Intensity ' in i]

pep_columns = [u'Sequence',u'Intensity'] + intensities

pep = pd.read_csv('peptides.txt',
                  sep = '\t', 
                  usecols = pep_columns,
                  index_col = 'Sequence')
pep.replace(0, np.nan, inplace = True)

mz_tol = 1 * 10**-6 # 10 ppm
rt_tol = 0.3 # min
L = []
RT = []
MZ = []
ratios = []
D = pd.DataFrame(columns=['DP',
                          'Base',
                          'Solubility',
                          'Experiment',
                          'substitution'])

gb = subs.groupby(['Solubility','Charge','DP Base Sequence','modified_sequence'])
substitution_index = 0
for i,j in gb:
    solubility, charge, base_sequence, modified_sequence = i
    ref_dp = j.loc[j['DP PEP'].argmin()]
    charge, dp_mz, base_sequence, dp_rt = ref_dp[['Charge',
                                        'm/z',
                                        'DP Base Sequence',
                                        'Calibrated retention time',
                                        ]]
    
    f1 = mf['Calibrated retention time'] < dp_rt + rt_tol
    f2 = mf['Calibrated retention time'] > dp_rt - rt_tol
    f3 = mf['m/z'] < dp_mz * (1 + mz_tol)
    f4 = mf['m/z'] > dp_mz * (1 - mz_tol)
    f5 = mf['Charge'] == charge
    
    m = mf[f1 & f2 & f3 & f4 & f5]
    l = len(m)

    if l>=1:
        substitution_index += 1
        d = pd.concat([m[intensities].sum(),pep.loc[base_sequence,intensities]], axis=1)
        d.columns=['DP', 'Base']
        d['substitution_index'] = substitution_index
        d['Solubility'] = solubility
        d['Experiment'] = d.index.map(lambda x: x[:-1])
        d['codon'] = ref_dp['codon'] 
        d['destination'] = ref_dp['destination'] 
        d['origin'] = ref_dp['origin'] 
        d['mispairing'] = ref_dp['mispairing'] 
        d['position'] = ref_dp['position']
        d['protein'] = ref_dp['protein']
        d['base_sequence'] = base_sequence
        d['modified_sequence'] = modified_sequence
        D = D.append(d)

D['Ratio'] = D['DP']/D['Base']
D['Sample'] = D.index.map(lambda x: x[-4:])
D['substitution'] = float('NaN')
D.loc[pd.notnull(D['codon']),'substitution'] = D[pd.notnull(D['codon'])].apply(lambda x: x['codon']+' to '+x['destination'],axis=1)
D.reset_index(inplace = True)
D.to_pickle('qSubs')