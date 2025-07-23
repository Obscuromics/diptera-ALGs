#!/usr/bin/env python3

import pandas as pd 
import sys

json_f = sys.argv[1]
json_df = pd.read_json(json_f)

accession = json_df.reports[0]['assembly_accession']

if 'GCA_' in accession:
        acc_str = 'genbank_accession'
if 'GCF_' in accession:
        acc_str = 'refseq_accession'

chromosomes = []
xs = []
ys = []

for sequence in json_df.reports:
    if (sequence['assigned_molecule_location_type'] == 'Chromosome' or sequence['assigned_molecule_location_type'] == 'Linkage Group') & (sequence['role']=='assembled-molecule'):
        chromosomes.append(sequence[acc_str])
        if sequence['chr_name'] == 'X':
            xs.append(sequence[acc_str])
        if sequence['chr_name'] == 'XL':
            xs.append(sequence[acc_str])
        if sequence['chr_name'] == 'XR':
            xs.append(sequence[acc_str])
        if sequence['chr_name'] == 'X1':
            xs.append(sequence[acc_str])
        if sequence['chr_name'] == 'X2':
            xs.append(sequence[acc_str])
        if sequence['chr_name'] == 'Y':
            ys.append(sequence[acc_str])

if xs:
    with open(f'{accession}.x.txt', 'w') as f:
        f.write('\n'.join([str(x) for x in xs]) + '\n')

if ys:
    with open(f'{accession}.y.txt', 'w') as f:
        f.write('\n'.join([str(y) for y in ys]) + '\n')

with open(f'{accession}.chromosomes.txt', 'w') as f:
    f.write('\n'.join([str(chromosome) for chromosome in chromosomes]) + '\n')

