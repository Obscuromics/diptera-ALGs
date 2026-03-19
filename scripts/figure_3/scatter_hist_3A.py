#!/usr/bin/env python3

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

palette = ["#169e73ff", "#e59d38ff", "#1573afff", "#f0e354ff","#60b5e1ff","#000000","#808080"]
alt_palette = ["#D3D3D3", "#D3D3D3", "#D3D3D3", "#D3D3D3","#D3D3D3","#000000","#D3D3D3"]

df = pd.read_csv('data/chromosomes_vs_ALGs.tsv', sep='\t')

df=df[df['species']!='Panorpa_germanica']

algs = ['d1','d2','d3','d4','d5','d6']
df['alg_total'] = df[algs].sum(axis=1)
df['alg_max_prop'] = df[algs].max(axis=1)/df['alg_total']
df['alg_submax_prop'] = df[algs].apply(lambda row: row.nlargest(2).values[-1],axis=1)/df['alg_total']

df['alg_ratio'] = df['alg_max_prop']/df['alg_submax_prop']
df.replace([np.inf, -np.inf], np.nan, inplace=True)
df['accept_max'] = df['alg_ratio']>1.5
df['alg_max'] = df[algs].idxmax(axis=1)
df.loc[df['accept_max']==False,'alg_max']='Unlabeled'
# % of chromosomes in the dataset
df['accept_max'].sum()/df.shape[0]
df['d6_prop'] = df['d6']/df['alg_total']

df['small_dot'] = np.where(df['d6_prop']>0.5, 1, 0)
print('Small dot chromosome size:')
print(df[df['small_dot']==1]['chromsome_size_b'].describe())
idx = df.groupby('species')['d6_prop'].transform(max) == df['d6_prop']
alg6_df = df[idx]
df.loc[idx,'alg6'] = 1
df.loc[~idx,'alg6'] = 0
df.to_csv('alg6.tsv', sep='\t')

sns.set(font_scale=2)
sns.set_style("white")

fig, axes = plt.subplot_mosaic(
    [["hist", "hist", "."],
     ["scat", "scat", "cbar"]],
    figsize=[8, 12],
    width_ratios=[1, 1, 0.08],  
    gridspec_kw={'hspace': 0.3},
)

ax0 = axes["hist"]
ax1 = axes["scat"]
cax = axes["cbar"]

ax0.set_ylabel('Count')
ax0.set_xticklabels(['0','0', '100', '200', '300', '400', '500'])
ax0.legend()
size_bins = np.arange(0, df['chromsome_size_b'].max(), 20_000_000)
ax0.hist(df[df['small_dot']==0]['chromsome_size_b'], color="#D3D3D3", bins=size_bins, alpha=0.5, label='ALGs 1-5')
ax0.hist(df[df['small_dot']==1]['chromsome_size_b'], color="#000000", bins=size_bins, alpha=0.5, label='ALG 6')
ax0.legend()

data = df[['busco_odb12_complete_count','chromsome_size_b','d6_prop']]

sns.scatterplot(data=data,
            x='chromsome_size_b',
            y='busco_odb12_complete_count',
            ax = ax1,
            hue = 'd6_prop',
            palette='viridis_r',
            s=50,
)

ax1.set_xlabel('Chromosome size (Mb)')
ax1.set_ylabel('Number of BUSCOs')

ax1.set_xticklabels(['0','0', '100', '200', '300', '400', '500'])

ax1.get_legend().remove()
    
scalarmappaple = cm.ScalarMappable(cmap=cm.viridis_r)
plt.colorbar(scalarmappaple, cax=cax, label='Proportion of BUSCOs on\n chromosome assigned to ALG6')

xmin = -20
xmax = df['chromsome_size_b'].max()

ax0.set_xlim(xmin, xmax+20)
ax1.set_xlim(xmin, xmax+20)

plt.savefig('figures/buscoVsize.png',dpi=200, bbox_inches='tight')
plt.savefig('figures/buscoVsize.svg',dpi=200, bbox_inches='tight')