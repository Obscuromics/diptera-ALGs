import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

palette = {'d1':"#169e73ff", 'd2':"#e59d38ff", 'd3':"#1573afff",
         "d4":"#f0e354ff", "d5":"#60b5e1ff", "d6":"black"}

palette_list = ["#169e73ff", "#e59d38ff", "#1573afff", "#f0e354ff", "#60b5e1ff", "black"]

busco_files = [
    ('Panorpa_germanica','OY783203.1','Panorpa_germanica.syngraph.buscos.tsv'), 
    ('Ctenocephalides_felis','NW_020538040.1','Ctenocephalides_felis.syngraph.buscos.tsv')
]
alg_file = 'ALGs_syngraph_diptera.tsv'
alg_df = pd.read_csv(alg_file, sep=' ', header=None)

busco_dfs = {}
for species, alg6_chrom, busco_file in busco_files:
    busco_df = pd.read_csv(busco_file, sep='\t',header=None)
    busco_df = busco_df.merge(alg_df, how='outer', on=0)
    busco_df.rename({0:'gene','1_x':'chromosome',2:'start',3:'stop', '1_y':'alg'}, axis=1, inplace=True)
    busco_df.dropna(inplace=True)
    busco_df = busco_df.groupby(['chromosome', 'alg']).size().reset_index().rename({0:'count'}, axis=1).sort_values(by=['chromosome','count'],ascending=False)
    busco_df['colour'] = busco_df['alg'].map(palette)
    busco_df['alg6_chrom'] = np.where(busco_df['chromosome']==alg6_chrom, 1, 0)
    
    busco_dfs[species] = busco_df

total = busco_dfs['Panorpa_germanica'].groupby('alg')['count'].sum().tolist()
sums = total+total
panorpa_data = busco_dfs['Panorpa_germanica'].groupby(['alg6_chrom', 'alg'])['count'].sum()/sums
panorpa_data = panorpa_data.reset_index()

total = busco_dfs['Ctenocephalides_felis'].groupby('alg')['count'].sum().tolist()
sums = total+total
flea_data = busco_dfs['Ctenocephalides_felis'].groupby(['alg6_chrom', 'alg'])['count'].sum()/sums
flea_data = flea_data.reset_index()

fig, axs = plt.subplots(ncols = 2, nrows=2, figsize=(16, 10))
sns.set(font_scale=2.5)
sns.set_style('white')

pa = axs[0][0]
fa = axs[0][1]
pr = axs[1][0]
fr = axs[1][1]


## UPPER LEFT
pa.set_yticklabels([])
pa.set_yticks([])
pa_ax1 = pa.twinx()
pa_ax1.barh(np.arange(6), panorpa_data[panorpa_data['alg6_chrom']==1]['count'],
                     align='center',
                     height=0.5,
                     color = palette_list
           )
pa_ax1.set_yticklabels(panorpa_data[panorpa_data['alg6_chrom']==1]['alg'])
pa_ax1.invert_xaxis()
pa.set_xticklabels([])
pa_ax1.set_xlim([1,0])

pa_ax1.set_yticks([])

pa.set_title('X chromosome (OY783203)', loc='right')

## LOWER LEFT
pr.set_yticklabels([])
pr.set_yticks([])
pr_ax1 = pr.twinx()
pr_ax1.barh(np.arange(6), panorpa_data[panorpa_data['alg6_chrom']==0]['count'],
                     align='center',
                     height=0.5,
                     color = palette_list
       )
pr_ax1.set_yticklabels(panorpa_data[panorpa_data['alg6_chrom']==1]['alg'])
pr_ax1.invert_xaxis()
pr_ax1.set_xlim([1,0])

pr_ax1.set_yticks([])

pr.set_title('autosomes (n = 20)', loc='right')


## UPPER RIGHT
fa.barh(np.arange(6), flea_data[flea_data['alg6_chrom']==1]['count'],
                     align='center',
                     height=0.5,
                     color = palette_list
       )
fa.set_yticklabels([])
fa.set_xlim([0,1])
fa.set_xticklabels([])
fa.set_title('X chromosome (NW_020538040.1)', loc='left')


## LOWER RIGHT
fr.barh(np.arange(6), flea_data[flea_data['alg6_chrom']==0]['count'],
                     align='center',
                     height=0.5,
                     color = palette_list
       )
fr.set_yticklabels([])
fr.set_xlim([0,1])

fr.set_title(r'autosomes (n $\approx$ 8)', loc='left')
fig.text(0.5, 0, '% of ALG', ha='center')
plt.tight_layout()


plt.savefig('outgroups.png')
plt.savefig('outgroups.svg')