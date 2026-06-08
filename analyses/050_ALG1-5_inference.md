### Running syngraph 

TODO: add running nextflow instructions

#### Instructions on modifying the syngraph run

```bash
module load syngraph/0.0.1-c3
```

The included busco files are in:

```bash
data/results/infer_alg_diptera/syngraph_busco_tables/
```

If you want to change the included files modify them in place, or to exclude them move them to:

```bash
data/results/infer_alg_diptera/excluded_syngraph_busco_tables/
```
and run

```bash
scripts/prep_syngraph_tree.py # TODO: THIS SCRIPT IS MISSING
```

which will create a new tree, modifying the iqtree output by removing the leaves which have files in /excluded_syngraph_busco_tables. Theres also a line you can change from 'None' to the outgroup clade file if you want and then it will reroot the tree.

Then you just have to rerun syngraph

```bash
syngraph build -d data/results/infer_alg_diptera/syngraph_busco_tables -m -o data/results/infer_alg_diptera/syngraph/diptera.pruned.syngraph_build

syngraph infer -g data/results/infer_alg_diptera/syngraph/diptera.pruned.syngraph_build.pickle -t data/results/infer_alg_diptera/diptera.pruned.newick -m 89 -r 2 -a quick -s Phlebotomus_papatasi -o data/results/infer_alg_diptera/syngraph/diptera.pruned.syngraph_infer

syngraph tabulate -g data/results/infer_alg_diptera/syngraph/diptera.pruned2.syngraph_infer.with_ancestors.pickle -o data/results/infer_alg_diptera/syngraph/diptera.pruned.syngraph_tabulate
```
If you run this as is the output files will overwrite the last analysis I did so could edit those.

#### Selection of species

```R
# cd ~/Projects/ALGs/syngraph

library(gsheet)
library(ape)

# tree <- read.tree("./data/results/infer_alg_diptera/iqtree/diptera.supermatrix.phy.treefile")
sams_tree <- read.tree('data/syngraph/diptera.pruned.syngraph_infer.min_dist.m165.further_pruned/diptera.pruned.syngraph_infer.min_dist.m165.newick.txt')

all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?usp=sharing", format='csv'),
                            stringsAsFactors = F, header = T, check.names = F)

# all_genome_data <- all_genome_data[all_genome_data[, 'TO ADD'] %in% c('KEEP', 'OUTGROUP'), ]
# species_in_the_tree <- sapply(all_genome_data[, 'species'], function(x){x %in% tree$tip.label})
species_in_sams_tree <- sapply(all_genome_data[, 'species'], function(x){x %in% sams_tree$tip.label})


species_with_decent_busco_and_allowed_to_be_used <- !(all_genome_data[, 'excluded_from_ALG_inference'] == "removed (B)" | all_genome_data[, 'excluded_from_ALG_inference'] == "removed (U)")
drosohilids <- all_genome_data[, 'family'] == 'Drosophilidae'
mosquitoes <- all_genome_data[, 'family'] %in% c('Culicidae', 'Chironomidae')

any(all_genome_data[, 'family'] == 'Cecidomyiidae')

rearranging <- all_genome_data[, 'family'] %in% c('Chironomidae', 'Glossinidae', 'Cecidomyiidae', 'Anisopodinae')


# table(all_genome_data[mosquitoes & species_in_sams_tree, 'family'])


# species_to_list <- species_in_the_tree[species_in_the_tree[, 'excluded_from_ALG_inference'] == '', 'species']
# species_to_list <- all_genome_data[drosohilids & species_in_the_tree & species_with_decent_busco_and_allowed_to_be_used, 'species']
species_to_list <- all_genome_data[rearranging, 'species']

write.table(data.frame(sp = species_to_list), file = 'data/pruning9_messy_remove_list.txt', col.names = F, quote = F, row.names = F)
```

now move all the species out, and back only those Julia used.

```bash
# mv ./data/results/infer_alg_diptera/syngraph_busco_tables/* ./data/results/infer_alg_diptera/excluded_syngraph_busco_tables/
mv ./data/results/infer_alg_diptera/excluded_syngraph_busco_tables/* ./data/results/infer_alg_diptera/syngraph_busco_tables/

# while read species; do echo $species; mv ./data/results/infer_alg_diptera/excluded_syngraph_busco_tables/"$species".syngraph.buscos.tsv ./data/results/infer_alg_diptera/syngraph_busco_tables/; done < data/pruning8_culicomoprha_keep_list.txt

while read species; do echo $species; mv ./data/results/infer_alg_diptera/syngraph_busco_tables/"$species".syngraph.buscos.tsv ./data/results/infer_alg_diptera/excluded_syngraph_busco_tables/; done < data/pruning9_messy_remove_list.txt
```

#### Kamil's attempt 3

The first attempt was to run syngraph with `-m 100`, that worked for the ALGs, but internal nodes were messed up. I will redo the run while excluding all those Julia excluded.

```bash
screen -r
bsub -G team360 -Is -n 16 -M 10240 -R "select[mem>10240] rusage[mem=10240]" bash -l
cd /data/tol/teams/jaron/lustre/users/se13/diptera_alg
module load syngraph/0.0.1-c3
```

```bash
 while read species; do echo echo $species; mv ./data/results/infer_alg_diptera/syngraph_busco_tables/"$species".syngraph.buscos.tsv ./data/results/infer_alg_diptera/excluded_syngraph_busco_tables/; done < data/pruning2_exclusion_list.txt
```

and run

```bash
./scripts/prep_syngraph_tree.py 
```

```bash
syngraph build -d data/results/infer_alg_diptera/syngraph_busco_tables -m -o data/results/infer_alg_diptera/syngraph/diptera.pruned.syngraph_build

syngraph infer -g data/results/infer_alg_diptera/syngraph/diptera.pruned.syngraph_build.pickle -t data/results/infer_alg_diptera/diptera.pruned.newick -m 100 -r 2 -a quick -s Bibio_marci -o data/results/infer_alg_diptera/syngraph/diptera.pruned3.syngraph_infer

syngraph tabulate -g data/results/infer_alg_diptera/syngraph/diptera.pruned3.syngraph_infer.with_ancestors.pickle -o data/results/infer_alg_diptera/syngraph/diptera.pruned3.syngraph_tabulate
```

#### Attempt 4 - a range of ms

```bash
while read species; do echo echo $species; mv ./data/results/infer_alg_diptera/syngraph_busco_tables/"$species".syngraph.buscos.tsv ./data/results/infer_alg_diptera/excluded_syngraph_busco_tables/; done < data/pruning4_exclusion_list.txt

./scripts/prep_syngraph_tree.py 

syngraph build -d data/results/infer_alg_diptera/syngraph_busco_tables -m -o data/results/infer_alg_diptera/syngraph/diptera.pruned.syngraph_build

for m in 85 90 95 100 105 110 115 120 125 130 135 140 145 150 155 160 165 170 175 180 185 190 195 200; do echo $m; 
syngraph infer -g data/results/infer_alg_diptera/syngraph/diptera.pruned.syngraph_build.pickle -t data/results/infer_alg_diptera/diptera.pruned.newick -m $m -r 2 -a quick -s Bibio_marci -o data/results/infer_alg_diptera/syngraph/mrange/diptera.pruned3.m"$m".syngraph_infer

syngraph tabulate -g data/results/infer_alg_diptera/syngraph/mrange/diptera.pruned3.m"$m".syngraph_infer.with_ancestors.pickle -o data/results/infer_alg_diptera/syngraph/mrange/diptera.pruned3.m"$m".syngraph_tabulate
done
```

These did not work at all. I have no idea what happened. I think I need to come back to basics. Can I regenerate Julia's results if I use the same genomes? Perhpas at the beginning I will remove those that we were asked to remove, but that's nearly identical dataset. None of the removed flies is likely going to change much.

#### Attempt 5 - getting back to species Julia selected

and now rerun syngraph

```bash
./scripts/prep_syngraph_tree.py 

syngraph build -d data/results/infer_alg_diptera/syngraph_busco_tables -m -o data/results/infer_alg_diptera/syngraph/diptera.juliaset.syngraph_build

for m in 87 100 150 200; do echo $m; 
syngraph infer -g data/results/infer_alg_diptera/syngraph/diptera.juliaset.syngraph_build.pickle -t data/results/infer_alg_diptera/diptera.pruned.newick -m $m -r 2 -a quick -s Bibio_marci -o data/results/infer_alg_diptera/syngraph/mrange/diptera.juliaset.m"$m".syngraph_infer > data/results/infer_alg_diptera/syngraph/mrange/diptera.juliaset.m"$m".log

syngraph tabulate -g data/results/infer_alg_diptera/syngraph/mrange/diptera.juliaset.m"$m".syngraph_infer.with_ancestors.pickle -o data/results/infer_alg_diptera/syngraph/mrange/diptera.juliaset.m"$m".syngraph_tabulate
done
```

This also failed for all `n1` nodes predicting 2, 2, 3 and 3 ALGs.

##### Just drosophilids

```bash
mv ./data/results/infer_alg_diptera/syngraph_busco_tables/* data/results/infer_alg_diptera/excluded_syngraph_busco_tables/
while read species; do echo $species; mv ./data/results/infer_alg_diptera/excluded_syngraph_busco_tables/"$species".syngraph.buscos.tsv ./data/results/infer_alg_diptera/syngraph_busco_tables/; done < data/pruning6_droshophilidae_keep_list.txt

./scripts/prep_syngraph_tree.py 
./syngraph build -d data/results/infer_alg_diptera/syngraph_busco_tables -m -o data/results/infer_alg_diptera/syngraph/drosophilidae.syngraph_build

./syngraph infer -g data/results/infer_alg_diptera/syngraph/drosophilidae.syngraph_build.pickle -t data/results/infer_alg_diptera/diptera.pruned.newick -m 100 -r 2 -a quick -s Drosophila_melanogaster -o data/results/infer_alg_diptera/syngraph/drosophilidae.syngraph_infer > data/results/infer_alg_diptera/syngraph/drosophilidae.syngraph.log

# syngraph tabulate -g data/results/infer_alg_diptera/syngraph/diptera.pruned3.syngraph_infer.with_ancestors.pickle -o data/results/infer_alg_diptera/syngraph/diptera.pruned3.syngraph_tabulate

##### 5 drosophilids

mv ./data/results/infer_alg_diptera/syngraph_busco_tables/* ./data/results/infer_alg_diptera/excluded_syngraph_busco_tables/

while read species; do echo $species; mv ./data/results/infer_alg_diptera/excluded_syngraph_busco_tables/"$species".syngraph.buscos.tsv ./data/results/infer_alg_diptera/syngraph_busco_tables/; done < data/pruning7_5_droshophilidae_keep_list.txt

./scripts/prep_syngraph_tree.py 

./syngraph build -d data/results/infer_alg_diptera/syngraph_busco_tables -m -o data/results/infer_alg_diptera/syngraph/drosophilidae_minimal5.syngraph_build

./syngraph infer -g data/results/infer_alg_diptera/syngraph/drosophilidae_minimal5.syngraph_build.pickle -t data/results/infer_alg_diptera/diptera.pruned.newick -m 100 -r 2 -a quick -s Drosophila_miranda -o data/results/infer_alg_diptera/syngraph/drosophilidae_minimal5.syngraph_infer 

```

#### Mosquitos

I am also adding Bibio marci as an outgorup;

```bash
./scripts/prep_syngraph_tree.py 

./syngraph build -d data/results/infer_alg_diptera/syngraph_busco_tables -m -o data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_build

./syngraph infer -g data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_build.pickle -t data/results/infer_alg_diptera/diptera.pruned.newick -m 160 -r 2 -a quick -s Bibio_marci -o data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_infer > data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_infer.log

./syngraph tabulate -g data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_infer.with_ancestors.pickle -o data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_infer.syngraph_tabulate

SYNGNWK=data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_infer.newick.txt
SYNGTAB=data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_infer.syngraph_tabulate.table.tsv

Rscript scripts/plot_syngraph_tree_with_nodes.R -i  data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_infer.newick.txt -o data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_infer.tree.pdf -c data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_infer.clusters.tsv

mkdir -p figures
Rscript ../../diptera-ALGs/scripts/20250718_tree2internal_nodes.R -t $SYNGNWK -o data/results/infer_alg_diptera/syngraph/culicomorpha_nodes

while read node; do 
    echo $node; 
    # node
    awk -v node="$node" 'NR==1 { for (i=1; i<=NF; i++){f[$i] = i} }{ if( $(f[node"_seq"]) != "NA" ){ print $1 "\t" $(f[node"_seq"])}}' $SYNGTAB > data/results/infer_alg_diptera/syngraph/culicomorpha_nodes/"$node"_asgn.tsv
done < data/results/infer_alg_diptera/syngraph/culicomorpha_nodes_internal_nodes.tsv 

# n55 is where the shit hits the fan
ALGs=data/diptera.pruned.syngraph_infer.min_dist_outgroups.further_pruned.m165_n3.tsv

mkdir -p figures/node_plots/
rm figures/node_plots/*
while read node; do 
    echo $node
    Rscript scripts/plot_node.R -n $node -a $ALGs -o figures/node_plots/$node
done < data/results/infer_alg_diptera/syngraph/culicomorpha_nodes_internal_nodes.tsv

```

#### Rearranging removal

```bash
# ./scripts/prep_syngraph_tree.py 

# ./syngraph build -d data/results/infer_alg_diptera/syngraph_busco_tables -m -o data/results/infer_alg_diptera/syngraph/no_rearrang3.syngraph_build

# ./syngraph infer -g data/results/infer_alg_diptera/syngraph/no_rearrang3.syngraph_build.pickle -t data/results/infer_alg_diptera/diptera.pruned.newick -m 165 -r 2 -a quick -s Bibio_marci -o data/results/infer_alg_diptera/syngraph/no_rearrang3.syngraph_infer > data/results/infer_alg_diptera/syngraph/no_rearrang3.syngraph_infer.log

./syngraph tabulate -g data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_infer.with_ancestors.pickle -o data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_infer.syngraph_tabulate
# SAM RUN FOR ME
data/diptera.no_woodgnat.mindist.m165.newick.txt
SYNGNWK=data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_infer.newick.txt
SYNGTAB=data/results/infer_alg_diptera/syngraph/culicomorpha.syngraph_infer.syngraph_tabulate.table.tsv
```

### Parsing results

Getting data from the cluster and preprocessing

```bash
scp -r farm:/lustre/scratch126/tol/teams/jaron/users/se13/diptera_alg/data/results/infer_alg_diptera/syngraph/diptera.pruned3.syngraph_infer.newick.txt data/syngraph
scp -r farm:/lustre/scratch126/tol/teams/jaron/users/se13/diptera_alg/data/results/infer_alg_diptera/syngraph/diptera.pruned3.syngraph_tabulate.table.tsv data/syngraph
# cd data/syngraph
# head -5 diptera.syngraph_tabulate.table.tsv > diptera.syngraph_tabulate.table_sample.tsv # sampling the header of the most important file
```

just to figure out what is in the file

```R
syngraph_tab <- read.table('diptera.pruned.syngraph_tabulate.table.tsv', header = T, sep = '\t')
syngraph_cols <- colnames(syngraph_tab)

syngraph_cols[grepl("n2_", syngraph_cols)]
which(grepl("n222_", syngraph_cols))
# 1238 1239 1240

which(grepl("n1_", syngraph_cols))
# 1316 1317 1318

# which(grepl("n0_", syngraph_cols))
# good nada
```

```bash
awk 'NR==1 { for (i=1; i<=NF; i++){f[$i] = i} }{ print $1 "\t" $(f["n1_seq"]) "\t" $(f["n2_seq"]) }' diptera.pruned.syngraph_tabulate.table.tsv > n1_2_assingments.txt
```

```R
sygraph <- read.table('n1_2_assingments.txt', header = T)

table(paste(sygraph[, 2], sygraph[, 3])) # 4 LGs, not 5
table(sygraph[, 3]) # 4 LGs, not 5

sygraph[, 'ALGs'] <- NA
# ALG 1
sygraph[which(sygraph[, 'n1_seq'] == 'n1_2' & sygraph[, 'n2_seq'] == 'n2_1'), 'ALGs'] <- "d1"

# ALG 2
sygraph[which(sygraph[, 'n1_seq'] == 'n1_1' & sygraph[, 'n2_seq'] == 'n2_3'), 'ALGs'] <- "d2"

# ALG 3
sygraph[which(sygraph[, 'n1_seq'] == 'n1_3' & sygraph[, 'n2_seq'] == 'n2_2'), 'ALGs'] <- "d3"

# ALG 4
sygraph[which(sygraph[, 'n1_seq'] == 'n1_4' & sygraph[, 'n2_seq'] == 'n2_4'), 'ALGs'] <- "d4"

sygraph_assigned <- sygraph[!is.na(sygraph[, 'ALGs']), ]

write.table(sygraph_assigned[, c('marker', 'ALGs')], file = 'syngraph.pruned.ALGs.inferred.tsv', col.names = F, row.names = F, quote = F)
```

`diptera.syngraph_infer.clusters.tsv` reports the number of inferred linkage groups;
Possible to plot on a tree, Arif has done long time ago;


### m = 100

```bash
awk 'NR==1 { for (i=1; i<=NF; i++){f[$i] = i} }{ print $1 "\t" $(f["n1_seq"]) "\t" $(f["n2_seq"]) }' diptera.pruned.100.syngraph_tabulate.table.tsv > pruned.m100.n1_2_assingments.txt
```

```R
sygraph <- read.table('pruned.m100.n1_2_assingments.txt', header = F)
colnames(sygraph) <- c('busco', 'n1_seq', 'n2_seq')

sort(table(paste(sygraph[, 2], sygraph[, 3]))) # 4 LGs, not 5
table(sygraph[, 3]) # 4 LGs, not 5

sygraph[, 'ALGs'] <- NA
# ALG 1
sygraph[which(sygraph[, 'n1_seq'] == 'n1_1' & sygraph[, 'n2_seq'] == 'n2_2'), 'ALGs'] <- "d1"

# ALG 2
sygraph[which(sygraph[, 'n1_seq'] == 'n1_3' & sygraph[, 'n2_seq'] == 'n2_5'), 'ALGs'] <- "d2"

# ALG 3
sygraph[which(sygraph[, 'n1_seq'] == 'n1_2' & sygraph[, 'n2_seq'] == 'n2_1'), 'ALGs'] <- "d3"

# ALG 4
sygraph[which(sygraph[, 'n1_seq'] == 'n1_4' & sygraph[, 'n2_seq'] == 'n2_4'), 'ALGs'] <- "d4"

# ALG 5
sygraph[which(sygraph[, 'n1_seq'] == 'n1_5' & sygraph[, 'n2_seq'] == 'n2_3'), 'ALGs'] <- "d5"

sygraph_assigned <- sygraph[!is.na(sygraph[, 'ALGs']), ]

write.table(sygraph_assigned[, c('busco', 'ALGs')], file = 'syngraph.pruned.100.ALGs.inferred.tsv', col.names = F, row.names = F, quote = F)
```

Looking at the intermediate nodes

```R
require('ape')
# require('ggtree')

tree <- read.tree("data/syngraph/diptera.pruned.syngraph_infer.newick.txt")

pdf('figures/syngraph_tree_with_node_names.pdf', height = 60, width = 30)
    plot(tree,  show.node.label = TRUE, cex=0.8, col= "blue")
dev.off()

# sort(as.numeric(substr(tree$node.label, 2, nchar(tree$node.label))))

write.table(data.frame(node = tree$node.label[tree$node.label != '']), 'data/syngraph/all_internal_nodes.tsv', quote = F, col.names = F, row.names = F)
```

### Excluding Chironomids and Tse tse - Final?

I will try to write it in the way that if we have to redo this, I will be able to regenerate all the things...

```bash
SYNGTAB=data/syngraph/diptera.no_cecid.syngraph_infer.min_dist.m165.syngraph_tabulate.table.tsv
SYNGNWK=data/syngraph/diptera.no_cecid.syngraph_infer.min_dist.m165.newick.txt
SYNGNODES=data/syngraph/diptera.no_cecid.syngraph_infer.min_dist.m165.all
# data/syngraph/diptera.pruned3.syngraph_infer.newick.txt

rm data/syngraph/node_assignments/* # clean previous node assignments
Rscript scripts/20250718_tree2internal_nodes.R -t $SYNGNWK -o $SYNGNODES
# this also plots the tree as "figures/syngraph_tree_with_nodes.pdf"

while read node; do 
    echo $node; 
    # node
    awk -v node="$node" 'NR==1 { for (i=1; i<=NF; i++){f[$i] = i} }{ if( $(f[node"_seq"]) != "NA" ){ print $1 "\t" $(f[node"_seq"])}}' $SYNGTAB > data/syngraph/node_assignments/"$node"_asgn.tsv
done < "$SYNGNODES"_internal_nodes.tsv

# $SYNGNWK
```

and now using `n1` and `n2` nodes, define ALGs1 - 5. The ALG 6 is added too, but that is described in the following .md file.

```bash
Rscript scripts/20250719_ALG_definition.R -o data/diptera.no_chiro.syngraph_infer.min_dist.m165_n3.tsv -n n3
Rscript scripts/20250719_ALG_definition.R -o data/diptera.no_chiro.syngraph_infer.min_dist.m165_n1.tsv -n n1
Rscript scripts/20250719_ALG_definition.R -o data/diptera.no_chiro.syngraph_infer.min_dist.m165_n1_n2.tsv
ALGs=data/diptera.no_ceci.syngraph_infer.min_dist.m165_n1.tsv
Rscript scripts/20250719_ALG_definition.R -o $ALGs -n n1
```

### Painting of internal nodes

```bash
rm figures/node_plots/*
while read node; do 
    echo $node
    Rscript scripts/plot_node.R -n $node -a $ALGs -o figures/node_plots/$node
done < data/syngraph/diptera.no_ceci.some_nodes.tsv
echo n14 n3 n7 n13 > data/syngraph/diptera.no_ceci.some_nodes.tsv
# data/syngraph/diptera.no_chiro.syngraph_infer.min_dist.m165.all_internal_nodes.tsv
```

```bash
for family in Tipulidae Bibionidae Culicidae Ptychopteridae Chironomidae Psychodidae; do Rscript scripts/ALG_painter_ME_building_blocks.R -o figures/wo_chironomids/"$family"_no_chiro_m165_n12 -a data/diptera.no_chiro.syngraph_infer.min_dist.m165_n1_n2.tsv -f $family; done
```

### no woodgnats, no plecia - ESEB final

We had to gradually remove all the branches of the tree where there were too many changes happening - the limit of detectable rearrangements is number of chromosomes. The larges clade removed are chironomidae, with remarkably rearranged genomes. We also had to removed quite a few bibionomorphans as there are not enough genomes. This is an inference that reconstructs  
- the ancestral 5 ALGs (with the ALG6 makes 6)
- moreless the same 5 ALGs as the common ancestor of bibionomorpha and other superfamilies
- the fragmented Brachyceran ALGs
-  

```bash
# ./syngraph tabulate -g ../../diptera-ALGs/data/syngraph/diptera.no_plecia.mindist.m165.with_ancestors.pickle -o ../../diptera-ALGs/data/syngraph/diptera.no_plecia.mindist.m165.with_ancestors_tabulate

RUN=diptera.no_plecia.mindist.m165
SYNGTAB=data/syngraph/"$RUN".with_ancestors_tabulate.table.tsv
SYNGNWK=data/syngraph/"$RUN".newick.txt
SYNGNODES=data/syngraph/"$RUN".all

Rscript scripts/20250718_tree2internal_nodes.R -t $SYNGNWK -o $SYNGNODES

while read node; do 
    echo $node; 
    # node
    awk -v node="$node" 'NR==1 { for (i=1; i<=NF; i++){f[$i] = i} }{ if( $(f[node"_seq"]) != "NA" ){ print $1 "\t" $(f[node"_seq"])}}' $SYNGTAB > data/syngraph/node_assignments/"$node"_asgn.tsv
done < $SYNGNODES$"_internal_nodes.tsv" 
```

Now I will also use the syngraph tree to filter species and guess the set of BUSCOs on the ALG6 as those that are present on at least 40% of the dots.

```bash
Rscript scripts/20250811_ALG6.R -t $SYNGNWK
```

This created `data/ALG6_BUSCOs.tsv`. Now using `n1` and `n2` nodes, define ALGs1 - 5. The ALG 6 will be added to it too.

```bash
Rscript scripts/20250719_ALG_definition.R -o data/"$RUN"_n3.tsv -n n3
Rscript scripts/20250719_ALG_definition.R -o data/"$RUN"_n1.tsv -n n1
Rscript scripts/20250719_ALG_definition.R -o data/"$RUN"_n1_n2.tsv
```

```bash
mkdir -p figures/ALG_testing

for family in Ceratopogonidae; do 
    Rscript scripts/ALG_painter_ME_building_blocks.R -o figures/ALG_testing/"$RUN"_"$family"_n1_n2 -a data/"$RUN"_n1_n2.tsv -f $family;
    Rscript scripts/ALG_painter_ME_building_blocks.R -o figures/ALG_testing/"$RUN"_"$family"_n1 -a data/"$RUN"_n1.tsv -f $family;
    Rscript scripts/ALG_painter_ME_building_blocks.R -o figures/ALG_testing/"$RUN"_"$family"_n3 -a data/"$RUN"_n3.tsv -f $family;
done

rm -r figures/ALG_testing 
```

These figures were showing painting of tips from both sides of the tree using different ways of defining the ALGs. In the end, the least noisy paintings are by the intersection of nodes n1 and n2. I will keep those as the ALG definition.

```bash
ALGs=data/"$RUN"_n1_n2.tsv
```

```bash
mkdir -p figures/node_plots/
rm figures/node_plots/*
while read node; do 
    echo $node
    Rscript scripts/plot_node.R -n $node -a $ALGs -o figures/node_plots/$node
done < $SYNGNODES$"_internal_nodes.tsv"


# mkdir -p figures/wo_chironomids/gnats
# for family in Anisopodidae Bibionidae Cecidomyiidae Sciaridae Mycetophilidae; do Rscript scripts/ALG_painter_ME_building_blocks.R -o figures/wo_chironomids/gnats/"$family"_no_chiro_m165_n12 -a data/diptera.no_chiro.syngraph_infer.min_dist.m165_n1_n2.tsv -f $family; done
```

```bash
cd data/syngraph/alg6
tail -n+2 diptera.dot.mindist.m1.syngraph_infer.reconstruction_order.tsv | cut -f 2 | sort | uniq > all_nodes
mkdir nodes
while read node; do 
    echo $node; 
    # node
    awk -v node="$node" 'NR==1 { for (i=1; i<=NF; i++){f[$i] = i} }{ if( $(f[node"_seq"]) != "NA" ){ print $1 "\t" $(f[node"_seq"])}}' diptera.dot.syngraph_tabulate.table.tsv > nodes/"$node"_asgn.tsv
done < all_nodes
cd ../../..
```

```R
ALG6_Julia <- read.table('data/ALG6_BUSCOs.tsv')

ALG6_syngraph <- read.table('data/syngraph/alg6/nodes/n1_asgn.tsv')

ALG6_syngraph[, 1] %in% ALG6_Julia[, 1]
ALG6_Julia[, 1] %in% ALG6_syngraph[, 1]
# good overlap, but far from perfect

ALGs <- read.table('tables/ALGs_syngraph_diptera.tsv')
row.names(ALGs) <- ALGs[, 1]

already_defined <- ALG6_syngraph[ALG6_syngraph[, 1] %in% ALGs[, 1], 1]
ALGs[already_defined, ][ALGs[already_defined, 2] != 'd6', ] 

#  162912at7147 d2
#  164810at7147 d1
#  60829at7147 d1
#  75965at7147 d1
#  61998at7147 d1
#  173520at7147 d1
#  245796at7147 d1

wo_old_ALG6 <- ALGs[ALGs[, 2] != 'd6', ]

sum(wo_old_ALG6[, 1] %in% ALG6_syngraph[, 1])
wo_old_previous_ALG6 <- wo_old_ALG6[!(wo_old_ALG6[, 1] %in% ALG6_syngraph[, 1]), ]

sum(ALG6_syngraph[, 1] %in% wo_old_ALG6[, 1])
wo_prev_asn_ALG6_syngraph <- ALG6_syngraph[!(ALG6_syngraph[, 1] %in% wo_old_ALG6[, 1]), ]

ALGs_syngraph <- rbind(wo_old_previous_ALG6, data.frame(V1 = c(ALG6_syngraph[, 1]), V2= 'd6'))
ALGs_syngraph_large_dominant <- rbind(wo_old_ALG6, data.frame(V1 = c(wo_prev_asn_ALG6_syngraph[, 1]), V2= 'd6'))

# write.table(ALGs_syngraph, file = 'tables/ALGs_syngraph_diptera_syngraph_ALG6.tsv', col.names = F, quote = F, row.names = F)
# write.table(ALGs_syngraph_large_dominant, file = 'tables/ALGs_syngraph_diptera_syngraph_ALG6_large_dominant.tsv', col.names = F, quote = F, row.names = F)
```

These 7 genes have conflicting definitions between ALG6 and ALG1-5 runs.

### Bootstraping

Files

```
/lustre/scratch126/tol/teams/jaron/users/se13/diptera_alg/data/results/infer_alg_diptera/syngraph_bootstrapping/results
/lustre/scratch126/tol/teams/jaron/users/se13/diptera_alg/scripts/bsub_light_bootstrap_syngraph.sh
/lustre/scratch126/tol/teams/jaron/users/se13/diptera_alg/data/results/infer_alg_diptera/diptera.pruned.newick

mkdir -p data/syngraph/bootstrap
scp farm:/lustre/scratch126/tol/teams/jaron/users/se13/diptera_alg/data/results/infer_alg_diptera/diptera.pruned.newick data/syngraph/bootstrap/
scp -r farm:/lustre/scratch126/tol/teams/jaron/users/se13/diptera_alg/data/results/infer_alg_diptera/syngraph_bootstrapping/results data/syngraph/bootstrap/
```

I killed the process after downloading ~50 of them. which is enough to explore.

Doublechecking the tree.

```R
library(ape)

bootstrap_tree <- read.tree('data/syngraph/bootstrap/diptera.pruned.newick')
alg_tree <- read.tree('data/syngraph/diptera.no_plecia.mindist.m165.newick.txt')

bootstrap_tree$tip.label == alg_tree$tip.label
# T
bootstrap_tree$node.label == alg_tree$node.label
# F, but the number is the same

head(bootstrap_tree$node.label)
# [1] "0"       "1100.00" "2100.00" "3100.00" "4100.00" "5100.00"
head(alg_tree$node.label)
# [1] ""    "n1"  "n3"  "n7"  "n13" "n21"
```

All good I think. 

Extracing all nodes 1 and 2 assignments.

```bash
ls results/* > replicates
mkdir -p n1n2_nodes

while read file; do
    replicate=$(echo $file | cut -f 2 -d '.')
    echo $replicate; 
    # node
    # for diptera ALGs
    # awk -v node="n1" 'NR==1 { for (i=1; i<=NF; i++){f[$i] = i} }{ if( $(f[node"_seq"]) != "NA" ){ print $1 "\t" $(f[node"_seq"])}}' $file > n1n2_nodes/bootstrap."$replicate".n1.tsv
    # awk -v node="n2" 'NR==1 { for (i=1; i<=NF; i++){f[$i] = i} }{ if( $(f[node"_seq"]) != "NA" ){ print $1 "\t" $(f[node"_seq"])}}' $file > n1n2_nodes/bootstrap."$replicate".n2.tsv
    # for brachycera ALGs (db)
    awk -v node="n21" 'NR==1 { for (i=1; i<=NF; i++){f[$i] = i} }{ if( $(f[node"_seq"]) != "NA" ){ print $1 "\t" $(f[node"_seq"])}}' $file > brachy_nodes/bootstrap."$replicate".n21.tsv
    # for schizophora ALGs (ds)
    # awk -v node="n133" 'NR==1 { for (i=1; i<=NF; i++){f[$i] = i} }{ if( $(f[node"_seq"]) != "NA" ){ print $1 "\t" $(f[node"_seq"])}}' $file > n1n2_nodes/bootstrap."$replicate".n133.tsv
    awk -v node="n12" 'NR==1 { for (i=1; i<=NF; i++){f[$i] = i} }{ if( $(f[node"_seq"]) != "NA" ){ print $1 "\t" $(f[node"_seq"])}}' $file > anoph_nodes/bootstrap."$replicate".n12.tsv

done < replicates
```

Check 
   1. do they ahve the same number of ALGs (% of cases) 

```R
library(pheatmap)

# ALGs_initial <- read.table('tables/ALGs_syngraph_diptera.tsv', col.names = c('busco', 'ALG'))
replicate <- 1


# dir('data/syngraph/bootstrap/n1n2_nodes/', pattern = 'n1.tsv')
# brachy_nodes/
# replicates <- sapply(strsplit(dir('n1n2_nodes/', pattern = 'n1.tsv'), '[.]'), function(x){ as.numeric(x[2]) } )
replicates <- 1:1000
replicate_asn <- list()

loadBrachyAsn <- function(rep){
    brachy_asn <- read.table(paste0('brachy_nodes/bootstrap.', rep, '.n21.tsv'), col.names = c('busco', 'n2'))
    brachy_asn[, 1] <- sapply(strsplit(brachy_asn[, 1], '[.]'), function(x){x[1]})
    return(brachy_asn)
}
replicate_brachy_asn <- lapply(replicates, loadBrachyAsn)

loadAnoAsn <- function(rep){
    ano_asn <- read.table(paste0('anoph_nodes/bootstrap.', rep, '.n12.tsv'), col.names = c('busco', 'n2'))
    ano_asn[, 1] <- sapply(strsplit(ano_asn[, 1], '[.]'), function(x){x[1]})
    return(ano_asn)
}
replicate_ano_asn <- lapply(replicates, loadAnoAsn)

for (replicate in replicates) {
    print(replicate)

    n1_asn <- read.table(paste0('n1n2_nodes/bootstrap.', replicate, '.n1.tsv'), col.names = c('busco', 'n1'))
    n2_asn <- read.table(paste0('n1n2_nodes/bootstrap.', replicate, '.n2.tsv'), col.names = c('busco', 'n2'))

    n1_asn[, 1] <- sapply(strsplit(n1_asn[, 1], '[.]'), function(x){x[1]})
    n2_asn[, 1] <- sapply(strsplit(n2_asn[, 1], '[.]'), function(x){x[1]})

    # any(rowSums(table(n1_asn[, 1], n1_asn[, 2]) > 0) > 1)
    # This means that no artificial replicates ever got different assignments

    asn <- data.frame(busco = unique(c(n1_asn[, 1], n2_asn[, 1])), 'n1' = 'unasn', 'n2' = 'unasn')
    row.names(asn) <- asn[, 1]

    asn[n1_asn[, 1], 'n1'] <- n1_asn[, 2]
    asn[n2_asn[, 1], 'n2'] <- n2_asn[, 2]

    # head(asn)

    # table(asn[, 'n1'], asn[, 'n2'])

    syngraph <- asn[asn[, 'n1'] != 'unasn' & asn[, 'n2'] != 'unasn', ]
    syngraph[, 'ALG'] <- paste(syngraph[, 'n1'], syngraph[, 'n2'], sep = '_')
    replicate_asn[[replicate]] <- syngraph[, c(1, 4)]

    # syngraph[, 'ALG'] <- NA

    # # # sorted table shows 5 corresponding ALGs + several genes that moved around, those we can ignore
    # table_of_groups <- sort(table(paste(syngraph[, 'n1'], syngraph[, 'n2'])), T) # 4 LGs, not 5
    # # print(table_of_groups)

    # if ((table_of_groups / sum(table_of_groups))[5] < 0.05 | (table_of_groups / sum(table_of_groups))[6] > 0.01){
    #     print("WARNING (too few, or too many chromosomes)")
    #     print(table_of_groups)
    # } else {
    #     dominant_alg_pairs <- strsplit(names(table_of_groups)[1:5], ' ')

    #     alg_groups <- list()
    #     for (ALG in 1:5){
    #         n1_name <- dominant_alg_pairs[[ALG]][1]
    #         n2_name <- dominant_alg_pairs[[ALG]][2]
    #         alg_rows <- which(syngraph[, 'n1'] == n1_name & syngraph[, 'n2'] == n2_name)

    #         alg_groups[[ALG]] <- syngraph[alg_rows, 'busco']
    #     }

    #     for (ALG in 1:5){
    #         alg_label <- names(sort(table(ALGs_initial[ALGs_initial[, 'busco'] %in% alg_groups[[ALG]], 'ALG']), T)[1])
            
    #         syngraph[alg_groups[[ALG]], 'ALG'] <- alg_label
    #         # print(paste(names(table_of_groups)[ALG], "assigned as", alg_label, "with", length(alg_groups[[ALG]]), "marker genes"))
    #     }

    #     write.table(syngraph[, c(1, 4)], file = paste0('data/syngraph/bootstrap/n1n2_nodes/bootstrap', replicate, '_overview.tsv'))
    # }
}

all_buscos <- unique(unlist(as.vector(sapply(replicate_brachy_asn, function(x){x[, 1]}))))
all_buscos <- unique(unlist(as.vector(sapply(replicate_ano_asn, function(x){x[, 1]}))))

# alg_matrix <- matrix(0, nrow = length(all_buscos), ncol = length(replicate_brachy_asn))
alg_matrix <- matrix(0, nrow = length(all_buscos), ncol = length(replicate_ano_asn))
row.names(alg_matrix) <- all_buscos

for ( i in 1:length(replicates)){
    print(i)
    alg_matrix[replicate_asn[[i]][, 'busco'], i] <- replicate_asn[[i]][, 'ALG']
}

for ( i in 1:length(replicate_brachy_asn)){
    print(i)
    alg_matrix[replicate_brachy_asn[[i]][, 1], i] <- replicate_brachy_asn[[i]][, 2]
}

for ( i in 1:length(replicate_ano_asn)){
    print(i)
    alg_matrix[replicate_ano_asn[[i]][, 1], i] <- replicate_ano_asn[[i]][, 2]
}

busco_similarity_matrix <- matrix(0, nrow = length(all_buscos), ncol = length(all_buscos))
number_of_buscos <- length(all_buscos)
for (b1i in 1:length(all_buscos)){
    print(b1i)
    b1 <- all_buscos[b1i]
    colocalised_buscos <- sapply(1:number_of_buscos, function(x){ sum(alg_matrix[b1, ] == alg_matrix[all_buscos[x], ] & alg_matrix[b1, ] != 0) } )
    busco_similarity_matrix[b1i, ] <- colocalised_buscos # makes it a distance
}

busco_similarity <- as.data.frame(busco_similarity_matrix)
colnames(busco_similarity) <- all_buscos
write.table(busco_similarity, 'busco_bootstrap_coocurence_anophe.tsv', col.names = T, row.names = F)
# busco_similarity <- read.table('tables/busco_bootstrap_coocurence_wo0.tsv')
# busco_similarity_matrix <- as.matrix(busco_similarity)

# , annotation_col = all_buscos
row.names(busco_similarity_matrix) <- all_buscos
pdf('busco_pheatmap_annot_annopth.pdf', width = 300, height = 300)
    pheatmap(busco_similarity_matrix)
dev.off()


# plot test
# busco_similarity_matrix_test <- matrix(runif(16, 0, 1), nrow = 4000, ncol = 4000)
# row.names(busco_similarity_matrix_test) <- paste0('gene', 1:4000)

# pdf('test.pdf', width = 500, height = 500)
#     pheatmap(busco_similarity_matrix_test)
# dev.off()

bootstrap <- read.table('tables/busco_bootstrap_coocurence_brachycera.tsv', header = T)
all_buscos <- sapply(strsplit(colnames(bootstrap), 'X'), function(x){ x[2]})
colnames(bootstrap) <- all_buscos

hist(bootstrap[, '206736at7147'], breaks = 10, main = '206736at7147', xlab = 'coocurance', ylab = 'number of BUSCOs')
hist(bootstrap[, '164810at7147'], breaks = 10, main = '164810at7147', xlab = 'coocurance', ylab = 'number of BUSCOs')
hist(bootstrap[, '60829at7147'], breaks = 10, main = '60829at7147', xlab = 'coocurance', ylab = 'number of BUSCOs')
hist(bootstrap[, '75965at7147'], breaks = 10, main = '75965at7147', xlab = 'coocurance', ylab = 'number of BUSCOs')
hist(bootstrap[, '61998at7147'], breaks = 10, main = '61998at7147', xlab = 'coocurance', ylab = 'number of BUSCOs')
hist(bootstrap[, '173520at7147'], breaks = 10, main = '173520at7147', xlab = 'coocurance', ylab = 'number of BUSCOs')
hist(bootstrap[, '245796at7147'], breaks = 10, main = '245796at7147', xlab = 'coocurance', ylab = 'number of BUSCOs')
# all these look solid!



BUSCOs <- length(all_buscos) # 4830 #4526
bootstrap_p <- function(n){ (1 - (1 / n))^n }
coocurence_p <- (1 - bootstrap_p(BUSCOs))^2
qbinom(1e-6, 1000, coocurence_p)

hist(unlist(alg1), breaks = 300, col = 'red')
hist(rbinom(length(alg1)^2, 1000, coocurence_p), add = T, breaks = 150)

assigned_to_an_alg_per_busco <- diag(as.matrix(bootstrap))
assigned_to_an_alg_per_busco <- c(rep(0, length = BUSCOs - length(assigned_to_an_alg_per_busco)), assigned_to_an_alg_per_busco)

pdf('figures/bootstrap/BUSCOs_assigned_in_an_ALG_Brachy.pdf')
    hist(assigned_to_an_alg_per_busco, breaks = 60, main = 'Frequency of BUSCOs assigned in an ALG', ylab = 'BUSCOs', xlab = 'Bootstrap replicate ALG asignment')
    lines(c((1 - bootstrap_p(BUSCOs)) * 1000, (1 - bootstrap_p(BUSCOs)) * 1000), c(0, 1000), lty = 3, col = 'orange', lwd = 3)
    legend('topleft', lty = 3, col = 'orange', lwd = 3, 'sampling expecation if assigned every time', bty = 'n')
dev.off()

frequenly_unconverging <- assigned_to_an_alg_per_busco < 327

reminding_ALGs <- bootstrap[!frequenly_unconverging, !frequenly_unconverging]
ALGs <- list()
i <- 1
while(nrow(reminding_ALGs) > 0){
    alg_subset <- reminding_ALGs[, 1] > 327
    if(sum(alg_subset) == 1){
        reminding_ALGs <- reminding_ALGs[!alg_subset, !alg_subset]
        next
    }
    ALGs[[i]] <- reminding_ALGs[alg_subset, alg_subset]
    print(paste("ALG:", i, "of", nrow(ALGs[[i]]), "markers"))
    print(paste("Range of co-occurances within of the group:", min(unlist(ALGs[[i]])), max(unlist(ALGs[[i]]))))
    print(paste("Maximal co-occurance outside of the group:", max(unlist(reminding_ALGs[alg_subset, !alg_subset]))))

    reminding_ALGs <- reminding_ALGs[!alg_subset, !alg_subset]
    i <- i + 1
}

hist(bootstrap[, '261225at7147'])
# alwyas with the same group of genes, but less often

bootstrap[bootstrap[, '15074at7147'] > 1, bootstrap[, '15074at7147'] > 1]
# strange, this one is never with any other gene?! how does that work? Ha, that will be a unique gene that is in one LG in culicomprpha and a different LG in the other branch

hist(unlist(alg2), breaks = 300, col = 'red')
sapply(ALGs, length)
# 493  709  840    2   16 1054  644   14   12   15   14    3
ALGs <- ALGs[1:3]

genes_in_ALGs <- lapply(ALGs, colnames)
instable_asn_buscos <- unlist(genes_in_ALGs[c(4, 5, 8:12)])
bigALGs <- genes_in_ALGs[c(1:3, 6:7)]

reminding_ALGs[, 1]

bigALGs_distance <- ALGs[c(1:3, 6:7)]
for (i in 1:5){
    diag(bigALGs_distance[[i]]) <- NA
    bigALGs_distance[[i]][upper.tri(bigALGs_distance[[i]])] <- NA
    print(range(bigALGs_distance[[i]], na.rm = T))
}


original_algs <- read.table('../../../tables/ALGs_syngraph_diptera.tsv')
original_algs <- read.table('../../../tables/ALGs_syngraph_brachycera.tsv')
row.names(original_algs) <- original_algs[, 1]

instable_asn_buscos[instable_asn_buscos %in% original_algs[, 1]]
original_algs[c("181571at7147", "19395at7147"), ]
# these three get removed from alg3

head(genes_in_ALGs)

ALG_labels <- sapply(1:9, function(x){table(original_algs[bigALGs[[x]], 2])})

names(bigALGs) <- names(ALG_labels)

bootstrap_ALGs <- do.call("rbind", lapply(1:9, function(x){ data.frame(busco = bigALGs[[x]], names(bigALGs)[x]) } ))

row.names(bootstrap_ALGs) <- bootstrap_ALGs[, 1]

ALG6_core_conflicts <- row.names(bootstrap_ALGs) %in% c('162912at7147', '164810at7147', '60829at7147', '75965at7147', '61998at7147', '173520at7147', '245796at7147') 

bootstrap_ALGs_cleaned <- bootstrap_ALGs[!ALG6_core_conflicts, ]
colnames(bootstrap_ALGs_cleaned) <- c('busco', 'alg')

ALG6_syngraph <- read.table('../../../data/syngraph/alg6/nodes/n1_asgn.tsv')
ALG6_syngraph_brachy <- read.table('../../../data/syngraph/alg6/nodes/n5_asgn.tsv')

ALG6_syngraph_brachy[, 1] %in% bootstrap_ALGs[, 1]

core_ALG6_conflicts <- ALG6_syngraph[, 1] %in% c('162912at7147', '164810at7147', '60829at7147', '75965at7147', '61998at7147', '173520at7147', '245796at7147')

colnames(bootstrap_ALGs) <- c('busco', 'alg')
final_ALGs <- rbind(bootstrap_ALGs, data.frame(busco = ALG6_syngraph_brachy[, 1], alg = 'db6'))
table(final_ALGs[, 2])
final_ALGs <- final_ALGs[order(final_ALGs[, 2]), ]
table(original_algs[, 2])
head(final_ALGs)

rownames(final_ALGs) <- final_ALGs[, 1]
buscos_in <- original_algs[!original_algs[, 1] %in% final_ALGs[, 1], ]
buscos_out <- final_ALGs[!final_ALGs[, 1] %in% original_algs[, 1], ]
final_ALGs <- final_ALGs[c(buscos_in, buscos_out), ]
head(final_ALGs)

write.table(final_ALGs, file = 'ALGs_brachy_bootstrapping.tsv', col.names = F, quote = F, row.names = F)

# 
library(RColorBrewer)
cols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(32)

pdf('figures/bootstrap/legend.pdf')
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), bty = 'n', axes = F, xlab = '', ylab = '')
for(i in 1:32){ 
    rect(0, 0 + (i - 1) * 1/32, 0.4, i * 1/32, col = cols[i], border = NA)
}
text(0.5, c(0, 200, 400, 600) / 687, labels = c(0, 200, 400, 600), cex = 2)
rect(0, 0 + 327 / 687, 0.4, 474 / 687, col = NA, lwd = 2)
 
dev.off()

```


   2. generate per-replicate ALG assignment (using n1 and n2 overlap);

Generate boostrapping heatmap - how stable is co-occurance of markers.

```
convert           \
   -verbose       \
   -density 30   \
   -trim          \
    busco_pheatmap.pdf      \
   -quality 20   \
   -flatten       \
   -sharpen 0x1.0 \
    busco_heatmap2.png

magick busco_pheatmap_annot2.pdf busco_pheatmap_annot2.png

busco_pheatmap_annot2.pdf

convert -verbose 1.pdf -resize 100% -quality 100 -flatten -sharpen 0x1.0 page-%03d.png
```

```R
BUSCOs <- 5067 #4526
bootstrap_p <- function(n){ (1 - (1 / n))^n }
coocurence_p <- (1 - bootstrap_p(BUSCOs))^2

hist(rbinom(3000, 1000, coocurence_p))

qbinom(c(1e-6, 1 - 1e-6), 1000, coocurence_p)
# 364 436
```

### family paints

```bash
line_index=1

while IFS= read -r line || [[ -n "$line" ]]; do
    echo $line_index $line

    if [ "$line_index" -gt 46 ]; then
        # 47 is the first nematoceran
        Rscript scripts/misc/ALG_painter_ME_building_blocks.R -f $line -a tables/ALGs_syngraph_diptera.tsv -o figures/paints/clean/"$line_index"_"$line" --tree-order --alg-sort
    else
        Rscript scripts/misc/ALG_painter_ME_building_blocks.R -f $line -a tables/ALGs_syngraph_brachycera.tsv -o figures/paints/clean/"$line_index"_"$line" --tree-order --alg-sort
    fi
        
    ((line_index++))
done < tables/families_ordered_by_tree.tsv


```