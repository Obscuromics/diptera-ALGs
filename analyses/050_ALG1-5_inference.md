### Running syngraph on the cluster

It runs fast, but it's more than what should be on the headnode, therefore it's better to do it in an interactive job.

```bash
bsub -G team360 -Is -n 16 -M 10240 -R "select[mem>10240] rusage[mem=10240]" bash -l
```

#### Sam's instructions for modifying syngraph run

```bash
cd /data/tol/teams/jaron/lustre/users/se13/diptera_alg
module load syngraph/0.0.1-c3
```
You need biopython somehow too

The included busco files are in:

```bash
./data/results/infer_alg_diptera/syngraph_busco_tables/
```

If you want to change the included files modify them in place, or to exclude them move them to:

```bash
./data/results/infer_alg_diptera/excluded_syngraph_busco_tables/
```
and run

```bash
./scripts/prep_syngraph_tree.py 
```

which will create a new tree, modifying the iqtree output by removing the leaves which have files in /excluded_syngraph_busco_tables. Theres also a line you can change from 'None' to the outgroup clade file if you want and then it will reroot the tree.

Then you just have to rerun syngraph

```bash
syngraph build -d data/results/infer_alg_diptera/syngraph_busco_tables -m -o data/results/infer_alg_diptera/syngraph/diptera.pruned.syngraph_build

syngraph infer -g data/results/infer_alg_diptera/syngraph/diptera.pruned.syngraph_build.pickle -t data/results/infer_alg_diptera/diptera.pruned.newick -m 89 -r 2 -a quick -s Phlebotomus_papatasi -o data/results/infer_alg_diptera/syngraph/diptera.pruned.syngraph_infer

syngraph tabulate -g data/results/infer_alg_diptera/syngraph/diptera.pruned2.syngraph_infer.with_ancestors.pickle -o data/results/infer_alg_diptera/syngraph/diptera.pruned.syngraph_tabulate
```
If you run this as is the output files will overwrite the last analysis I did so could edit those.

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

### generalising

I will try to write it in the way that if we have to redo this, I will be able to regenerate all the things...

```bash
SYNGTAB=data/syngraph/diptera.pruned3.syngraph_tabulate.table.tsv
SYNGNWK=data/syngraph/diptera.pruned3.syngraph_infer.newick.txt

rm data/syngraph/node_assignments/* # clean previous node assignments
Rscript scripts/20250718_tree2internal_nodes.R -t $SYNGNWK -o data/syngraph/all
# this also plots the tree as "figures/syngraph_tree_with_nodes.pdf"

while read node; do 
    echo $node; 
    # node
    awk -v node="$node" 'NR==1 { for (i=1; i<=NF; i++){f[$i] = i} }{ if( $(f[node"_seq"]) != "NA" ){ print $1 "\t" $(f[node"_seq"])}}' $SYNGTAB > data/syngraph/node_assignments/"$node"_asgn.tsv
done < data/syngraph/all_internal_nodes.tsv
```

and now using `n1` and `n2` nodes, define ALGs1 - 5. The ALG 6 is added too, but that is described in the following .md file.

```bash
Rscript scripts/20250719_ALG_definition.R -o data/syngraph.pruned3.100.ALGs.inferred.tsv
```


### Painting of internal nodes

```bash
ALGs=data/ALG_assignments_pruned2_m100.tsv
while read node; do 
    echo $node
    Rscript scripts/plot_node.R -n $node -a $ALGs -o data/syngraph/node_plots/$node
done < data/syngraph/all_internal_nodes.tsv
```