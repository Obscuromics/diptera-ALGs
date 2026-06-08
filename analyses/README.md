There are various analyses documents. I will most likely collide them into one long README in the end. This is just to have it more streamlined for the moment.

The prefixes are made only to order them. Large numbers are chosen to allow for things to be inserted before, between or after.



```R
Belo <- read.table('/Users/kj11/Projects/curation/idBacElon/busco_reduced.tsv')
ALGs <- read.table('tables/ALGs_syngraph_brachycera.tsv', row.names = 1)

Belo[, 'ALG'] <- ALGs[Belo[, 1], 1]

data.frame(scf = names(table(Belo[Belo[, 'ALG'] == 'db6', 3])), ALG6 = as.numeric(table(Belo[Belo[, 'ALG'] == 'db6', 3])))



```

## Muscidae

```R
muscids <- read.table('table.tsv', sep = '\t', header = T)

buscos <- lapply(paste0('BUSCO_runs/', muscids[, 1], '.tsv'), read.table, col.names = c('busco', 'chr', 'from', 'to'))

scf_buscos <- lapply(buscos, function(x){table(x[, 2])})
names(scf_buscos) <- muscids[, 1]

head(buscos[[1]][, 'chr'])

x <- buscos[[1]]

ALG6_buscos <- lapply(buscos, function(x){x[grepl('SUPER_6', x[, 'chr']), 1]})

unplaced_buscos <- lapply(buscos, function(x){x[!grepl('SUPER', x[, 'chr']), 1]})

names(unplaced_buscos) <- muscids[, 1]

muscids[5, ]


ALGs <- read.table('/Users/kj11/Projects/diptera-ALGs/tables/ALGs_syngraph_brachycera.tsv', row.names = 1)

ALGs["322552at7147", ]
```