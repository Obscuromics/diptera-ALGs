


```R
require(gsheet)
require(ape)
dmel_busco <- read.table('data/busco_tables/Drosophila_melanogaster.syngraph.buscos.tsv', header = F)
syngraph_tree <- read.tree('data/syngraph/diptera.pruned.syngraph_infer.newick.txt')
# 336 tips

table(dmel_busco[, 2])
# NC_004353.4 is the dot

dot_buscos <- dmel_busco[dmel_busco[, 2] == 'NC_004353.4', 1]

all_busco_files <- dir('data/busco_tables/')

# all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?usp=sharing", format='csv'),
#                             stringsAsFactors = F, header = T, check.names = F)
# all_genome_data <- all_genome_data[all_genome_data[, 'TO ADD'] == 'KEEP', ]

all_busco_files %in% paste0(syngraph_tree$tip.label, ".syngraph.buscos.tsv") 

files_used_for_syngraph <- paste0(syngraph_tree$tip.label, ".syngraph.buscos.tsv") 

busco_file <- files_used_for_syngraph[[4]]

# getTheSizeOfTheMostDottedCh <- function(busco_file){
#     busco_tab <- read.table(paste0('data/busco_tables/', busco_file), header = F)
#     busco_per_ch <- table(busco_tab[, 2])

#     just_dot <- busco_tab[busco_tab[, 1] %in% dot_buscos, ]
#     if (nrow(just_dot) == 0 ){
#         chr_name = NA
#         buscos = 0
#         dot_buscos = 0
#     } else {
#         dot_buscos_per_ch <- table(just_dot[, 2])
#         chr_name <- names(which.max(dot_buscos_per_ch))
#         buscos <- as.numeric(busco_per_ch[names(which.max(dot_buscos_per_ch))])
#         dot_buscos <- max(dot_buscos_per_ch)
#     }
#     return(data.frame('chr_name' = chr_name, 'buscos' = buscos, 'dot_buscos' = dot_buscos))
# }

# dot_sizes <- sapply(files_used_for_syngraph, getTheSizeOfTheMostDottedCh)
# dot_tab <- t(dot_sizes)

# hist(as.numeric(dot_tab[, 'buscos']), breaks = 80)
# hist(as.numeric(dot_tab[, 'dot_buscos']), breaks = 20)

# plot(as.numeric(dot_tab[, 'buscos']) ~ as.numeric(dot_tab[, 'dot_buscos']), ylim = c(0, 300))

getTheDot <- function(busco_file){
    busco_tab <- read.table(paste0('data/busco_tables/', busco_file), header = F)
    busco_per_ch <- table(busco_tab[, 2])

    just_dot <- busco_tab[busco_tab[, 1] %in% dot_buscos, ]
    if (nrow(just_dot) != 0 ){
        dot_buscos_per_ch <- table(just_dot[, 2])
        chr_name <- names(which.max(dot_buscos_per_ch))
        buscos <- as.numeric(busco_per_ch[names(which.max(dot_buscos_per_ch))])
        dot_buscos <- max(dot_buscos_per_ch)
        if (buscos < 150){
            return(chr_name)
        }
    }
    return(NA)
}

dot_chromosomes <- sapply(files_used_for_syngraph, getTheDot)
dot_overview <- data.frame(species = names(dot_chromosomes), dot = as.character(dot_chromosomes))

mean(!is.na(dot_overview[, 2]))
# 0.7172619; 241

with_a_dot <- dot_overview[!is.na(dot_overview[, 2]), ]
row.names(with_a_dot) <- with_a_dot[, 1]

busco_file <- with_a_dot[1, 1]
get_buscos_of_a_chromosome <- function(busco_file){
    busco_tab <- read.table(paste0('data/busco_tables/', busco_file), header = F)
    the_dot <- with_a_dot[busco_file, 'dot']
    busco_tab[busco_tab[, 2] == the_dot, 1]
}

buscos_on_dots <- lapply(with_a_dot[, 1], get_buscos_of_a_chromosome)

species_with_dot_buscos <- table(unlist(buscos_on_dots))
species_with_dot_buscos <- sort(species_with_dot_buscos, decreasing = T)
# names(species_with_dot_buscos) <- 1:407

barplot(head(species_with_dot_buscos, 130) / 241, xlab = 'BUSCOs ordered', ylab = 'Shared by # species')

ALG6 <- data.frame(V1 = names(species_with_dot_buscos)[(species_with_dot_buscos / 241) > 0.4], V2 = 'd6')

ALGs <- read.table('data/syngraph/syngraph.pruned.100.ALGs.inferred.tsv')

any(ALG6[, 1] %in% ALGs[, 1])

all_ALGs <- rbind(ALGs, ALG6)
write.table(all_ALGs, file = 'data/ALG_assignments_pruned_m100.tsv', quote = F, sep = '\t', col.names = F, row.names = F)
```

### Add the previsouly inferred ALG6 to newer assignments

```R
ALGs_pruned <- read.table('data/ALG_assignments_pruned_m100.tsv', col.names = c('busco', 'alg_pruned'))
ALGs1to5_pruned2 <- read.table('data/syngraph.pruned2.100.ALGs.inferred.tsv', col.names = c('busco', 'alg_pruned2'))
# is.na(ALGs1to5_pruned2[, 2])

dim(ALGs1to5_pruned2)
dim(ALGs_pruned)

both <- merge(ALGs_pruned, ALGs1to5_pruned2)
table(paste(both[, 2], both[, 3]))

ALG6_buscos <- ALGs_pruned[ALGs_pruned[, 'alg_pruned'] == 'd6', 'busco']

any(ALG6_buscos %in% ALGs1to5_pruned2[, 'busco'])
# all good, none of them is assigned already

ALGs_incl_the_dot_pruned2 <- rbind(ALGs1to5_pruned2, data.frame('busco' = ALG6_buscos, 'alg_pruned2' = 'd6'))

ordered_ALGs <- ALGs_incl_the_dot_pruned2[order(ALGs_incl_the_dot_pruned2[, 2]), ]

write.table(ordered_ALGs, file = 'data/ALG_assignments_pruned2_m100.tsv', col.names = F, row.names = F, quote = F, sep = '\t')

```