suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(gsheet))
suppressPackageStartupMessages(library(ape))

parser <- ArgumentParser()
parser$add_argument("-t", "--tree",  
    dest="t", help="Path to the tree that was used for the syngraph run")

args <- parser$parse_args()

dmel_busco <- read.table('data/busco_tables/Drosophila_melanogaster.syngraph.buscos.tsv', header = F)
# syngraph_tree <- read.tree('data/syngraph/diptera.pruned.syngraph_infer.newick.txt')
syngraph_tree <- read.tree(args$t) # 'data/syngraph/diptera.no_plecia.mindist.m165.newick.txt'

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
# 0.741433; 238

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

pdf('figures/BUSCOs_on_the_dot.pdf')
plot(as.numeric(head(species_with_dot_buscos, 130) / length(buscos_on_dots)), xlab = 'BUSCOs ordered by % of dots carrying the gene', ylab = 'Shared by # species', pch = 20)
lines(c(0, 200), c(0.4, 0.4), lty = 2)
dev.off()

ALG6 <- data.frame(V1 = names(species_with_dot_buscos)[(species_with_dot_buscos / length(buscos_on_dots)) > 0.4], V2 = 'd6')

# ALGs <- read.table('data/syngraph/syngraph.pruned.100.ALGs.inferred.tsv')

# any(ALG6[, 1] %in% ALGs[, 1])

# all_ALGs <- rbind(ALGs, ALG6)
write.table(ALG6, file = 'data/ALG6_BUSCOs.tsv', quote = F, sep = '\t', col.names = F, row.names = F)