suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(gsheet))
suppressPackageStartupMessages(library(ape))

parser <- ArgumentParser()
parser$add_argument("-p", "--phylo-tree",  
    dest="p", help="Path to the tree that was used for the generating phylogeny")

args <- parser$parse_args()

# phylo_tree_file <- args$p
phylo_tree_file <- 'data/diptera.supermatrix.phy.treefile'
phylo_tree <- read.tree(phylo_tree_file)
print('Are all the files present?')
all(paste0(phylo_tree$tip.label, ".syngraph.buscos.tsv")  %in% dir('data/busco_tables/'))

buscos_of_all_analysed <- paste0("data/busco_tables/", phylo_tree$tip.label, ".syngraph.buscos.tsv") 
# dot_overview will be my final table
dot_overview <- data.frame(species = phylo_tree$tip.label)

ALGs <- read.table('data/diptera.no_plecia.mindist.m165_n1_n2.tsv')
ALG6_genes <- ALGs[ALGs[, 2] == 'd6', 1]

# read.table(buscos_of_all_analysed[1])

getTheDot <- function(busco_file){
    busco_tab <- read.table(busco_file, header = F)
    asn_busco_tab <- busco_tab[busco_tab[, 1] %in% ALGs[, 1], ]
    chromosomes_tab <- data.frame(table(asn_busco_tab[, 2]), ALG6_buscos = 0)
    
    colnames(chromosomes_tab) <- c('chr', 'buscos', 'ALG6_buscos')
    rownames(chromosomes_tab) <- chromosomes_tab[, 1]

    just_dot <- busco_tab[busco_tab[, 1] %in% ALG6_genes, ]

    if (nrow(just_dot) != 0 ){
        dot_buscos_per_ch <- table(just_dot[, 2])
        chromosomes_tab[names(dot_buscos_per_ch), 'ALG6_buscos'] <- as.numeric(dot_buscos_per_ch)
        chromosomes_tab[, 'dot_frac'] <- round(chromosomes_tab[, 'ALG6_buscos'] / chromosomes_tab[, 'buscos'], 4)
        highest_fract <- max(chromosomes_tab[,'dot_frac'])
        names(highest_fract) <- chromosomes_tab[which.max(chromosomes_tab[,'dot_frac']), 'chr']
        # buscos <- as.numeric(busco_per_ch[names(which.max(dot_buscos_per_ch))])
        # dot_buscos <- max(dot_buscos_per_ch)
        if (highest_fract > 0.5){
            return(names(highest_fract))
        }
    }
    return(NA)
}

dots <- sapply(buscos_of_all_analysed, getTheDot)

dot_overview[, 'ALG6_chr'] <- as.vector(dots)
# 243 has a dot
# 71.26%

all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?gid=1940964825#gid=1940964825", format='csv'),
                            stringsAsFactors = F, header = T, check.names = F)
row.names(all_genome_data) <- all_genome_data[, 'chromosome']

assigned_dots <- which(!is.na(dot_overview[, 2]))
dot_overview[assigned_dots, 'size'] <- all_genome_data[dot_overview[assigned_dots, 2], 'chromsome_size_b']
dot_overview[assigned_dots, 'buscos'] <- all_genome_data[dot_overview[assigned_dots, 2], 'busco_odb12_complete_count']

write.table(dot_overview, file = 'data/species_and_ALG6_chromosomes.tsv', quote = F, sep = '\t', col.names = T, row.names = F)