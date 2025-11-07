suppressPackageStartupMessages(library('gsheet'))
suppressPackageStartupMessages(library('argparse'))

source('scripts/20250620_colour_pal.R')


parser <- ArgumentParser()
parser$add_argument("-o", "--output",  
    dest="o", help="Base of the output name (_all.pdf and _nematocera.pdf will be attached)")
parser$add_argument("-a", "--alg-set",  
    help="File with BUSCO to ALG assignments", default = "tables/ALGs_syngraph_diptera.tsv")

args <- parser$parse_args()

all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM", format='csv'),
                            stringsAsFactors = F, header = T, check.names = F)
all_genome_data <- all_genome_data[all_genome_data[, 'TO ADD'] %in% c('KEEP'), ]

species2suborder <- all_genome_data[, c('species', 'suborder')]
row.names(species2suborder) <- species2suborder[, 'species']

busco_asn_file <- args$a #'tables/ALGs_syngraph_diptera.tsv'
# read

ALG_data <- read.table(busco_asn_file, col.names = c('busco', 'ALG'))
ALG_sizes <- table(ALG_data[, 'ALG'])
row.names(ALG_data) <- ALG_data[, 'busco']

# make matrix of all BUSCOs

# busco_files <- dir('data/busco_tables/')

process_species_data <- function(species_name){
    busco_on_chr <- read.table(paste0('data/busco_tables/', species_name, '.syngraph.buscos.tsv'), col.names = c('busco', 'chr', 'f', 't'))[, c(1,2)]
    busco_on_chr[, 'alg'] <- ALG_data[busco_on_chr[, 'busco'], 'ALG']

    chr_busco_matrix <- table(busco_on_chr[, 'alg'], busco_on_chr[, 'chr'])
    apply(chr_busco_matrix, 1, max) / ALG_sizes
}

all_species_matrix <- sapply(species2suborder[, 1], process_species_data)
ALGs <- row.names(all_species_matrix)

# 'c(bottom, left, top, right)'
#           which gives the number of lines of margin to be specified on
#           the four sides of the plot.  The default is 'c(5, 4, 4, 2) +
#           0.1'.

if (length(ALGs) == 6){
    grid_type <- c(2, 3)
} else {
    grid_type <- c(3, 4)
}
cex = 1.4

pdf(paste0(args$o, '_all.pdf'))
# pdf('figures/ALG_stability_histograms_all_species.pdf')
    par(mfrow = grid_type, mar = c(3,3,2,1))

    for (i in 1:nrow(all_species_matrix)){
        hist(all_species_matrix[i, ], breaks = 40, main = NA, col = pal[ALGs[i]], xlim = c(0,1), border = NA, ylab = NA, xlab = NA, cex.axis = cex)
        legend('topleft', ALGs[i], bty = 'n', cex = cex, text.font = 2)
    }

dev.off()

pdf(paste0(args$o, '_nematocera.pdf'))
# pdf('figures/ALG_stability_histograms_Nematocera.pdf')
    par(mfrow = grid_type, mar = c(3,3,2,1))

    for (i in 1:nrow(all_species_matrix)){
        hist(all_species_matrix[i, species2suborder[, 'suborder'] == 'Nematocera'], breaks = 20, main = NA, col = pal[ALGs[i]], xlim = c(0,1), border = NA, ylab = NA, xlab = NA, cex.axis = cex)
        legend('topleft', ALGs[i], bty = 'n', cex = cex, text.font = 2)
    }

dev.off()

pdf(paste0(args$o, '_brachycera.pdf'))
# pdf('figures/ALG_stability_histograms_Nematocera.pdf')
    par(mfrow = grid_type, mar = c(3,3,2,1))

    for (i in 1:nrow(all_species_matrix)){
        hist(all_species_matrix[i, species2suborder[, 'suborder'] == 'Brachycera'], breaks = 40, main = NA, col = pal[ALGs[i]], xlim = c(0,1), border = NA, ylab = NA, xlab = NA, cex.axis = cex)
        legend('topleft', ALGs[i], bty = 'n', cex = cex, text.font = 2)
    }

dev.off()

# xlab = paste0('Highest proportion of the ALG ', i, ' on a single chromosome'), 

# busco_per_chr_list <- split(busco_on_chr[, 1], busco_on_chr[, 'chr'])
# ALG_data[busco_per_chr_list[[1]], 'ALG']


# length(busco_per_chr_list)



# what is the set of genomes we want do analyse here?

# read the genomes one by one
#        mark coocurance of genes onto matrix

