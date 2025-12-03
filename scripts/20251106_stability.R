suppressPackageStartupMessages(library('gsheet'))
suppressPackageStartupMessages(library('argparse'))
suppressPackageStartupMessages(library('ape'))

source('scripts/20250620_colour_pal.R')


parser <- ArgumentParser()
parser$add_argument("-o", "--output",  
    dest="o", help="Base of the output name (_all.pdf and _nematocera.pdf will be attached)")
parser$add_argument("-a", "--alg-set",  
    help="File with BUSCO to ALG assignments", default = "tables/ALGs_syngraph_diptera.tsv")

args <- parser$parse_args()
# args$o <- "figures/ALG_stability_histograms"
# args$a <- "tables/ALGs_syngraph_diptera.tsv"

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
taxonomy_colors <- c("Nematocera" = "#67d7dbff", "Brachycera" = "#069c24ff", "Schizophora" = "#bb3fa0ff")
colnames(all_genome_data)

nematocera_species <- all_genome_data[all_genome_data[, 'suborder'] == "Nematocera", 'species']
tree <- read.tree("data/syngraph/diptera.no_plecia.mindist.m165.newick.txt")
# 'n133'
# get all descendants of n133
schizophora_species <- extract.clade(tree, node = "n133")$tip.label

brachycera_species <- all_genome_data[all_genome_data[, 'suborder'] == "Brachycera", 'species']
brachycera_species <- brachycera_species[!brachycera_species %in% schizophora_species]

# length(nematocera_species) + length(schizophora_species) + length(brachycera_species)
# 340 
all_species_matrix <- all_species_matrix[, 1:340] # removing Panorpa

nematocera_species_matrix <- all_species_matrix[, nematocera_species]
species_list <- list(nematocera_species, brachycera_species, schizophora_species)

pdf(paste0(args$o, '_all.pdf'), width = 12, height = 8)
# pdf('figures/ALG_stability_histograms_all_species.pdf')
    par(mfrow = grid_type, mar = c(3,3,2,1))

    for (i in 1:nrow(all_species_matrix)){
        # hist(all_species_matrix[i, ], breaks = 40, main = NA, col = pal[ALGs[i]], xlim = c(0,1), border = NA, ylab = NA, xlab = NA, cex.axis = cex)
        # fixed_bin_histogram(all_species_matrix[i, ], breaks = 40, main = NA, xlim = c(0,1), border = NA, ylab = NA, xlab = NA, cex.axis = cex)
        plot(density(all_species_matrix[i, schizophora_species],na.rm=T), col = 'black', lwd = 4, lty = 1, type = 'l', main = '', xlab = NA, ylab = NA, xlim = c(0,1))
        lines(density(all_species_matrix[i, schizophora_species],na.rm=T), col = taxonomy_colors["Schizophora"], lwd = 3, lty = 1)
        
        lines(density(all_species_matrix[i, nematocera_species],na.rm=T), col = 'black', lwd = 4, lty = 1)
        lines(density(all_species_matrix[i, nematocera_species],na.rm=T), col = taxonomy_colors["Nematocera"], lwd = 3, lty = 1)

        lines(density(all_species_matrix[i, brachycera_species],na.rm=T), col = 'black', lwd = 4, lty = 1)
        lines(density(all_species_matrix[i, brachycera_species],na.rm=T), col = taxonomy_colors["Brachycera"], lwd = 3, lty = 1)
        

        legend('topleft', paste('ALG', i), bty = 'n', cex = cex, text.font = 2, fill = pal[ALGs[i]])
        if (i == 3){
            legend('topright', legend = c('Nematocera', 'Brachycera', 'Schizophora'), col = taxonomy_colors, lwd = 2, bty = 'n', cex = cex, lty = 1)
        }
    }

dev.off()

xlim <- c(0,1)
bins <- 40
ax <- pretty(xlim, n = bins) # Make a neat vector for the breakpoints
# histograms[[3]]$counts <- histograms[[3]]$counts

pdf(paste0(args$o, '_stacked_plot.pdf'), width = 12, height = 8)
par(mfrow = grid_type, mar = c(3,3,2,1))

    for (i in 1:nrow(all_species_matrix)){

        histograms <- lapply(1:3, function(x){ hist(all_species_matrix[i, species_list[[x]]], breaks = ax, plot = F) })

    # plot stacked histogram
        histograms[[1]]$counts <- histograms[[1]]$counts + histograms[[2]]$counts + histograms[[3]]$counts
        histograms[[2]]$counts <- histograms[[2]]$counts + histograms[[3]]$counts

        for (clad in 1:3){
            add <- ifelse(clad == 1, F, T)
            plot(histograms[[clad]], col = taxonomy_colors[clad], add = add, main = '', xlab = '', ylab = '', border = F)
        }

        # 
        legend('topleft', paste('ALG', i), bty = 'n', cex = cex, text.font = 2, fill = pal[ALGs[i]])
        if (i == 3){
            legend('topright', legend = c('Nematocera', 'Brachycera', 'Schizophora'), fill = taxonomy_colors, bty = 'n', cex = cex) # , lwd = 2, lty = 1
        }
    }
dev.off()

pdf(paste0(args$o, '_nematocera.pdf'))
# pdf('figures/ALG_stability_histograms_Nematocera.pdf')
    par(mfrow = grid_type, mar = c(3,3,2,1))

    for (i in 1:nrow(all_species_matrix)){
        hist(all_species_matrix[i, nematocera_species], breaks = 20, main = NA, col = pal[ALGs[i]], xlim = c(0,1), border = NA, ylab = NA, xlab = NA, cex.axis = cex)
        legend('topleft', ALGs[i], bty = 'n', cex = cex, text.font = 2)
    }

dev.off()

pdf(paste0(args$o, '_brachycera.pdf'), width = 12, height = 8)
# pdf('figures/ALG_stability_histograms_Nematocera.pdf')
    par(mfrow = grid_type, mar = c(3,3,2,1))

    for (i in 1:nrow(all_species_matrix)){
        hist(all_species_matrix[i, brachycera_species], breaks = 20, main = NA, col = pal[ALGs[i]], xlim = c(0,1), border = NA, ylab = NA, xlab = NA, cex.axis = cex)
        legend('topleft', ALGs[i], bty = 'n', cex = cex, text.font = 2)
    }

dev.off()


pdf(paste0(args$o, '_schizophora.pdf'), width = 12, height = 8)
# pdf('figures/ALG_stability_histograms_Nematocera.pdf')
    par(mfrow = grid_type, mar = c(3,3,2,1))

    for (i in 1:nrow(all_species_matrix)){
        hist(all_species_matrix[i, schizophora_species], breaks = 10, main = NA, col = pal[ALGs[i]], xlim = c(0,1), border = NA, ylab = NA, xlab = NA, cex.axis = cex)
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

