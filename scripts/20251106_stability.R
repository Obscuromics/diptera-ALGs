suppressPackageStartupMessages(library('gsheet'))
suppressPackageStartupMessages(library('argparse'))
suppressPackageStartupMessages(library('ape'))
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library("patchwork"))

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

busco_asn_file <- args$a

ALG_data <- read.table(busco_asn_file, col.names = c('busco', 'ALG'))
ALG_sizes <- table(ALG_data[, 'ALG'])
row.names(ALG_data) <- ALG_data[, 'busco']

process_species_data <- function(species_name){
    busco_on_chr <- read.table(paste0('data/busco_tables/', species_name, '.syngraph.buscos.tsv'), col.names = c('busco', 'chr', 'f', 't'))[, c(1,2)]
    busco_on_chr[, 'alg'] <- ALG_data[busco_on_chr[, 'busco'], 'ALG']

    chr_busco_matrix <- table(busco_on_chr[, 'alg'], busco_on_chr[, 'chr'])
    apply(chr_busco_matrix, 1, max) / ALG_sizes
}

all_species_matrix <- sapply(species2suborder[, 1], process_species_data)
ALGs <- row.names(all_species_matrix)

if (length(ALGs) == 6){
    grid_type <- c(2, 3)
} else {
    grid_type <- c(3, 4)
}
cex = 1.4
taxonomy_colors <- c("Nematocera" = "#67d7dbff", "Brachycera" = "#069c24ff", "Schizophora" = "#bb3fa0ff")

nematocera_species <- all_genome_data[all_genome_data[, 'suborder'] == "Nematocera", 'species']
tree <- read.tree("data/syngraph/diptera.no_plecia.mindist.m165.newick.txt")
schizophora_species <- extract.clade(tree, node = "n133")$tip.label

brachycera_species <- all_genome_data[all_genome_data[, 'suborder'] == "Brachycera", 'species']
brachycera_species <- brachycera_species[!brachycera_species %in% schizophora_species]

all_species_matrix <- all_species_matrix[, 1:340] #removing Panorpa

nematocera_species_matrix <- all_species_matrix[, nematocera_species]
species_list <- list(nematocera_species, brachycera_species, schizophora_species)

stability_heatplot <- function(species_list, colour){
    n_bins <- 10
    breaks <- seq(0, 1, length.out = n_bins + 1)

    heat_df <- map_dfr(1:nrow(all_species_matrix), function(i) {
      x <- all_species_matrix[i, species_list]
      x <- x[is.finite(x) & x >= 0 & x <= 1]
      h <- hist(
        x,
        breaks = breaks,
        plot = FALSE
      )
      
      tibble(
        ALG = ALGs[i],
        bin = (breaks[-1] + breaks[-length(breaks)]) / 2,
        density = h$density / max(h$density)
      )
    })

    ggplot(heat_df, aes(x = bin, y = ALG, fill = density)) +
      geom_tile() +
      scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(9, colour),
      limits = c(0, 1),
      name = "Density"
      )+
      scale_x_continuous(limits = c(0, 1)) +
      labs(x = "Stability", y = NULL) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank()
      )
}

compact_stability_heatplot <- function(species_list, colour){
  n_bins <- 10
  breaks <- seq(0, 1, length.out = n_bins + 1)
  bin_centers <- (breaks[-1] + breaks[-length(breaks)]) / 2

  heat_df <- map_dfr(seq_len(nrow(all_species_matrix)), function(i) {
    x <- all_species_matrix[i, species_list]
    x <- x[is.finite(x) & x >= 0 & x <= 1]

    if (length(x) == 0) {
      return(tibble(
        ALG = ALGs[i],
        bin = bin_centers,
        density = 0
      ))
    }

    h <- hist(x, breaks = breaks, plot = FALSE)

    tibble(
      ALG = ALGs[i],
      bin = bin_centers,
      density = h$density / max(h$density)
    )
  })

  heat_df <- heat_df %>%
  mutate(
    ALG_num = as.numeric(factor(ALG, levels = unique(ALGs)))
  )


  ggplot(heat_df, aes(x = bin, y = ALG, fill = density)) +
    geom_tile(height = 1) +   # ↓ tighter rows
    scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(9, colour),
      limits = c(0, 1),
      name = "Density"
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Stability", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank(),
      panel.spacing = unit(0, "lines")
    )
}


p1 <- compact_stability_heatplot(nematocera_species, "Blues")+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())+
        ggtitle("Nematocera")
p2 <- compact_stability_heatplot(brachycera_species, "Greens")+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())+
        ggtitle("Brachycera")
p3 <- compact_stability_heatplot(schizophora_species, "Purples")+
      ggtitle("Schizophora")

# make a stack with each of the 3 groups coloured by the fig 2 colouring
# repeat for alpha, beta, and gamma algs
# refine

plt_all <- (p1 | p2 | p3) +
  plot_layout(nrow = 3, heights = c(0.8, 0.8, 0.8)) &
  theme(legend.position = "none")

ggsave("compact_alpha_heatmap.pdf", plot=plt_all, dpi=600, width = 7, height = 15)

