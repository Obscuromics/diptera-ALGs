# Figure 2A - phyligenetic tree with families

library(ape)
library(ggtree)
library(gsheet)
library(ggplot2)
library(dplyr)
library(tidyr)
library(phytools)

tree <- read.tree('data/diptera.supermatrix.phy.treefile')
all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?usp=sharing", format='csv'),
                            stringsAsFactors = F, header = T, check.names = F)
genome_data <- all_genome_data[all_genome_data[, 'TO ADD'] == 'KEEP', ]
species_family_tab <- genome_data[, c('species', 'family')]

one_species_per_family_df <- species_family_tab %>%
  distinct(family, .keep_all = TRUE)

tips_to_keep <- c('Panorpa_germanica', one_species_per_family_df$species)
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, tips_to_keep))

sorted_one_species_per_family_df <- one_species_per_family_df %>%
  mutate(species = factor(species, levels = pruned_tree$tip.label)) %>%
  arrange(species)

tip_labels <- c(sorted_one_species_per_family_df$family, 'Panorpidae')
pruned_tree$tip.label <- tip_labels

rooted_tree <- reroot(pruned_tree, node.number = which(tip_labels == 'Panorpidae'), position =  node.depth.edgelength(pruned_tree)[which(tip_labels == 'Panorpidae')]/4)

tree_plt <- ggtree(rooted_tree, layout = "roundrect", size = 1.5) + 
  geom_tiplab(size = 10, align = TRUE, offset = 0.05, linesize = 1, linetype = "dotted") +
  geom_rootedge() +
  xlim(NA, 1.5) +
  theme(plot.margin = margin(10, 10, 10, 10),
        legend.position = "right")

ggsave("figures/tree_families_only_wide.png", plot = tree_plt, width = 30, height = 30, dpi = 200)
ggsave("figures/tree_families_only_wide.svg", plot = tree_plt, width = 30, height = 30, dpi = 600)
