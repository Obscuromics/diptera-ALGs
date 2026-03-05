library(ggtree)
library(gsheet)
library(ggplot2)
library(dplyr)
library(tidyr)
library(phytools)

# ./data/syngraph/ARCHIVE/diptera.no_chiro.syngraph_infer.min_dist.m165.newick.txt
tree <- read.tree('data/syngraph/diptera.no_plecia.mindist.m165.newick.txt')
# plot(tree)

all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?usp=sharing", format='csv'),
                            stringsAsFactors = F, header = T, check.names = F)
genome_data <- all_genome_data[all_genome_data[, 'TO ADD'] == 'KEEP', ]
species_family_tab <- genome_data[, c('species', 'family')]

one_species_per_family_df <- species_family_tab %>%
  distinct(family, .keep_all = TRUE)

length(tree$tip.label)

tips_to_keep <- c(one_species_per_family_df$species)
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, tips_to_keep))

plot(pruned_tree)
pruned_tree$tip.label

rownames(species_family_tab) <- species_family_tab[, 1]
pruned_tree$tip.label <- species_family_tab[pruned_tree$tip.label, 2]

# sorted_one_species_per_family_df <- one_species_per_family_df %>%
#   mutate(species = factor(species, levels = pruned_tree$tip.label)) %>%
#   arrange(species)

# tip_labels <- c(sorted_one_species_per_family_df$family)
# tip_labels

# plot(pruned_tree,  show.node.label = TRUE)

tree_plt <- ggtree(pruned_tree, layout = "roundrect", size = 1.5) + 
  geom_tiplab(size = 10, align = TRUE, offset = 0.05, linesize = 1, linetype = "dotted") +
  geom_nodelab(size = 10) +
  xlim(NA, 1.5) +
  theme(plot.margin = margin(10, 10, 10, 10),
        legend.position = "right")

ggsave("figures/tree_families_syngraph_nodelabels.png", plot = tree_plt, width = 30, height = 30, dpi = 200)
ggsave("figures/tree_families_syngraph_nodelabels.svg", plot = tree_plt, width = 30, height = 30, dpi = 600)