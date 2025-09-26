# Figure 2A - phyligenetic tree with families
################################################################################
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
################################################################################
#root <- "/Users/ab66/Documents/sanger_work/diptera/diptera-ALGs/"
root <- paste0(getwd(), "/")
################################################################################
# A - tree with families
diptera_info <- read.table(paste0(root, "data/names_acc_fam_2.csv"), sep = ",", header = TRUE)
colnames(diptera_info) <- c("label", "num", "acc", "fam")

# read tree file
tree <- read.tree(paste0(root, "data/315_supermatrix.phy.treefile"))
names <- read.table(paste0(root, "data/acc_names_match.txt"), sep = ",")

# tips, die da bleiben sollen
all_fam_sex_vector <- readLines(paste0(root, "data/keep_one_per_family.csv"))
mecop <- c("GCA_963678705") # mecoptera outgroup

# prune tree
keep_tips <- c(mecop, all_fam_sex_vector)
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, keep_tips))

# change name of tips to families
tip_labels_df <- data.frame(label = pruned_tree$tip.label)
tip_labels_df <- left_join(tip_labels_df, diptera_info[, c("acc", "fam")], by = c("label" = "acc"))
pruned_tree$tip.label <- tip_labels_df$fam

tree_plt <- ggtree(pruned_tree, layout = "roundrect", size = 0.8) + 
  #geom_tree(size = 0.8) +
  geom_tiplab(size = 10, align = TRUE, offset = 0.05, linesize = 1, linetype = "solid") +
  xlim(NA, 1.5) +
  theme(plot.margin = margin(10, 10, 10, 10),
        legend.position = "right")
#tree_plt

ggsave(paste0(root, "figures/tree_families_only_wide.png"), plot = tree_plt, width = 30, height = 30, dpi = 600)
ggsave(paste0(root, "figures/tree_families_only_wide.svg"), plot = tree_plt, width = 30, height = 30, dpi = 600)

################################################################################
#sam edit
#this is taking the first species in each family as a family representative

library(ape)
library(ggtree)
library(gsheet)
library(ggplot2)
library(dplyr)
library(tidyr)
library(phytools)

tree <- read.tree('diptera.supermatrix.phy.treefile')
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
