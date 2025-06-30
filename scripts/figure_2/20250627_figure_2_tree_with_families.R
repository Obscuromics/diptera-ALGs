# Figure 2A - phyligenetic tree with families
################################################################################
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
################################################################################
root <- "/Users/ab66/Documents/sanger_work/diptera/diptera-ALGs/"
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