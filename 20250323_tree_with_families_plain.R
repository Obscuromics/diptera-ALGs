# plot a simple tree with 52 families plus mecoptera
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
################################################################################
names_acc_fam_2 <- read.table("C:/Users/julia/Documents/Kamil/busco_fulltable/names_acc_fam_2.csv", sep = ",", header = TRUE)
colnames(names_acc_fam_2) <- c("names", "num", "acc", "fam")
head(names_acc_fam_2)
diptera_info <- names_acc_fam_2
colnames(diptera_info) <- c("label", "num", "acc", "fam")
head(diptera_info)

tree_path <- "C:/Users/julia/Documents/Kamil/trees/315_supermatrix.phy.treefile"
tree <- read.tree(tree_path)
names <- read.table("C:/Users/julia/Documents/Kamil/trees/acc_names_match.txt", sep = ",")
# tips, die da bleiben sollen
all_fam_sex <- read.table("C:/Users/julia/Documents/Kamil/busco_fulltable/keep_one_per_family.csv", sep = ",", header = FALSE)
mecop <- c("GCA_963678705") # mecoptera outgroup
# prune tree
all_fam_sex_vector <- all_fam_sex$V1
keep_tips <- c(mecop, all_fam_sex_vector)
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, keep_tips))

# change name of tips to families
tip_labels_df <- data.frame(label = pruned_tree$tip.label)
tip_labels_df <- left_join(tip_labels_df, diptera_info[, c("acc", "fam")], by = c("label" = "acc"))
#tip_labels_df <- left_join(tip_labels_df, diptera_info[, c("label", "fam")], by = "label")
tip_labels_df
pruned_tree$tip.label <- tip_labels_df$fam

# plot tree
t <- ggtree(pruned_tree, layout = "roundrect", size = 0.8) + 
  #geom_tree(size = 0.8) +
  geom_tiplab(size = 10, align = TRUE, offset = 0.05, linesize = 1, linetype = "solid") +
  xlim(NA, 1.5) +
  theme(plot.margin = margin(10, 10, 10, 10),
        legend.position = "right")
t
ggsave("tree_families_only_wide.png", plot = t, width = 30, height = 30, dpi = 600)

