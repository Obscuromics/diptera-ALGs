
################################################################################
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
################################################################################
tree_path <- "C:/Users/julia/Documents/Kamil/trees/315_supermatrix.phy.treefile"
tree <- read.tree(tree_path)
names <- read.table("C:/Users/julia/Documents/Kamil/trees/acc_names_match.txt", sep = ",")

all_fam_sex <- read.table("C:/Users/julia/Documents/Kamil/busco_fulltable/keep_one_per_family.csv", sep = ",", header = FALSE)
mecop <- c("GCA_963678705")
################################################################################
names_acc_fam_2 <- read.table("C:/Users/julia/Documents/Kamil/busco_fulltable/names_acc_fam_2.csv", sep = ",", header = TRUE)
colnames(names_acc_fam_2) <- c("names", "num", "acc", "fam")
head(names_acc_fam_2)
diptera_info <- names_acc_fam_2
colnames(diptera_info) <- c("label", "num", "acc", "fam")
head(diptera_info)

all_fam_sex <- merge(all_fam_sex, names_acc_fam_2[, c("acc", "fam")], 
                     by.x = "V1", by.y = "acc", all.x = TRUE)
head(all_fam_sex)
keep_fams <- c("Drosophilidae",
               "Empididae",
               "Syrphidae",
               "Ptychopteridae",
               "Bibionidae",
               "Culicidae",
               "Chironomidae",
               "Hippoboscidae",
               "Anthomyiidae",
               "Stratiomyidae",
               "Tephritidae",
               "Cecidomyiidae",
               "Muscidae",
               "Diopsidae")
all_fam_sex <- all_fam_sex[all_fam_sex$fam %in% keep_fams, ]

all_fam_sex_vector <- all_fam_sex$V1
keep_tips <- c(mecop, all_fam_sex_vector)
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, keep_tips)) # works up to here

tip_labels_df <- data.frame(label = pruned_tree$tip.label)
tip_labels_df <- left_join(tip_labels_df, diptera_info[, c("acc", "fam")], by = c("label" = "acc"))
pruned_tree$tip.label <- tip_labels_df$fam
pruned_tree$tip.label[pruned_tree$tip.label == "Hippoboscidae"] <- "Glossinidae"
pruned_tree$tip.label[pruned_tree$tip.label == "Panorpidae"] <- "Mecoptera"

# cladogram
pruned_tree <- ladderize(pruned_tree, right = FALSE)

plot_cladogram <- function(tree) {
  pastel_green <- "palegreen3"
  dark_yellow <- "darkgoldenrod2"
  pastel_green <- "black"
  dark_yellow <- "black"
  gray_color <- "gray50"
  
  pastel_green_tips <- c("Ptychopteridae", "Bibionidae", "Culicidae", "Chironomidae", "Cecidomyiidae")
  panorpidae_tip <- "Mecoptera"
  
  tip_colors <- rep("black", length(tree$tip.label))
  label_colors <- rep("black", length(tree$tip.label))
  
  for (tip in pastel_green_tips) {
    tip_index <- which(tree$tip.label == tip)
    label_colors[tip_index] <- pastel_green
  }
  
  for (tip in tree$tip.label) {
    if (!(tip %in% pastel_green_tips) && tip != panorpidae_tip) {
      tip_index <- which(tree$tip.label == tip)
      label_colors[tip_index] <- dark_yellow
    }
  }
  
  # Setze das Label von Panorpidae auf grau
  label_colors[tree$tip.label == panorpidae_tip] <- gray_color
  
  edge_colors <- rep("black", length(tree$edge[, 1]))
  
  for (tip in pastel_green_tips) {
    tip_index <- which(tree$tip.label == tip)
    edge_colors[tree$edge[, 2] == tip_index] <- "black"
  }
  
  for (tip in tree$tip.label) {
    if (!(tip %in% pastel_green_tips) && tip != panorpidae_tip) {
      tip_index <- which(tree$tip.label == tip)
      edge_colors[tree$edge[, 2] == tip_index] <- "black"
    }
  }
  
  panorpidae_index <- which(tree$tip.label == panorpidae_tip)
  edge_colors[tree$edge[, 2] == panorpidae_index] <- gray_color
  
  plot(tree, type = "cladogram", direction = "rightward", 
       use.edge.length = FALSE, no.margin = TRUE, 
       cex = 2.1, edge.width = 2, tip.color = label_colors, edge.color = edge_colors)
}

png("cladogram.png", width = 10, height = 15, units = "in", res = 600)
plot_cladogram(pruned_tree)
dev.off()


