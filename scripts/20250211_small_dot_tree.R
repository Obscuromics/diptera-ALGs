# plot a tree with 289 species
# where color of branches indicates fate of ALG 6
# possible fates: present (p), fused (f), dispersed (s), on unplaced scaffolds (a)
################################################################################
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
################################################################################
# tree information for minus_26 (5 ALGs at n3, -m 89 -r 2 -a quick)
tree2 <- read.tree("C:/Users/julia/Documents/Kamil/trees/syngraph.infer.m89.2.quick.newick_minus26.txt")
names <- read.table("C:/Users/julia/Documents/Kamil/trees/acc_names_match.txt", sep = ",")
################################################################################
# exchange accession number for species names on the tree
tip_labels <- tree2$tip.label
name_mapping <- setNames(names$V2, names$V1)
new_tip_labels <- ifelse(tip_labels %in% names(name_mapping), 
                         name_mapping[tip_labels], 
                         tip_labels)
tree2$tip.label <- new_tip_labels
################################################################################
# make a ggtree
tt <- ggtree(tree2)
################################################################################
# visualize small dot presence

# I have manually/ qualitatively assessed fate of ALG 6 in each genome
# by looking at paintings of genomes in colours of ALGs
file_path <- "C:/Users/julia/Documents/Kamil/txt/small_dot_status.csv"
small_dot_defin <- read.csv(file_path, header = FALSE, sep = ",")
head(small_dot_defin)

# fusion of information on small dot to tree data object
tt$data <- tt$data %>%
  left_join(small_dot_defin, by = c("label" = "V1")) %>%
  rename(small_dot = V2)
tt$data$small_dot[is.na(tt$data$small_dot)] <- "i"
################################################################################
s <- ggtree(tree2) +
  geom_tree(data = tt$data, aes(color = small_dot), size = 0.5) +
  geom_tiplab(data = tt$data, size = 0.5, align = TRUE, offset = 0.05, 
              linetype = "dotted", linesize = 0.05) +
  xlim(NA, 2) +
  scale_color_manual(
    #values = c("s" = "blue", "f" = "darkgreen", "a" = "red", "p" = "darkgrey"),
    values = c("s" = "blue", "f" = "darkgreen", "a" = "red", "p" = "black"),
    name = "Fate of the small dot chromosome",
    labels = c("s" = "scattered across chromosomes",
               "f" = "fused to a specific chromosome",
               "a" = "absent or on unplaced scaffolds",
               "p" = "present",
               "i" = "internal branch")
  ) +
  theme(legend.position = c(0.8, 0.7),  # Position legend inside the plot
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))

s

ggsave("tree_small_d_black.png", plot = s, width = 10, height = 8, dpi = 600)
################################################################################
# check how many families are included in genomes labeled as small dot "p" present
names_acc_fam_2 <- read.csv("C:\\Users\\julia\\Documents\\Kamil\\busco_fulltable\\names_acc_fam_2.csv", 
               sep = ",", header = TRUE, stringsAsFactors = FALSE)
head(names_acc_fam_2)
colnames(names_acc_fam_2) <- c("names", "num", "acc", "fam")

df <- data.frame(label = tt$data$label, small_dot = tt$data$small_dot) %>%
  left_join(names_acc_fam_2, by = c("label" = "names"))

# check number of families (fam) or species (names)
# that are present (p), fused (f) etc
num_families <- df %>%
  filter(small_dot == "p") %>%
  distinct(fam) %>%
  nrow()
print(num_families) # 37
################################################################################


