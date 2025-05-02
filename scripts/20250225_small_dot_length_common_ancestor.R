# infer the size of the small dot equivalent in the common ancestor
################################################################################
library(dplyr)
library(stringr)
library(ggplot2)
library(tibble)
library(ape)
library(ggtree)
################################################################################
# from 20250222_small_dot_length_for_tree.R
################################################################################
# small dot status
file_path <- "C:/Users/julia/Documents/Kamil/txt/small_dot_status.csv"
small_dot_defin <- read.csv(file_path, header = FALSE, sep = ",")
# make a copy to change
small_dot_defin_cp <- small_dot_defin
colnames(small_dot_defin_cp) <- c("label", "status","position")
head(small_dot_defin_cp)

# lengths of small dots
# manually removed 2 entries
# GCA_963971445.1.tsv	OZ020531.1	76	1	43148219 # Ornithomya_chloropus, not small dot
# GCA_949715645.1.tsv	OX454530.1	2	1	8431506    # Portevinia_maculate, two small dots, check with curation
file_path <- "C:/Users/julia/Documents/Kamil/txt/filtered_potential_small_dots_with_length.txt"
small_d_table <- read.delim(file_path, header = TRUE, sep = "\t")
head(small_d_table)

# more info for merging
names_acc_fam_2 <- read.table("C:/Users/julia/Documents/Kamil/busco_fulltable/names_acc_fam_2.csv", sep = ",", header = TRUE)
head(names_acc_fam_2)

# preparing table with lengths to merge
small_d_table <- small_d_table %>%
  mutate(Accession = str_remove(Genome, "\\..*\\.tsv"))
# add column with label (names)
small_d_table <- small_d_table %>%
  left_join(names_acc_fam_2 %>% select(Accession, label), by = "Accession")
head(small_d_table)

# merge length to status info
small_dot_defin_cp <- small_dot_defin_cp %>%
  left_join(small_d_table %>% select(label, Chromosome_Length), by = "label")
head(small_dot_defin_cp)
na_count <- sum(is.na(small_dot_defin_cp$Chromosome_Length)) # 76, therefore warning message
# manuelly adjust to NA, something is wrong with this assembly
# Tephritidae that is absolutely mental (GCA_027475135.1)
small_dot_defin_cp <- small_dot_defin_cp %>%
  mutate(Chromosome_Length = if_else(label == "Bactrocera_correcta", NA_real_, Chromosome_Length))
head(small_dot_defin_cp)
# finished preparing correct data frame
################################################################################
# preparing tree
# tree information for minus_26 (5 ALGs at n3, -m 89 -r 2 -a quick)
tree2 <- read.tree("C:/Users/julia/Documents/Kamil/trees/syngraph.infer.m89.2.quick.newick_minus26.txt")
names <- read.table("C:/Users/julia/Documents/Kamil/trees/acc_names_match.txt", sep = ",")
cluster <- read.table("C:/Users/julia/Documents/Kamil/trees/syngraph.infer.m89.2.quick.clusters_minus26.tsv")
################################################################################
# exchange accession number for species names on the tree
tip_labels <- tree2$tip.label
name_mapping <- setNames(names$V2, names$V1)
new_tip_labels <- ifelse(tip_labels %in% names(name_mapping), 
                         name_mapping[tip_labels], 
                         tip_labels)
tree2$tip.label <- new_tip_labels
################################################################################
# remove tips where the small dot is not present
length(tree2$tip.label)
tips_to_drop <- small_dot_defin$V1[small_dot_defin$V2 != "p" | is.na(small_dot_defin$V2)]
tree2_pruned <- drop.tip(tree2, tips_to_drop)
length(tree2_pruned$tip.label) # 212 out of 289 left, dropped 76

# prepare chromosome_length
chromosome_lengths <- setNames(small_dot_defin_cp$Chromosome_Length, small_dot_defin_cp$label)
chromosome_lengths <- chromosome_lengths[names(chromosome_lengths) %in% tree2_pruned$tip.label]

# reconstruct length of chromosome at internal nodes
estimate_internal_lengths <- function(tree, tip_lengths) {
  internal_lengths <- ace(tip_lengths, tree, method = "pic")$ace
  return(internal_lengths)
}
# apply function to tree
internal_lengths <- estimate_internal_lengths(tree2_pruned, chromosome_lengths)
internal_lengths

# plot tree with lengths at internal nodes
png("ALG_6_length_common_ancestor.png", width = 2100, height = 3000, res = 800)
plot(tree2_pruned, show.tip.label = TRUE, cex = 0.09, direction = "rightwards")
nodelabels(round(internal_lengths, 2), cex = 0.1, frame = "none", col = "blue", adj = c(0.8, -1.2))
dev.off()

