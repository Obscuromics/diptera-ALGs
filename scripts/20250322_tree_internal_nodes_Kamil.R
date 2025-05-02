# plot tree with internal nodes and species coloured by number of
# chromosomes or ancestral linkage groups
################################################################################
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
################################################################################
# minus_26, final tree
tree2 <- read.tree("C:/Users/julia/Documents/Kamil/trees/syngraph.infer.m89.2.quick.newick_minus26.txt")
names <- read.table("C:/Users/julia/Documents/Kamil/trees/acc_names_match.txt", sep = ",")
cluster <- read.table("C:/Users/julia/Documents/Kamil/trees/syngraph.infer.m89.2.quick.clusters_minus26.tsv")

# version before final tree, 9 diptera and 2 outgroup species removed
# from all_11 (minus 9 most problematic genomes)
#tree2 <- read.tree("C:/Users/julia/Documents/Kamil/trees/syngraph.infer.m89.2.quick.newick.txt")
#names <- read.table("C:/Users/julia/Documents/Kamil/trees/acc_names_match.txt", sep = ",")
#cluster <- read.table("C:/Users/julia/Documents/Kamil/trees/syngraph.infer.m89.2.quick.clusters.tsv")

# version with all species I have ever used
#tree2 <- read.tree("C:/Users/julia/Documents/Kamil/trees/syngraph.infer.m80.2.quick.newick.txt")
#names <- read.table("C:/Users/julia/Documents/Kamil/trees/acc_names_match.txt", sep = ",")
#cluster <- read.table("C:/Users/julia/Documents/Kamil/trees/syngraph.infer.m80.2.quick.clusters.tsv")
################################################################################
# exchange accession number for species names at tips
tip_labels <- tree2$tip.label
name_mapping <- setNames(names$V2, names$V1)
new_tip_labels <- ifelse(tip_labels %in% names(name_mapping), 
                         name_mapping[tip_labels], 
                         tip_labels)
tree2$tip.label <- new_tip_labels
################################################################################
# exchange accession number for species names in the cluster table
name_mapping <- setNames(names$V2, names$V1)
new_cluster_V2 <- ifelse(cluster$V2 %in% names(name_mapping), 
                         name_mapping[cluster$V2], 
                         cluster$V2)

cluster$V2 <- new_cluster_V2
cluster <- cluster[,-1]
colnames(cluster) <- c("node", "chromosome_number")
################################################################################
# add missing values (outgrous and root node) to the cluster table
# because outgroups are missing and common node of in and outgroups too

# this if for the final version
neue_werte_u <- data.frame(Spalte1 = c("Panorpa_germanica"),
                           Spalte2 = c(21))
neue_werte_o <- data.frame(Spalte1 = c(""),
                           Spalte2 = c(0))

# in the earliest version there are two more outgroups that have to be adjusted
#neue_werte_u <- data.frame(Spalte1 = c("Limnephilus_lunatus", "Panorpa_germanica", "Danaus_plexippus"),
                           #Spalte2 = c(13,21,32))
#neue_werte_o <- data.frame(Spalte1 = c(""),
                           #Spalte2 = c(0))

colnames(neue_werte_u) <- c("node", "chromosome_number")
colnames(neue_werte_o) <- c("node", "chromosome_number")
cluster2 <- rbind(cluster, neue_werte_u)
cluster3 <- rbind(neue_werte_o, cluster2)

# extract tips and nodes seperatly from the cluster table
filtered_cluster3 <- cluster3[nchar(as.character(cluster3[[1]])) < 5, ]
names_cluster3 <- cluster3[nchar(as.character(cluster3[[1]])) > 4, ]
################################################################################
# make a ggtree
tt <- ggtree(tree2, layout = "roundrect")

# generate df from both parts of cluster table in the right order to attach it to the trees data attribute
# check if this is really the right order!
#keep both columns for separate matching
combined_column <- rbind(names_cluster3, filtered_cluster3)
new_df <- data.frame(combined_column)
names(new_df)[1] <- "label"
tt$data <- left_join(tt$data, new_df, by = "label")

# attach the column of the data frame to the data object from the ggtree
# replace the value 0 at the root with 3 so that it does not show in the legend
tt$data <- tt$data %>%
  mutate(chromosome_number = ifelse(chromosome_number == 0, 3, chromosome_number))
################################################################################
p <- tt + 
  geom_tippoint(aes(color = factor(chromosome_number)), size = 0.5) +
  geom_nodepoint(aes(color = factor(chromosome_number)), size = 0.5) +
  geom_tiplab(aes(color = factor(chromosome_number)),
              size = 0.5, align = TRUE, offset = 0.05,
              linetype = "dotted", linesize = 0.05,
              show.legend = FALSE ) + xlim(NA, 2) +
  labs(color = "Chromosome Number") +
  #theme(legend.text = element_text(size = 5),
  #legend.title = element_text(size = 5),
  #legend.position = c(0.13, 0.6)) +
  guides(color = guide_legend(override.aes = list(size = 3)))
p

# save the tree as an image
ggsave("tree_trial_final.png", plot=p, width=10, height=8, dpi=600)
################################################################################
# plot tree with names at internal nodes
rr <- tt + geom_tiplab(size = 0.5, align = TRUE)
tail(rr$data)
tt <- tt + geom_text2(aes(subset = !isTip, label = label, color = "red"), size = 1)
tt
ggsave("tree_m71_q2_node_names.png", plot=rr, width=10, height=8, dpi=600)
################################################################################