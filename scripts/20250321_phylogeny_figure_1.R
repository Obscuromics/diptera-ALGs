# phylogeny with 312 diptera and 1 mecoptera species
# excluded species marked in red
# separately plot traits: BUSCO completeness, genome size, haploid chromosome number
################################################################################
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(phytools)
################################################################################
tree <- read.tree("C:/Users/julia/Documents/Kamil/trees/315_supermatrix.phy.treefile")
names <- read.table("C:/Users/julia/Documents/Kamil/trees/acc_names_match.txt", sep = ",")

# change accessions to species names
tip_labels <- tree$tip.label
name_mapping <- setNames(names$V2, names$V1)
new_tip_labels <- ifelse(tip_labels %in% names(name_mapping), 
                         name_mapping[tip_labels], 
                         tip_labels)
tree$tip.label <- new_tip_labels

# remove two unused outgroups
tip_to_remove <- c("Danaus_plexippus",
                   "Limnephilus_lunatus")
tree_clean <- drop.tip(tree, tip_to_remove)

# color 24 species that were excluded
colour_tips <- c("Myopa_testacea",
                 "Melieria_omissa",
                 "Drosophila_pseudoobscura",
                 "Eupeodes_latifasciatus",
                 "Syrphus_vitripennis",
                 "Drosophila_helvetica",
                 "Drosophila_lowei",
                 "Baccha_elongata",
                 "Pseudolycoriella_hygida",
                 "Drosophila_miranda",
                 "Drosophila_athabasca",
                 "Drosophila_affinis",
                 "Anopheles_parensis",
                 "Anopheles_longipalpis",
                 "Anopheles_rivulorum",
                 "Drosophila_nebulosa",
                 "Drosophila_gunungcola",
                 "Drosophila_willistoni",
                 "Anopheles_vaneedeni",
                 "Scaeva_pyrastri",
                 "Eupeodes_luniger",
                 "Meiosimyza_platycephala",
                 "Melieria_crassipennis",
                 "Eupeodes_corollae")

# prepare scale (assume that flies diverged from last common ancestor 260 Mio
# years ago)
#max_depth <- max(nodeHeights(tree_clean))
#scaling_factor <- 260 / max_depth
#tree_clean$edge.length <- tree_clean$edge.length * scaling_factor

# plot tree
#tree_show <- ggtree(tree_clean) +
  #geom_tree(aes(color = label %in% colour_tips), size = 0.5) +
  #geom_tiplab(aes(color = "black"), align = TRUE, offset = 0.05, 
              #linetype = "dotted", linesize = 0.05, size = 1.2) +
  #scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  #xlim(NA, 500) +
  #theme(legend.position = "none") +
  #geom_treescale(x = 80, y = -10, width = 50)


tree_show <- ggtree(tree_clean) +
  geom_tree(aes(color = label %in% colour_tips), size = 0.5) +
  geom_tiplab(aes(color = "black"), align = TRUE, offset = 0.05, 
              linetype = "dotted", linesize = 0.05, size = 1.2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  theme(legend.position = "none") 

tree_show


# safe image
ggsave("tree_figure_1.png", plot=tree_show, width=10, height=13, dpi=600)
################################################################################
# plot traits
# extract correct order of species
data <- tree_show$data

sorted_data_tips_desc <- data %>%
  filter(isTip == TRUE) %>%
  arrange(desc(y))

#write.table(sorted_data_tips_desc$label, file = "labels_sorted.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
################################################################################
# plot the three traits
details <- read.table("C:/Users/julia/Documents/Kamil/txt/20250318_master_sheet_Diptera - data_figure_1.tsv", sep = "\t", header = FALSE)
details$group <- (cumsum(details$V2 != dplyr::lag(details$V2, default = details$V2[1])) %% 2) + 1
head(details)

# plot haploid chromosome number
fig1_hap <- ggplot(details, aes(x = V5, y = factor(V4, levels = rev(unique(V4))), fill = factor(group))) +
  geom_bar(stat = "identity") +
  xlab("V5") +
  ylab("V4 (Reversed)") +
  scale_fill_manual(values = c("1" = "gray30", "2" = "gray60")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        #panel.grid.major.x = element_line(color = "gray", size = 0.5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

fig1_hap
ggsave("figure_1_haploidnumber.png", plot=fig1_hap, width=2, height=13, dpi=600)

# plot BUSCO completeness
fig1_BUS <- ggplot(details, aes(x = V7, y = factor(V4, levels = rev(unique(V4))), fill = factor(group))) +
  geom_bar(stat = "identity") +
  xlab("V7") +
  ylab("V4 (Reversed)") +
  scale_fill_manual(values = c("1" = "gray30", "2" = "gray60")) +
  scale_x_continuous(breaks = c(0, 50, 100)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        # panel.grid.major.x = element_line(color = "gray", size = 0.5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")
fig1_BUS
ggsave("figure_1_BUSCO_complete.png", plot=fig1_BUS, width=1, height=13, dpi=600)

# plot genome size
fig1_size <- ggplot(details, aes(x = V6, y = factor(V4, levels = rev(unique(V4))), fill = factor(group))) +
  geom_bar(stat = "identity") +
  xlab("V6 (Millionen)") +
  ylab("V4 (Reversed)") +
  scale_fill_manual(values = c("1" = "gray30", "2" = "gray60")) +
  scale_x_continuous(
    breaks = seq(0, max(details$V6, na.rm = TRUE), by = 500e6),
    labels = scales::label_number(scale = 1e-6)
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        # panel.grid.major.x = element_line(color = "gray", size = 0.5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")
fig1_size
ggsave("figure_1_size.png", plot=fig1_size, width=2, height=13, dpi=600)
################################################################################
# plot BUSCO completeness for removed families
# for supplementary figures
removed_genomes <- read.table("C:/Users/julia/Documents/Kamil/txt/removed.txt", 
                              header = FALSE, sep = "\t", stringsAsFactors = FALSE)
head(removed_genomes)
details <- read.table("C:/Users/julia/Documents/Kamil/txt/20250318_master_sheet_Diptera - data_figure_1.tsv", sep = "\t", header = FALSE)
details$group <- (cumsum(details$V2 != dplyr::lag(details$V2, default = details$V2[1])) %% 2) + 1
head(details)

details$V9_clean <- sub("\\.\\d+$", "", details$V9)
filtered_details <- details[details$V9_clean %in% removed_genomes$V1, ]
filtered_details <- filtered_details[match(removed_genomes$V1, filtered_details$V9_clean), ]
head(filtered_details)

fig1_BUS_rem <- ggplot(filtered_details, aes(x = V7, y = factor(V4, levels = rev(unique(V4))))) +
  geom_bar(stat = "identity") +
  xlab("V7") +
  ylab("V4 (Reversed)") +
  #scale_fill_manual(values = c("1" = "gray30", "2" = "gray60")) +
  scale_x_continuous(breaks = c(0, 50, 100)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        # panel.grid.major.x = element_line(color = "gray", size = 0.5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 36))
fig1_BUS_rem
ggsave("figure_1_BUSCO_complete_rem.png", plot=fig1_BUS_rem, width=3, height=20, dpi=600)
################################################################################











