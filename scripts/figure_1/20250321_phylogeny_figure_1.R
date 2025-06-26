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
#root <- "/Users/ab66/Documents/sanger_work/diptera/diptera-ALGs/"
root <- getwd()
tree <- read.tree(paste0(root, "data/315_supermatrix.phy.treefile"))
names <- read.table(paste0(root, "data/acc_names_match.txt"), sep = ",")

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

tree_show <- ggtree(tree_clean)
  #geom_tree(aes(color = label %in% colour_tips), size = 0.5) +
  #geom_tiplab(aes(color = "black"), align = TRUE, offset = 0.05, 
  #            linetype = "dotted", linesize = 0.05, size = 1.2) +
  #scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  #theme(legend.position = "none") 
tree_show

# safe image
#ggsave(paste0(root, "figures/tree_figure_1.png"), plot=tree_show, width=10, height=13, dpi=600)

################################################################################
# extract information about the genomes
require('gsheet')
all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?usp=sharing", format='csv'),
                            stringsAsFactors = F, header = T, check.names = F)
genome_data <- all_genome_data[all_genome_data[, 'TO ADD'] == 'KEEP', ]

species_family_tab <- genome_data[, c('species', 'family', 'genome_size', 'num_chrom')]
colnames(species_family_tab) <- c("label", "family", 'genome_size', 'num_chrom')

# extract correct order of species
data <- tree_show$data
sorted_data_tips_desc <- data %>% filter(isTip == TRUE) %>% arrange(-y)

# add family names
sorted_data_tips_desc <- left_join(sorted_data_tips_desc, species_family_tab, by = "label")
sorted_data_tips_desc <- sorted_data_tips_desc %>% filter(!is.na(family))
sorted_data_tips_desc <- as.data.frame(sorted_data_tips_desc)

################################################################################
# plot chromosome number and genome size

sorted_data_tips_desc$group <- (cumsum(
  sorted_data_tips_desc$family != dplyr::lag(sorted_data_tips_desc$family, default = sorted_data_tips_desc$family[1])) %% 2) + 1

# plot haploid chromosome number
fig1_hap <- ggplot(sorted_data_tips_desc, aes(x = num_chrom, y = factor(y, levels = rev(unique(y))), fill = factor(group))) +
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
        legend.position = "none") +
  annotate(xmin = 6, xmax = 7, )

#fig1_hap

#ggsave(paste0(root, "figures/figure_1_haploidnumber.png"), plot=fig1_hap, width=2, height=13, dpi=600)

# plot genome size
fig1_size <- ggplot(sorted_data_tips_desc, aes(x = genome_size, y = factor(y, levels = rev(unique(y))), fill = factor(group))) +
  geom_bar(stat = "identity") +
  xlab("V6 (Millionen)") +
  ylab("V4 (Reversed)") +
  scale_fill_manual(values = c("1" = "gray30", "2" = "gray60")) +
  scale_x_continuous(
    breaks = seq(0, max(sorted_data_tips_desc$genome_size, na.rm = TRUE), by = 500e6),
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

#fig1_size

#ggsave(paste0(root, "figures/figure_1_size.png"), plot=fig1_size, width=2, height=13, dpi=600)

################################################################################

# arrange three plots together
plt_all <- ggarrange(tree_show + fig1_size + fig1_hap, nrow = 1)
ggsave(paste0(root, "figures/fig1_first_three.png"), plot=plt_all, dpi=600)
ggsave(paste0(root, "figures/fig1_first_three.svg"), plot=plt_all, dpi=600)

################################################################################
# plot species to family connector

# reverse the order
sorted_data_tips_desc <- sorted_data_tips_desc %>% arrange(y)

number_of_species <- nrow(sorted_data_tips_desc)
families <- unique(sorted_data_tips_desc$family)
number_of_families <- length(families)

species_y <- 1:number_of_species / number_of_species
names(species_y) <- sorted_data_tips_desc$label

family_y <- 1:number_of_families / number_of_families
names(family_y) <- families

# code fgrid# code from https://stackoverflow.com/a/32100956/29281675
curveMaker <- function(x1, y1, x2, y2, ...){
  curve( plogis( x, scale = 0.08, loc = (x1 + x2) /2 ) * (y2-y1) + y1, 
         x1, x2, add = TRUE, ...)
}

pdf(paste0(root, 'figures/figure_1_sp_family_connectors.pdf'), width = 4, height = 10)

plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = F, xlab = '', ylab = '')
text(0.60, family_y, families, cex = 0.50, pos = 4)


for ( i in 1:nrow(sorted_data_tips_desc)){ 
  curveMaker(0.05, species_y[sorted_data_tips_desc[i, 'label']], 0.6, family_y[sorted_data_tips_desc[i, 'family']])
}

dev.off()

#write.table(sorted_data_tips_desc$label, file = "labels_sorted.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

################################################################################