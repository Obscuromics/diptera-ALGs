# phylogeny with 312 diptera and 1 mecoptera species
# excluded species marked in red
# separately plot traits: BUSCO completeness, genome size, haploid chromosome number
################################################################################
require('ape')
require('ggtree')
require('ggplot2')
require('dplyr')
require('phytools')
require('gsheet')
require("ggpubr")
require("viridis")
require("patchwork")
require("stringr")
require("tibble")
require('plyr')

################################################################################
#root <- "/Users/ab66/Documents/sanger_work/diptera/diptera-ALGs/"
# root <- paste0(getwd(), "/")
# names <- read.table("data/acc_names_match.txt", sep = ",")
# root <- "/Users/ab66/Documents/sanger_work/diptera/diptera-ALGs/"
#root <- paste0(getwd(), "/")
# tree <- read.tree(paste0(root, "data/diptera.supermatrix.phy.treefile"))
# names <- read.table(paste0(root, "data/acc_names_match.txt"), sep = ",")

# change accessions to species names
# tip_labels <- tree$tip.label
# name_mapping <- setNames(names$V2, names$V1)
# new_tip_labels <- ifelse(tip_labels %in% names(name_mapping), 
#                          name_mapping[tip_labels], 
#                          tip_labels)
# tree$tip.label <- new_tip_labels

# remove two unused outgroups
# tip_to_remove <- c("Danaus_plexippus",
#                    "Limnephilus_lunatus")
# tree_clean <- drop.tip(tree, tip_to_remove)

# TODO: colour species excluded from the analysis
# color 24 species that were excluded
# colour_tips <- c("Myopa_testacea",
#                  "Melieria_omissa",
#                  "Drosophila_pseudoobscura",
#                  "Eupeodes_latifasciatus",
#                  "Syrphus_vitripennis",
#                  "Drosophila_helvetica",
#                  "Drosophila_lowei",
#                  "Baccha_elongata",
#                  "Pseudolycoriella_hygida",
#                  "Drosophila_miranda",
#                  "Drosophila_athabasca",
#                  "Drosophila_affinis",
#                  "Anopheles_parensis",
#                  "Anopheles_longipalpis",
#                  "Anopheles_rivulorum",
#                  "Drosophila_nebulosa",
#                  "Drosophila_gunungcola",
#                  "Drosophila_willistoni",
#                  "Anopheles_vaneedeni",
#                  "Scaeva_pyrastri",
#                  "Eupeodes_luniger",
#                  "Meiosimyza_platycephala",
#                  "Melieria_crassipennis",
#                  "Eupeodes_corollae")

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

#tree <- read.tree("data/diptera.supermatrix.phy.treefile")
tree <- read.tree("../diptera_family_tree/diptera.supermatrix.phy.treefile")

all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?usp=sharing", format='csv'),
                            stringsAsFactors = F, header = T, check.names = F)
genome_data <- all_genome_data[all_genome_data[, 'TO ADD'] == 'KEEP', ]
genome_sizes <- all_genome_data[all_genome_data[, 'TO ADD'] == 'KEEP', c('species', 'genome_size')]

scale_factor <- 8
outgroup <- "Panorpa_germanica"
tip_labels <- tree$tip.label

rooted_tree <- reroot(tree, node.number = which(tip_labels == outgroup), position =  node.depth.edgelength(tree)[which(tip_labels == outgroup)] / scale_factor)
rooted_tree$node.label <- c(1:rooted_tree$Nnode)

# Genome size reconstruction

#manual_to_remove <- c('Anopheles_coustani')
#genome_sizes[genome_sizes$species=='Anopheles_coustani',]$genome_size <- NA
#genome_sizes[genome_sizes$species=='Dexiosonum_caninum',]$genome_size <- NA
#genome_sizes <- add_row(genome_sizes, 'species'='Dexiosoma_caninum', 'genome_size'=NA)
#genome_sizes <- add_row(genome_sizes, 'species'='Dexiosonum_caninum', 'genome_size'=NA)
#species_to_remove <- c(outersect(genome_sizes[complete.cases(genome_sizes), 'species'], rooted_tree$tip.label), manual_to_remove)
#dropped_tree <- drop.tip(rooted_tree, species_to_remove)
#genome_sizes <- genome_sizes %>% filter(species %in% dropped_tree$tip.label)

genome_sizes <- deframe(genome_sizes)
fit <- ace(genome_sizes, rooted_tree, method="pic", type='continuous')

#ace_df <- data.frame(fit$ace)
#colnames(ace_df) <- 'genome_sizes'
#node_vals <-  dplyr::bind_rows(ace_df ,data.frame(genome_sizes))
#node_vals$species <- rownames(node_vals)
s#izes <- deframe(node_vals[,c('species','genome_sizes')])
r#ooted_tree$sizes <- sizes

# edges
#edge_df <- as.data.frame(rooted_tree$edge)
#colnames(edge_df) <- c("node", "child")
#edge_df <- edge_df |>
#  mutate(value = sizes[node])
#edge_df <- edge_df[,c('node','value')]


#tree_show <- ggtree(rooted_tree) %<+% edge_df + aes(color=value) +
#              theme_tree2() + 
#              scale_x_continuous(labels = abs) + 
#              xlab('')+
#              scale_color_viridis_c()

#tree_show <- revts(tree_show)

#ggsave("figures/figure_1_tree.png", plot=tree_show, width=10, height=13, dpi=600)

#nodes
tree_show <- ggtree(rooted_tree)+ 
              theme_tree2() + 
              scale_x_continuous(labels = abs) + 
              xlab('')
              #geom_nodepoint(aes(color = sizes), size = 1) +
              #scale_color_viridis_c()+
              #theme(legend.position = "none")

tree_show <- revts(tree_show)
  #geom_tree(aes(color = label %in% colour_tips), size = 0.5) +
  #geom_tiplab(aes(color = "black"), align = TRUE, offset = 0.05, 
  #            linetype = "dotted", linesize = 0.05, size = 1.2) +
  #scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  #theme(legend.position = "none") 

# safe image
#ggsave("figures/figure_1_tree.png", plot=tree_show, width=10, height=13, dpi=600)

################################################################################
# extract information about the genomes



species_family_tab <- genome_data[, c('species', 'family', 'genome_size', 'n_chromosomes_in_fasta')]
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

# set chromosome number and genome size of the outgroup to 0 to remove from the plot, 
# but keep it so the plots align well with the tree

sorted_data_tips_desc$group <- (cumsum(
  sorted_data_tips_desc$family != dplyr::lag(sorted_data_tips_desc$family, default = sorted_data_tips_desc$family[1])) %% 2) + 1

sorted_data_tips_desc$genome_size[sorted_data_tips_desc$label == outgroup] <- 1
sorted_data_tips_desc$num_chrom[sorted_data_tips_desc$label == outgroup] <- 1
sorted_data_tips_desc$group[sorted_data_tips_desc$label == outgroup] <- "4"

# plot genome size
fig1_size <- ggplot(sorted_data_tips_desc, aes(x = genome_size, y = factor(y, levels = rev(unique(y))), fill = factor(group))) +
  geom_bar(stat = "identity") +
  xlab("Genome size (Mbp)") +
  ylab("V4 (Reversed)") +
  scale_fill_manual(values = c("1" = "gray30", "2" = "gray60", "4"="white")) +
  scale_x_continuous(
    breaks = seq(0, max(sorted_data_tips_desc$genome_size, na.rm = TRUE), by = 1000e6),
    labels = scales::label_number(scale = 1e-6)
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        #panel.grid.major.x = element_line(color = "gray", size = 0.5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x.bottom = element_line(size=0.5),
        # axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        legend.position = "none")

#ggsave("figures/figure_1_size.png", plot=fig1_size, width=2, height=13, dpi=600)

# plot haploid chromosome number
#sorted_data_tips_desc$group[sorted_data_tips_desc$num_chrom==6] <- "3"
#sorted_data_tips_desc <- add_row(sorted_data_tips_desc, 'label'='Panorpa_germanica','num_chrom'=2, 'genome_size'=100e6, 'group'="4")
fig1_hap <- ggplot(sorted_data_tips_desc, aes(x = num_chrom, y = factor(y, levels = rev(unique(y))), fill = factor(group))) +
  geom_bar(stat = "identity") +
  xlab("Chromosomes") +
  ylab("V4 (Reversed)") +
  scale_fill_manual(values = c("1" = "gray30", "2" = "gray60", "4"="white")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        #panel.grid.major.x = element_line(color = "gray", size = 0.5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x.bottom = element_line(size=0.5),
        #axis.title.x = element_blank(),
        legend.position = "none")

#fig1_hap
#ggsave("figures/figure_1_haploidnumber.png", plot=fig1_hap, width=2, height=13, dpi=600)

################################################################################

# arrange three plots together
plt_all <- (tree_show | fig1_size | fig1_hap) +
            plot_layout(widths = c(3, 1, 1))
print('Generating figures/fig1_first_three.svg and figures/fig1_first_three.png')
ggsave("figures/fig1_first_three.png", plot=plt_all, dpi=600, width = 6, height = 10)
ggsave("figures/fig1_first_three.svg", plot=plt_all, dpi=600, width = 6, height = 10)

################################################################################
# plot species to family connector

# reverse the order
sorted_data_tips_desc <- sorted_data_tips_desc %>% filter(!label == outgroup) %>% arrange(y) 
# table(sorted_data_tips_desc[, 'family'])
number_of_species <- nrow(sorted_data_tips_desc)
families <- unique(sorted_data_tips_desc$family)
number_of_families <- length(families)

species_step <- (1 / number_of_species)
species_y <- (1:number_of_species / number_of_species) - (species_step / 2)
names(species_y) <- sorted_data_tips_desc$label

family_step <- (1 / number_of_families)
family_y <- (1:number_of_families / number_of_families) - (family_step / 2)
names(family_y) <- families

# Inpired by code here: https://stackoverflow.com/a/32100956/29281675
curve_manual <- function(x1, y1, x2, y2, scale = 0.08, plot = F, ...){
  line_resolution = 1001
  curve_list <- list()
  curve_list[['x']] <- seq(x1, x2, len = line_resolution)
  curve_list[['y']] <- plogis( seq(0, 1, len = line_resolution), scale = scale, loc = 0.5 ) * (y2 - y1) + y1
  if ( plot ){
    lines(curve_list[['x']], curve_list[['y']], ...)
  }
  return(curve_list)
}

# species richness
diptera_org <- read.table("tables/diptera_org.tsv", header = FALSE)
colnames(diptera_org) <- c("family", "count")
count2color <- data.frame(count = sort(unique(diptera_org$count)), color = viridis(length(unique(diptera_org$count))))
diptera_org <- left_join(diptera_org, count2color, by = "count")
diptera_org <- left_join(diptera_org, count(sorted_data_tips_desc, 'family'), by='family')
diptera_org$prop <- diptera_org$freq / diptera_org$count
prop2color <- data.frame(prop = sort(unique(diptera_org$prop)), prop_color = viridis(length(unique(diptera_org$prop))))
diptera_org <- left_join(diptera_org, prop2color, by = "prop")

sorted_data_tips_desc <- left_join(sorted_data_tips_desc, diptera_org, by = "family")
sorted_data_tips_desc$color[is.na(sorted_data_tips_desc$count)] <- "grey60"

print('Generating figures/figure_1_sp_family_connectors.pdf')
pdf('figures/TEST_figure_1_sp_family_connectors.pdf', width = 4, height = 10)

plot(NULL, xlim = c(0, 2), ylim = c(0, 1), axes = F, xlab = '', ylab = '')
text(1, family_y, families, cex = 0.70)#, pos = 4)

x1 <- 0.25
x2 <- 0.65

x11 <- 1.75
x22 <- 1.35

curviness <- 0.10
for ( i in 1:nrow(sorted_data_tips_desc)){
  fam <- families[i]
  col <- sorted_data_tips_desc$color[sorted_data_tips_desc$family == fam]
  border <- sorted_data_tips_desc$color[sorted_data_tips_desc$family == fam]
  all_species <- sorted_data_tips_desc[sorted_data_tips_desc[, 'family'] == fam, 'label']
  ys1 <- max(species_y[all_species])
  ys2 <- min(species_y[all_species])
  yf <- family_y[fam]
  top_curve <- curve_manual(x1, ys1 + (species_step / 6), x2, yf, scale = curviness)
  bot_curve <- curve_manual(x1, ys2 - (species_step / 6), x2, yf, scale = curviness)
  polygon(c(top_curve[['x']], rev(bot_curve[['x']])), c(top_curve[['y']], rev(bot_curve[['y']])), 
          col = col, border = col)
  # rect(0.37, ys2 - (species_step / 2), 0.39, ys1 + (species_step / 2), col = col[1], border = NA) # , bty = 'n'
}

for ( i in 1:nrow(sorted_data_tips_desc)){
  fam <- families[i]
  col <- sorted_data_tips_desc$prop_color[sorted_data_tips_desc$family == fam]
  border <- sorted_data_tips_desc$color[sorted_data_tips_desc$family == fam]
  all_species <- sorted_data_tips_desc[sorted_data_tips_desc[, 'family'] == fam, 'label']
  ys1 <- max(species_y[all_species])
  ys2 <- min(species_y[all_species])
  yf <- family_y[fam]
  top_curve <- curve_manual(x11, ys1 + (species_step / 6), x22, yf, scale = curviness)
  bot_curve <- curve_manual(x11, ys2 - (species_step / 6), x22, yf, scale = curviness)
  polygon(c(top_curve[['x']], rev(bot_curve[['x']])), c(top_curve[['y']], rev(bot_curve[['y']])), 
          col = col, border = col)
  # rect(0.37, ys2 - (species_step / 2), 0.39, ys1 + (species_step / 2), col = col[1], border = NA) # , bty = 'n'
}

dev.off()

#write.table(sorted_data_tips_desc$label, file = "labels_sorted.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

################################################################################