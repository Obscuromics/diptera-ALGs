require('ggtree')
require('dplyr')
require('ggplot2')

dot_overview <- read.table('data/species_and_ALG6_chromosomes.tsv')
phylo_tree_file <- 'data/diptera.supermatrix.phy.treefile'

tree <- read.tree(phylo_tree_file)
tree_show <- ggtree(tree)
data <- tree_show$data
sorted_data_tips_desc <- data %>% filter(isTip == TRUE) %>% arrange(-y)

# add family names
colnames(dot_overview)[1] <- "label"
sorted_data_tips_desc <- left_join(sorted_data_tips_desc, dot_overview, by = "label")
sorted_data_tips_desc <- as.data.frame(sorted_data_tips_desc)

################################################################################
# plot chromosome number and genome size

# set chromosome number and genome size of the outgroup to 0 to remove from the plot, 
# but keep it so the plots align well with the tree

outgroup <- "Panorpa_germanica"
sorted_data_tips_desc$genome_size[sorted_data_tips_desc$label == outgroup] <- 0
sorted_data_tips_desc$num_chrom[sorted_data_tips_desc$label == outgroup] <- 0

# sorted_data_tips_desc$group <- (cumsum(
#   sorted_data_tips_desc$family != dplyr::lag(sorted_data_tips_desc$family, default = sorted_data_tips_desc$family[1])) %% 2) + 1

fig1_size <- ggplot(sorted_data_tips_desc, aes(x = size, y = factor(y, levels = rev(unique(y))))) +
  geom_bar(stat = "identity") +
  xlab("V6 (Millionen)") +
  ylab("V4 (Reversed)") +
  scale_x_continuous(
    breaks = seq(0, max(sorted_data_tips_desc$size, na.rm = TRUE), by = 500e6),
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

pdf('figures/dot_sizes.pdf', height = 30, width = 3)
    fig1_size
dev.off()

pdf('figures/dot_sizes_base.pdf', height = 26, width = 3)
    barplot(rev(sorted_data_tips_desc[, 'size']) / 1e6, horiz = T, col = 'black')
dev.off()


