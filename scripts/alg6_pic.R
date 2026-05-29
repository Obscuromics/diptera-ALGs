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
require('reshape2')
source('scripts/20250620_colour_pal.R')

estimate_internal_lengths <- function(tree, tip_lengths) {
  internal_lengths <- ace(tip_lengths, tree, method = "pic")$ace
  return(internal_lengths)
}

tree <- read.tree("data/syngraph/diptera.no_plecia.mindist.m165.newick.txt") #tree
chromosome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?gid=1940964825#gid=1940964825", format='csv'),
                            stringsAsFactors = F, header = T, check.names = F)

alg6_chrs <- chromosome_data %>% 
  filter(alg6_resembling_chromosome == 'YES')

pruned_tree <- drop.tip(
  tree,
  setdiff(tree$tip.label, alg6_chrs$species)
)

chromosome_lengths <- setNames( alg6_chrs$chromosome_size_b, alg6_chrs$species)
chromosome_lengths <- chromosome_lengths[pruned_tree$tip.label]

internal_lengths <- estimate_internal_lengths(pruned_tree, chromosome_lengths)
head(internal_lengths)
