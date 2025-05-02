install.packages("ape")
install.packages("ggtree")
install.packages("ggplot2")
install.packages("dplyr")
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
################################################################################
# 08/11/2024

# modify the 315 tree

#load the tree
tree_path <- "C:/Users/julia/Documents/Kamil/trees/315_supermatrix.phy.treefile"
tree <- read.tree(tree_path)

# 11 worst genomes, least BUSCOS on assembled chromosomes
tip_to_remove <- c("GCA_024703675",
                   "GCA_018908115",
                   "GCF_018135715", # outgroup
                   "GCA_016254315",
                   "GCA_016170015",
                   "GCA_037356375",
                   "GCA_917563855", # outgroup
                   "GCA_029228625",
                   "GCA_016170025",
                   "GCF_025200985",
                   "GCF_018902025")

tree_minus11 <- drop.tip(tree, tip_to_remove)

output_file <- "C:/Users/julia/Documents/Kamil/trees/tree_minus11.treefile"

# Baum als .treefile speichern
write.tree(tree_minus11, file = output_file)
# in bash:
# 03_syngraph_all_11

################################################################################
# 11/11/2024

# reintegrate 2 outgroups
# and 4 others
# 1779 GCF_018135715 Dannaus plexippus !! outgroup: Lepidoptera
# 2418 GCA_016254315 Anopheles parensis
# 2433 GCA_016170015 Anopheles longipalpis
# 2457 GCA_037356375 Drosophila_affinis
# 2557 GCA_917563855 Limnephilus_lunatus outgroup: Trichoptera
# 2584 GCA_029228625 Pseudolycoriella hygida

# remove 2 Psychodidae
# GCF_024334085.1 Lutzomyia_longipalpis
# GCF_024763615.1 Phlebotomus_papatasi

tree_path <- "C:/Users/julia/Documents/Kamil/trees/315_supermatrix.phy.treefile"
tree <- read.tree(tree_path)

tip_to_remove <- c("GCA_024703675",
                   "GCA_018908115",
                   #"GCF_018135715",
                   #"GCA_016254315",
                   #"GCA_016170015",
                   #"GCA_037356375",
                   #"GCA_917563855",
                   #"GCA_029228625",
                   "GCA_016170025",
                   "GCF_025200985",
                   "GCF_018902025",
                   "GCF_024334085",
                   "GCF_024763615")
GCA_024703675
GCA_018908115
GCA_016170025
GCF_025200985
GCF_018902025
GCF_024334085
GCF_024763615

tree_minus7 <- drop.tip(tree, tip_to_remove)

output_file <- "C:/Users/julia/Documents/Kamil/trees/tree_minus7.treefile"

# Baum als .treefile speichern
write.tree(tree_minus7, file = output_file)

################################################################################
# baum, in dem der last common ancestor nur 1 ALG hat
# beide Knoten davor aber 5
# Mecoptera outgroup wurde als ingroup betrachtet
# entfernen! vlt ist das das Problem

tree_path <- "C:/Users/julia/Documents/Kamil/trees/tree_minus7.treefile"
tree <- read.tree(tree_path)

tip_to_remove <- c("GCA_963678705")

tree_minus8 <- drop.tip(tree, tip_to_remove)

output_file <- "C:/Users/julia/Documents/Kamil/trees/tree_minus8.treefile"

# Baum als .treefile speichern
write.tree(tree_minus8, file = output_file)
################################################################################
# 19/11/2024
# remove all tips with split chromosomes from the tree tree_minus8

#load the tree
tree_path <- "C:/Users/julia/Documents/Kamil/trees/tree_minus8.treefile"
tree <- read.tree(tree_path)
tree
# 11 worst genomes, least BUSCOS on assembled chromosomes
tip_to_remove <- c("GCF_000001215",
                   "GCF_016746395",
                   "GCF_017639315",
                   "GCF_037355615",
                   "GCF_016746365",
                   "GCF_003369915",
                   "GCF_004382195",
                   "GCF_004382145",
                   "GCA_943734665",
                   "GCF_017562075",
                   "GCF_025231255",
                   "GCF_009650485",
                   "GCF_016746245",
                   "GCA_014170255",
                   "GCF_014743375",
                   "GCF_023558435",
                   "GCF_023558535",
                   "GCA_024703675", # removed previously
                   "GCF_016746235",
                   "GCA_963583835",
                   "GCA_004354385",
                   "GCA_016170035",
                   "GCA_016254315",
                   "GCA_016170015")

tree_minus_split <- drop.tip(tree, tip_to_remove)

output_file <- "C:/Users/julia/Documents/Kamil/trees/tree_minussplit.treefile"

# Baum als .treefile speichern
write.tree(tree_minus_split, file = output_file)
# in bash:
# 03_syngraph_all_11
################################################################################
# remove 18 more species from all_11
# all_29
# they seem to cause distrubance, desited to remove them by looking at tree
tree_path <- "C:/Users/julia/Documents/Kamil/trees/tree_minus11.treefile"
tree <- read.tree(tree_path)

tip_to_remove <- c("GCA_951217065",
                   "GCA_963662105",
                   "GCA_008121215",
                   "GCA_963969585",
                   "GCA_008121275",
                   "GCF_003369915",
                   "GCF_009870125",
                   "GCA_002237135",
                   "GCA_964056695", 
                   "GCA_963971225",
                   "GCA_963668005",
                   "GCA_920105625",
                   "GCA_949629155",
                   "GCF_945859685",
                   "GCA_920104205",
                   "GCA_951509635",
                   "GCA_905146935",
                   "GCA_958431115")

tree_minus_29 <- drop.tip(tree, tip_to_remove)

output_file <- "C:/Users/julia/Documents/Kamil/trees/tree_minus29.treefile"

# Baum als .treefile speichern
write.tree(tree_minus_29, file = output_file)
################################################################################

# do the same as in minus_29 but keep certain species as removal of them removes whole families
# 
tree_path <- "C:/Users/julia/Documents/Kamil/trees/tree_minus11.treefile"
tree <- read.tree(tree_path)

tip_to_remove <- c("GCA_951217065",
                   "GCA_963662105",
                   "GCA_008121215",
                   "GCA_963969585",
                   "GCA_008121275",
                   "GCF_003369915",
                   "GCF_009870125",
                   #"GCA_002237135", Diopsidae
                   #"GCA_964056695", Ulidiidae
                   "GCA_963971225",
                   "GCA_963668005",
                   #"GCA_920105625", Clusiidae
                   "GCA_949629155",
                   "GCF_945859685",
                   "GCA_920104205",
                   "GCA_951509635",
                   "GCA_905146935",
                   "GCA_958431115")

tree_minus_26 <- drop.tip(tree, tip_to_remove)

output_file <- "C:/Users/julia/Documents/Kamil/trees/tree_minus26.treefile"

# Baum als .treefile speichern
write.tree(tree_minus_26, file = output_file)
################################################################################







