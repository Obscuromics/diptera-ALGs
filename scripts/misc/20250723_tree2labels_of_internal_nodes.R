suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(gsheet))


parser <- ArgumentParser()
# parser$add_argument("-f", "--family", 
    # help="The name of the dipteran family to plot (first letter capitalised; ex.: Sciaridae)")
parser$add_argument("-o", "--output",  
    dest="o", help="Base of the output name (.png will be attached)")
parser$add_argument("-t", "--tree-file",  
    help="tree in nwk format")
# parser$add_argument("-l", "--level",  
#     help="What level the nodes should be extracted.")

args <- parser$parse_args()

tree <- read.tree(args$t)
# tree <- read.tree("data/syngraph/diptera.no_chiro.syngraph_infer.min_dist.m165.newick.txt")
# tree <- read.tree("data/syngraph/diptera.juliaset.newick")

tips_nodes <- c(tree$tip.label, tree$node.label)

all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?usp=sharing", format='csv'),
                            stringsAsFactors = F, header = T, check.names = F)

all_genome_data <- all_genome_data[all_genome_data[, 'TO ADD'] %in% c('KEEP', 'OUTGROUP'), ]
species_in_the_tree <- sapply(all_genome_data[, 'species'], function(x){x %in% tree$tip.label})

tree_genome_data <- all_genome_data[species_in_the_tree, ]
# SHOULD GIVE STATS HERE ABOUT THIS FILTERING
# SHOULD TEST HERE FOR TIPS MISSING IN THE TABLE

# tree_genome_data[, c('species', 'family')]
all_brachycera <- tree_genome_data[tree_genome_data[, 'suborder'] == "Brachycera", 'species']

if(is.monophyletic(tree, all_brachycera)){
    tips_nodes[getMRCA(tree, all_brachycera)]
} else {
    print('Brachycera are not monophyletic')
}

# all_Asiloidea <- tree_genome_data[tree_genome_data[, 'superfamily'] == 'Asiloidea', 'species']
# if(is.monophyletic(tree, all_Asiloidea)){
#     tips_nodes[getMRCA(tree, all_Asiloidea)]
# } else {
#     print('Asiloidea are not monophyletic')
# }
tips_nodes[nodepath(tree, 1, 3)]
tips_nodes[getMRCA(tree, tips_nodes[c(1,3)])]
all_Syrphoidea <- tree_genome_data[tree_genome_data[, 'superfamily'] == "Syrphoidea", 'species']

# if(is.monophyletic(tree, all_Syrphoidea)){
#     tips_nodes[getMRCA(tree, all_brachycera)]
# } else {
#     print('Brachycera are not monophyletic')
# }

all_sciaridae <- tree_genome_data[tree_genome_data[, 'family'] == "Sciaridae", 'species']
tips_nodes[getMRCA(tree, all_sciaridae)]
extract.clade(tree, "n99")

# sort(as.numeric(substr(tree$node.label, 2, nchar(tree$node.label))))
extract.clade(tree, "n42")


write.table(data.frame(node = tree$node.label[tree$node.label != '']), paste0(args$o,"_internal_nodes.tsv"), quote = F, col.names = F, row.names = F)


fake_tree <- read.tree('data/syngraph/fake_tree.nwk')

fake_tree <- read.tree(text = '((Anastrepha_ludens:0.00499,Anastrepha_obliqua:0.00455)n463:0.0363,(Anomoia_purmunda:0.03743,Philophylla_caesio:0.04317)n464:0.01056);')

c(fake_tree$tip.label, fake_tree$node.label)[getMRCA(fake_tree, c("Anastrepha_ludens", "Anastrepha_obliqua"))]

