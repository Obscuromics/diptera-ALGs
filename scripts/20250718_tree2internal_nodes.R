suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(argparse))


parser <- ArgumentParser()
# parser$add_argument("-f", "--family", 
    # help="The name of the dipteran family to plot (first letter capitalised; ex.: Sciaridae)")
parser$add_argument("-o", "--output",  
    dest="o", help="Base of the output name (.png will be attached)")
parser$add_argument("-t", "--tree-file",  
    help="tree in nwk format")

args <- parser$parse_args()

tree <- read.tree(args$t) #read.tree("data/syngraph/diptera.pruned2.syngraph_infer.newick.txt")

pdf('figures/syngraph_tree_with_nodes.pdf', height = 60, width = 20)
    plot(tree,  show.node.label = TRUE)
dev.off()

# sort(as.numeric(substr(tree$node.label, 2, nchar(tree$node.label))))

write.table(data.frame(node = tree$node.label[tree$node.label != '']), paste0(args$o,"_internal_nodes.tsv"), quote = F, col.names = F, row.names = F)