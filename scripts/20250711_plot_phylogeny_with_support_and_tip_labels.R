require('ape')
# require('ggtree')

tree <- read.tree("data/diptera.supermatrix.phy.treefile")
rooted_tree <- root(tree, "Panorpa_germanica", resolve.root = T)
rooted_tree$node.label[tree$node.label == "100"] <- ""

pdf('figures/complete_tree_with_confidence.pdf', height = 60, width = 20)
    plot(rooted_tree,  show.node.label = TRUE)
dev.off()
