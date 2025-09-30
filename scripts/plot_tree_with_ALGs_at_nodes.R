suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ape))
# suppressPackageStartupMessages(library(ggplot2))

parser <- ArgumentParser()

# parser$add_argument("-a", "--alg", 
#     help="The name of the ALG definitions to be used.")
# parser$add_argument("-o", "--output",  
#     dest="o", help="Base of the output name (.png will be attached)")

args <- parser$parse_args()
source('scripts/20250620_colour_pal.R')

# these need to be made into arguments in the end
treefile <- 'data/syngraph/diptera.no_plecia.mindist.m165.newick.txt'
rearrangement_file <- 'data/syngraph/diptera.no_plecia.mindist.m165.rearrangements.tsv'

# treefile <- 'data/syngraph/diptera.no_plecia.mindist.m165.newick.txt'
busco_asn_file <- 'tables/ALGs_syngraph_diptera.tsv' # args$a
output <- 'figures/syngraph_tree_of_changes.pdf'
# simplified; if 0 means not simplified, if a number it is the number of pieces it gets simplfied to
simplified <- 0
just_internal <- TRUE

directory_with_nodes <- paste0("data/syngraph/node_assignments/")
all_node_files <- dir(directory_with_nodes)

ALG_data <- read.table(busco_asn_file, col.names = c('busco', 'chrQ'))

########## Loading states at all the nodes
# testing: input_file <- paste0("data/syngraph/node_assignments/", all_node_files[1])
load_node <- function(node_file){
    node_lgs <- read.table(node_file, header = T, col.names = c('busco', 'chrNode'))
    # needs ALG_data
    node_lgs_by_algs <- merge(ALG_data, node_lgs)
    dat <- data.frame(table(node_lgs_by_algs$chrQ,node_lgs_by_algs$chrNode))
    names(dat) <- c("ALG","node","Count")
    return(dat)
}

all_node_asn <- lapply(paste0(directory_with_nodes, all_node_files), load_node)
names(all_node_asn) <- sapply(strsplit(all_node_files, '_'), function(x){x[1]})
# all_node_asn is a list node -> busco assignments gruped by assigments

######### Processing rearrangement file

all_rearrangements <- read.table(rearrangement_file, sep = '\t')
nodes_with_changes <- list()
nodes_with_changes[names(table(all_rearrangements[, 2]))] <- TRUE

# if(is.null(nodes_with_changes[['n84']])){ print('well hellp')}
# is.null(nodes_with_changes[['n583']])

######### Loading and plotting the tree
tree <- read.tree(treefile)
rooted_tree <- tree # root(tree, "Panorpa_germanica", resolve.root = T)
rooted_tree$node.label[tree$node.label == "100"] <- ""

# read.tree(treefile)
# fakt_tree <- rcoal(10)
# tree_plot <- plot(tree, type = 'n')
# tree_plot$xx

# fakt_tree$node.label <- tree$node.label[2:11]

# plot(fakt_tree,  show.node.label = TRUE)

pdf(output, height = 60, width = 20)

plot(rooted_tree,  show.node.label = TRUE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

chwidth <- max(lastPP$x.lim) / 100
chheight <- 0.65
chgap <- max(lastPP$x.lim) / 400
y_smallstep <- 0.04

for (node_i in c(length(rooted_tree$tip.label) + 2):length(lastPP$xx)){
    node <- c(rooted_tree$tip.label, rooted_tree$node.label)[node_i]
    if ( is.null(nodes_with_changes[[node]]) ){
        print(paste("Skipping", node, "no changes."))
        next
    }
    node_asn <- all_node_asn[[node]]
    x <- lastPP$xx[node_i]
    y <- lastPP$yy[node_i]
    lgs_to_plot <- levels(node_asn[, 'node'])

    yfrom <- y + y_smallstep
    yto <- y + y_smallstep + chheight
    for (chr in 1:length(lgs_to_plot)){
        xfrom <- x - (chr * (chwidth + chgap))
        xto <- x - ((chr - 1) * (chwidth + chgap) + chgap)

        bar_subset <- node_asn[node_asn[, 'node'] == lgs_to_plot[chr], ]
        # making sure the order in the plot will be consisent
        bar_subset <- bar_subset[order(bar_subset[, 'ALG'], decreasing = F), ]
        if (simplified == 0){
            columns_to_plot <- y + y_smallstep + (chheight * c(0, cumsum(bar_subset[, "Count"])) / sum(bar_subset[, "Count"]))
        }

        for(i in 1:nrow(bar_subset)){
            rect(xfrom, columns_to_plot[i], 
                xto, columns_to_plot[i + 1], 
                col = pal[as.character(bar_subset[i, 'ALG'])], bty = 'n')
        }
    }
}

dev.off()