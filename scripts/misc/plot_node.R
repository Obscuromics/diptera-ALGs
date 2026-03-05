suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))

parser <- ArgumentParser()
parser$add_argument("-n", "--node", 
    help="The name of the node to plot.")
parser$add_argument("-a", "--alg", 
    help="The name of the ALG definitions to be used.")
parser$add_argument("-o", "--output",  
    dest="o", help="Base of the output name (.png will be attached)")


args <- parser$parse_args()
source('scripts/20250620_colour_pal.R')
names(pal)[1:6] <- paste0('d', 1:6)

# node_to_plot <- "n14"
node_to_plot <- args$n
busco_asn_file <- args$a # 'data/diptera.no_plecia.mindist.m165_n1_n2.tsv' # 'Diptera_ALG_new.tsv'
# output_file <- 'data/syngraph/node_plots/n1.png'
output_file <-  paste0(args$o, '.png')
input_file <- paste0("data/syngraph/node_assignments/", node_to_plot, "_asgn.tsv")

ALG_data <- read.table(busco_asn_file, col.names = c('busco', 'chrQ'))
node_lgs <- read.table(input_file, header = T, col.names = c('busco', 'chrNode'))

node_lgs_by_algs <- merge(ALG_data, node_lgs)

# <- table(node_lgs_by_algs[, 'chrNode'])

# dat <- data.frame(table(df$Fruit,df$Bug))
# names(dat) <- c("Fruit","Bug","Count")
dat <- data.frame(table(node_lgs_by_algs$chrQ,node_lgs_by_algs$chrNode))
names(dat) <- c("ALG","node","Count")

nodes_to_plot <- levels(dat[, 'node'])
barwidth <- 0.4

png(output_file)
    plot(NA, xlim = c(1 - barwidth, length(nodes_to_plot) + barwidth), ylim = c(0, 1), bty = 'n', axes = F, xlab = '', ylab = '')

    for (node in 1:length(nodes_to_plot)){
        xpos <- node
        bar_subset <- dat[dat[, 'node'] == nodes_to_plot[node], ]
        # making sure the order in the plot will be consisent
        bar_subset <- bar_subset[order(bar_subset[, 'ALG'], decreasing = F), ]
        columns_to_plot <- c(0, cumsum(bar_subset[, "Count"])) / sum(bar_subset[, "Count"])

        for(i in 1:nrow(bar_subset)){
            rect(xpos - barwidth, columns_to_plot[i], 
                xpos + barwidth, columns_to_plot[i + 1], 
                col = pal[as.character(bar_subset[i, 'ALG'])], bty = 'n')
        }
    }
dev.off()