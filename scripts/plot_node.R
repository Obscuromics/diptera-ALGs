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

# node_to_plot <- "n94"
node_to_plot <- args$n
busco_asn_file <- args$a # 'data/ALG_assignments_pruned_m100.tsv' # 'Diptera_ALG_new.tsv'
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

png(output_file)
    ggplot(data=dat, aes(x=node, y=Count, fill=ALG)) + geom_bar(stat="identity") + scale_fill_manual(values=pal)
dev.off()