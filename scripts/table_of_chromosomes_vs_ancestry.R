
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument("-a", "-alg", 
    help="The name of the ALG definitions to be used.")
parser$add_argument("-s", "-species", 
    help="Species to be tabularised.")

args <- parser$parse_args()

source('scripts/20250620_colour_pal.R')

# "Atylotus_latistriatus.syngraph.buscos.tsv"
species_file <- paste0('data/busco_tables/',  args$s, '.syngraph.buscos.tsv')
#species_file_y <- paste0('data/busco_tables/',  args$s, '.syngraph.y.buscos.tsv')

buscos <- read.table(species_file, sep = '\t')
#buscos <- rbind(read.table(species_file, sep = '\t'), read.table(species_file_y, sep = '\t'))
# buscos <- read.table(species_file, sep = '\t')
ALGs <- read.table(args$a)

row.names(ALGs) <- ALGs[, 1]

buscos[, 'ALG'] <- ALGs[buscos[, 1], 2]

sum(is.na(buscos[, 'ALG'])) # filt BUSCOs wo ALG assignment 
buscos <- buscos[!is.na(buscos[, 'ALG']), ]

table(buscos[, 2], buscos[, 'ALG'])