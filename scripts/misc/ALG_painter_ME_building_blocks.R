# colour genomes with ALGs
# bins of 20 BUSCOs
################################################################################
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(gsheet))
suppressPackageStartupMessages(library(ggtree))

################################################################################
# function to make chromosomes coherent without steps at the ends
clean_data <- function(data) {
  cleaned_values <- c()
  count <- 0
  for (i in seq_along(data$chrQ.y)) {
    value <- data$chrQ.y[i]
    if (value == "100") {
      if (count < 20 && count > 0 && cleaned_values[length(cleaned_values)] != "100") {
        cleaned_values <- head(cleaned_values, -count)
        count <- 0
      }
      count <- 0
    } else {
      count <- count + 1
    }
    cleaned_values <- c(cleaned_values, value)
    if (count == 20) {
      count <- 0
    }
  }
  data.frame(chrQ.y = cleaned_values)
}

# Agument: Take a family name, and create a painting of that using the last
# Requierements
#  - file with BUSCO assignments
#  - direcotory with all the BUSCO runs
        # -> should crosscheck if the files are there
#  - ???

parser <- ArgumentParser()
parser$add_argument("-f", "--family", 
    help="The name of the dipteran family to plot (first letter capitalised; ex.: Sciaridae)", default = "")
parser$add_argument("-l", "--list_of_species", 
    help="A name of a text file with a species list (one per line, underscores between genus and species name)", default = "")
parser$add_argument("-o", "--output",  
    dest="o", help="Base of the output name (.png will be attached)")
parser$add_argument("-a", "--alg-set", dest = "a", 
    help="File with BUSCO to ALG assignments", default = 'data/ALG_assignments_pruned2_m100.tsv')
parser$add_argument("-p", "--paint_by",  
    help="What reference should be used for painting", default = 'ALGs')
parser$add_argument("-s", "--subsample", default = 0,
    help="If there is more than -s genomes selected, subsample to this number")
parser$add_argument("--keep-order", action="store_true", default=FALSE,
    dest="keep", help="Use the order of species as in the list/table.")
parser$add_argument("--tree-order", action="store_true", default=FALSE,
    dest="tree_order", help="Order species according to the tree (instead of alphabetically).")
parser$add_argument("--alg-sort", action="store_true", default=FALSE,
    dest="alg_sort", help="Order chromsomes by their ALG content.")

args <- parser$parse_args()

family <- args$f
species_list_file <- args$l 
busco_asn_file <- args$a # 'Diptera_ALG_new.tsv'
output_file <-  paste0(args$o, '.png')
output_tree <-  paste0(args$o, '_tree.png')
chromosome_files_dir <- 'data/diptera_chromosome_files/'
busco_files_dir <- "data/busco_tables/"

if (family == "" && species_list_file == ""){
  stop("Either a list of species, or family needs to be specified.")
}

################################################################################
ALG_data <- read.table(busco_asn_file, col.names = c('busco', 'chrQ'))
# ALG_data <- read.table("tables/ALGs_syngraph_diptera.tsv", col.names = c('busco', 'chrQ'))
# head(ALG_data)
rownames(ALG_data) <- ALG_data[, 'busco']

################################################################################
# all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?usp=sharing", format='csv'),
#                             stringsAsFactors = F, header = T, check.names = F)
all_genome_data <- read.table('tables/supplementary_table_1.tsv', sep = '\t', header = T, check.names = F)

all_genome_data <- all_genome_data[all_genome_data[, 'TO ADD'] %in% c('KEEP', 'OUTGROUP'), ]


print(paste("Loaded ", nrow(all_genome_data), "genomes."))

family_table <- all_genome_data
# family <- "Conopidae"

if (args$tree_order){
  source("scripts/supp_genome_paintings_per_family/genome_paintings_per_family.R")
  species_to_plot <- all_family_dfs[[family]][, 'species']

  row.names(family_table) <- family_table[, 'species']

  not_in_tree <- family_table[grepl(family, family_table[, 'family']), 'species']
  not_in_tree <- not_in_tree[!not_in_tree %in% species_to_plot]

  family_table <- family_table[c(species_to_plot, not_in_tree), ]

} else {
  if ( nchar(family) != 0 ){
    family_table <- family_table[grepl(family, family_table[, 'family']), ]
    print(paste("Plotting ", nrow(family_table), " members of", family, 'family'))
  }

  if ( nchar(species_list_file) != 0 ){
    species_to_plot <- read.table(species_list_file, header = F)[, 1]
    
    row.names(family_table) <- family_table[, 'species']

    family_table <- family_table[species_to_plot, ]
    print(paste("Subsetting to ", nrow(family_table), " species in the specified list"))
  }

  if (nrow(family_table) > args$s & args$s != 0){
    family_table <- family_table[sample(1:nrow(family_table), args$s), ]
  }

  if ( args$keep | args$tree_order){
    family_table <- family_table[nrow(family_table):1, ]
  } else {
    family_table <- family_table[order(family_table[, 'species'], decreasing = TRUE), ]
  }
}



accesions_to_plot <- family_table[, 'accession']
species_to_plot <- family_table[, 'species']
accesions_wo_suffix <- sapply(strsplit(accesions_to_plot, '[.]'), function(x){x[1]})
# busco_files <- paste0(accesions_wo_suffix, '.tsv')
chrom_files <- paste0(accesions_to_plot, '.chromosomes.txt')
busco_files <- paste0(species_to_plot, '.syngraph.buscos.tsv')

################################################################################
source("scripts/20250620_colour_pal.R")

new_colnames <- c("busco", "chrQ", "Qstart", "Qend")

#### 

available_buscos <- dir(busco_files_dir)
available_chrom_files <- dir(chromosome_files_dir)

present_busco_files <- busco_files %in% available_buscos
present_chromosome_files <- chrom_files %in% available_chrom_files

keep_following <- present_busco_files & present_chromosome_files

if(any(!keep_following)){
  print(paste('Skipping ', paste(accesions_to_plot[!keep_following]), ' - missing busco file!'))
}

busco_files <- busco_files[keep_following]
chrom_files <- chrom_files[keep_following]
family_table <- family_table[keep_following, ]

################################################################################

plot_list <- list()
for (j in 1:length(busco_files)) {
  
  busco_file_path <- file.path(busco_files_dir, busco_files[j])
  chrom_file_path <- file.path(chromosome_files_dir, chrom_files[j])
  
  # Load the chromosome data
  chrom_data <- read.table(chrom_file_path, sep = "\t", header = F, stringsAsFactors = FALSE)

  # Load the busco data
  busco_data <- read.table(busco_file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(busco_data) <- new_colnames
  busco_data <- busco_data[busco_data[, 2] %in% chrom_data[, 1], ] # Remove BUSCOs on non-chromosmal scaffolds

  # Match ALG identity of BUSCOs to genome of interest
  # if ( args$p == 'muller' ){
  #   busco_data_ALG <- merge(busco_data, muller_blocks[, c("busco", "chrQ")], by = "busco", all.x = TRUE)
  # } else {
    busco_data_ALG <- merge(busco_data, ALG_data[, c("busco", "chrQ")], by = "busco", all.x = TRUE)
  # }
  busco_data_ALG <- busco_data_ALG[!is.na(busco_data_ALG$chrQ.y), ] # Remove unassigned BUSCOs
  
  # print("In BUSCOs in the genome")
  # print(table(busco_data_ALG[, "chrQ.y"]))
  # print("In BUSCOs assigned to ALGs")
  # print(table(ALG_data[, c("chrQ")]))
  # print("Proportion")
  # print(round(table(busco_data_ALG[, "chrQ.y"]) / table(ALG_data[, c("chrQ")]), 2))
  
  
  

  if ( args$alg_sort ){
    per_chr_weighted_order <- colSums(table(busco_data_ALG[, "chrQ.y"],  busco_data_ALG[, "chrQ.x"] ) * 1:length(table(busco_data_ALG[, "chrQ.y"]))) / colSums(table(busco_data_ALG[, "chrQ.y"],  busco_data_ALG[, "chrQ.x"] ))
    chr_order <- names(per_chr_weighted_order)[order(per_chr_weighted_order)]
  } else {
    # Sort data by order of chromosomes
    # chromosomes <- unique(busco_data_ALG$chrQ.x)
    # chromosomes <- sort(chromosomes)
    # chromosomes <- chromosomes[order(chrom_data$order)]
    # chrom_data <- chrom_data[order(chrom_data$order), ]

    chr_order <- chrom_data[, 1]
  }

  busco_data_ALG_sorted <- busco_data_ALG[order(factor(busco_data_ALG$chrQ.x, levels = chr_order)), ]
  
  # Insert 80 empty BUSCOs between chromosomes
  busco_data_ALG_sorted <- busco_data_ALG %>%
    mutate(chrQ.x = factor(chrQ.x, levels = chr_order)) %>%
    arrange(chrQ.x, Qstart)
  # remove start and end column
  busco_data_ALG_sorted <- busco_data_ALG_sorted %>%
    select(-Qstart, -Qend)
  # insert 80 empty BUSCOs between chromosomes
  result <- list()
  for (i in 1:(nrow(busco_data_ALG_sorted) - 1)) {
    result <- append(result, list(busco_data_ALG_sorted[i, ]))
    if (busco_data_ALG_sorted$chrQ.x[i] != busco_data_ALG_sorted$chrQ.x[i + 1]) {
      new_rows <- data.frame(
        busco = rep(100, 80),
        chrQ.x = rep(busco_data_ALG_sorted$chrQ.x[i + 1], 80),
        chrQ.y = rep(100, 80)
      )
      result <- append(result, list(new_rows))
    }
  }
  result <- append(result, list(busco_data_ALG_sorted[nrow(busco_data_ALG_sorted), ]))
  busco_data_ALG_final <- do.call(rbind, result)
  # head(busco_data_ALG_final, 100)
  # tail(busco_data_ALG_final, 100)
  # remove everything but chrQ.y
  busco_data_ALG_final <- busco_data_ALG_final %>%
    select(-busco, -chrQ.x)
  # head(busco_data_ALG_final, 100)
  
  # clean edges of chromosomes
  cleaned_data <- clean_data(busco_data_ALG_final)
  busco_data_ALG_final = cleaned_data
  
  # Create bins
  value_counts <- table(busco_data_ALG_final$chrQ.y)
  bin_vector <- vector("integer", length = length(busco_data_ALG_final$chrQ.y))
  
  current_bin <- 1
  counter <- 1
  
  for (i in 1:length(busco_data_ALG_final$chrQ.y)) {
    if (i > 1 && busco_data_ALG_final$chrQ.y[i] == "100" && busco_data_ALG_final$chrQ.y[i-1] != "100") {
      current_bin <- current_bin + 1
      counter <- 1
    }
    
    bin_vector[i] <- current_bin
    counter <- counter + 1
    
    if (counter > 20) {
      current_bin <- current_bin + 1
      counter <- 1
    }
  }
  
  busco_data_ALG_final$bin <- bin_vector
  df_aggregated <- as.data.frame(table(busco_data_ALG_final$bin, busco_data_ALG_final$chrQ.y))
  colnames(df_aggregated) <- c("bin", "value", "count")
  # if(args$p == 'muller'){
  #   pal <- block_pal
  # } else {
    pal <- alg_pal
  # }
  df_aggregated$fill_color <- factor(df_aggregated$value, levels = names(pal))
  
  # Generate the plot
  p <- ggplot(df_aggregated, aes(x = factor(bin), y = count, fill = value)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = pal) +
    labs(x = NULL, y = NULL, title = family_table[j, 'species']) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.title = element_text(size=16)
    )
  
  # Append the plot to the plot_list
  plot_list[[j]] <- p
}

plot_list <- rev(plot_list)
################################################################################
# grid.arrange(grobs = plot_list, ncol = 1)
ggsave(output_file, plot = grid.arrange(grobs = plot_list, ncol = 1), width = 15, height = min(2 * length(plot_list), 40))
################################################################################

# library(tidyverse)
# library(ape)
# library(ggtree)

# p_clean <- ggtree(family_trees[[family]], aes(x, y)) + theme_tree()

# p <- ggtree(family_trees[[family]], aes(x, y), size=2) +
#   geom_tiplab() + theme_tree()

# ggsave(output_tree, plot = p, width = 5, height = min(2 * length(plot_list), 40))

# facet_plot(
#   p,
#   panel = "Paint",
#   data = tip_data,
#   geom_col,
#   mapping = aes(x = value),  # <- no y aesthetic here
#   width = 0.8
# )
# png(output_tree) # , width = 15, height = min(2 * length(plot_list), 40)

#   ggplot(family_trees[[family]], aes(x, y)) + geom_tree() + theme_tree()
  

# dev.off()