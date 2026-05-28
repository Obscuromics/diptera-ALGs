
#_______________________________________________________________________________
library(ape)
library(dplyr)
library(ggtree)
library(zoo)
#_______________________________________________________________________________

tree_file <- "/Users/jg40/Documents/diptera_centromeres/tree/supermatrix.phy.treefile"
tree_diptera <- read.tree(tree_file)

diptera_species <- "/Users/jg40/Documents/diptera_centromeres/tree/20260119_all_species.tsv"
diptera_species <- read.delim(diptera_species, sep = "\t", header = TRUE)

#________________________________species________________________________________

keep_species <- c(
  "Leptarthrus_brevirostris", #1
  "Thereva_nobilitata", #2
  "Empis_livida", #3
  "Nephrocerus_scutellatus", #9
  "Trigonometopus_frontalis", #10
  "Coelopa_pilipes", #11
  "Dryomyza_anilis", #12
  "Conops_quadrifasciatus", #13
  "Pherbina_coryleti", #14 with cen
  "Anthomyza_gracilis", #15
  "Suillia_variegata", #17
  "Neria_commutata", #19
  "Bactrocera_dorsalis", #21 with cen
  "Megamerina_dolium", #23
  "Nemopoda_nitidula", #25
  "Phortica_okadai", #28
  "Drosophila_melanogaster", #28
  "Crataerina_pallida", #29
  "Hydrotaea_cyrtoneurina", #31
  "Scathophaga_stercoraria", #32
  "Eustalomyia_histrio", #33
  "Sarcophaga_rosellei", #34 with cen
  "Epicampocera_succincta", #35
  "Pollenia_amentaria", #36 with cen
  "Phyto_melanocephala", #37 with cen
  "Bellardia_pandia", #38 with cen
  "Stomorhina_lunata" #39 with cen
)

#_______________________________buid_tree_rotate_tips___________________________

common_species <- intersect(tree_diptera$tip.label, keep_species)
tree_overview <- keep.tip(tree_diptera, common_species)
overview_tips <- tree_overview$tip.label
tree_overview$tip.label <- gsub("_", " ", tree_overview$tip.label)

tree_depth <- max(node.depth.edgelength(tree_overview))
plot(
  tree_overview,
  cex = 1.3,
  edge.width = 1.9,
  align.tip.label = TRUE,
  no.margin = TRUE,
  x.lim = c(0, tree_depth * 3.0)
)

p2 <- ggtree(tree_overview, layout = "roundrect") +
  geom_tiplab(size = 6, align = TRUE, fontface = "italic") +
  xlim(0, tree_depth * 3.0) +
  theme_tree()
p2 <- rotate(p2, 38)

p2

tip_order <- p2$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label)

tip_order_raw <- gsub(" ", "_", tip_order)

tree_overview_accession_ggtree <- data.frame(
  species   = tip_order,
  accession = diptera_species$accession[
    match(tip_order_raw, diptera_species$species)
  ],
  stringsAsFactors = FALSE
)

tree_overview_accession_ggtree_rev <- tree_overview_accession_ggtree[nrow(tree_overview_accession_ggtree):1, ]
rownames(tree_overview_accession_ggtree_rev) <- NULL

#_____________________generate_raw_data_for_plots_______________________________

file_path_chrom_names <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/chromosome_names_max12_long_filtered_renamed.tsv"
dir_path_busco <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes"

chromosome_name_len <- read.table(
  file_path_chrom_names,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)
colnames(chromosome_name_len) <- c("chromosome", "length")

# count chromosomes per accession
chrom_df <- chromosome_name_len
is_header <- grepl("^chromosome_GC[AF]_", chrom_df$chromosome)
chrom_df$accession <- NA
chrom_df$accession[is_header] <- sub("^chromosome_", "", chrom_df$chromosome[is_header])
chrom_df$accession <- na.locf(chrom_df$accession)
chrom_counts <- chrom_df[!is_header & chrom_df$chromosome != "", ]
chrom_counts$accession_short <- sub("\\.\\d+$", "", chrom_counts$accession)
tree_overview_accession_ggtree_rev$accession_short <- sub("\\.\\d+$", "", tree_overview_accession_ggtree_rev$accession)

chrom_counts_summary <- aggregate(chromosome ~ accession_short, data = chrom_counts, FUN = length)
colnames(chrom_counts_summary)[2] <- "n_chromosomes"

tree_overview_accession_ggtree_rev$n_chromosomes <- chrom_counts_summary$n_chromosomes[
  match(tree_overview_accession_ggtree_rev$accession_short, chrom_counts_summary$accession_short)
]
tree_overview_accession_ggtree_rev$accession_short <- NULL

# add a second D. melanogaster at the top for better visibility
new_row <- data.frame(
  species = "Drosophila pelanogaster",
  accession = "GCF_000001215.9",
  n_chromosomes = 4
)

tree_overview_accession_ggtree_rev <- rbind(
  new_row,
  tree_overview_accession_ggtree_rev
)

#_____________________________chrom_order_list__________________________________

chrom_order <- list(
  GCF_000001215.9 = c(2,5,3,6,4,1), # 2x D. melanogaster
  GCF_000001215.4 = c(2,5,3,6,4,1),
  GCA_034638295.1 = c(1,3,4,5,2,6),
  GCA_949710015.1 = c(3,1,5,4,2,6),
  GCA_958296145.1 = c(2,3,5,4,1,6),
  GCA_040938175.1 = c(2,4,5,3,1,6),
  GCA_949748255.1 = c(3,4,5,2,1,6),
  GCA_930367235.1 = c(2,3,5,4,1,6),
  GCA_932526305.1 = c(5,4,3,2,1,6),
  GCA_943735925.1 = c(3,2,5,4,1,6),
  GCA_941918925.1 = c(2,4,5,3,1,6),
  GCA_916048285.2 = c(4,2,5,3,1,6),
  GCA_933228675.1 = c(4,2,5,3,1,6), 
  GCA_963854835.1 = c(3,5,1,4,2,6),
  GCA_964194425.1 = c(5,4,3,1,2,6),
  GCF_023373825.1 = c(5,1,3,4,2,6),
  GCA_949127995.1 = c(5,3,1,4,2,6),
  GCA_963457695.1 = c(5,2,3,4,1,6),
  GCA_964263875.1 = c(4,5,3,1,2,6),
  GCA_943735915.1 = c(4,1,3,5,2,6),
  GCA_949752815.1 = c(5,4,1,3,2,6),
  GCA_951804985.1 = c(1,2,4,5,3,6),
  GCA_947389925.1 = c(5,1,3,4,2,6),
  GCA_964659535.1 = c(3,5,2,4,1,6),
  GCA_947095585.1 = c(3,4,2,5,1,6),
  GCA_963932195.1 = c(4,5,2,3,1,6),
  GCA_963855945.1 = c(4,5,1,3,2,6),
  GCA_964036015.1 = c(4,5,2,3,1,6)
)

genomes_order <- tree_overview_accession_ggtree_rev$accession
tsv_files_busco <- paste0(genomes_order, "_reduced.tsv")

#_______________________________________________________________________________
# generate hard-coded R list
#chrom_order_hardcoded <- lapply(
  #seq_len(nrow(tree_overview_accession_ggtree_rev)), 
  #function(i) {
    #seq_len(tree_overview_accession_ggtree_rev$n_chromosomes[i])
  #}
#)
#names(chrom_order_hardcoded) <- tree_overview_accession_ggtree_rev$accession
#cat("chrom_order <- list(\n")
#for(acc in names(chrom_order_hardcoded)){
  #cat("  ", acc, " = c(", paste(chrom_order_hardcoded[[acc]], collapse = ","), "),\n", sep="")
#}
#cat(")\n")

#____________________________load_ALG_data______________________________________

file_path_ALG <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/ALGs_syngraph_brachycera.tsv"
palette_ALG <- c("db1a" = "#25d8a0ff", "db1b" = "#168a65ff",
                 "db2a" = "#fccf8f", "db2b" = "#f29717ff", 
                 "db3a" = "#8d96e5", "db3b" = "#005990ff",
                 "db4a" = "#f0e354ff", "db4b" = "#e2d119ff",
                 "db5" =  "#60b5e1ff", "db6" = "black")

ALGs_syngraph_diptera <- read.table(
  file_path_ALG,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)
colnames(ALGs_syngraph_diptera) <- c("busco", "ALG")

#_______________________execute_functions_______________________________________

genomic_data <- prepare_genome_gene_objects(
  file_path_chrom_names = file_path_chrom_names,
  dir_path_busco = dir_path_busco,
  tsv_files_busco = tsv_files_busco,
  chrom_order = chrom_order
)

my_genomes     <- genomic_data$my_genomes
my_genes       <- genomic_data$my_genes
my_gene_names  <- genomic_data$my_gene_names

#___________________________flip_chromosomes____________________________________

file_path_chrom_names <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/chromosome_names_max12_long_filtered_renamed.tsv"
chrom_names <- read.table(
  file_path_chrom_names,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

chroms_to_flip  <- list(
  GCA_034638295.1 = c(1,2),
  GCF_000001215.4 = c(),
  GCF_000001215.9 = c(),
  GCA_933228675.1 = c(1,4,5), # cen S. lunata
  GCA_916048285.2 = c(2,4), # cen B. pandia
  GCA_941918925.1 = c(1,5,3), # cen P. melanocephala
  GCA_943735925.1 = c(1,3,5), # cen P. amentaria
  GCA_932526305.1 = c(1,5),
  GCA_930367235.1 = c(1,3,2), # cen S. rosellei
  GCA_949748255.1 = c(2,4),
  GCA_040938175.1 = c(1,5,2),
  GCA_958296145.1 = c(2,3,4,5),
  GCA_949710015.1 = c(1,2,3),
  GCA_963854835.1 = c(1,3),
  GCA_964194425.1 = c(1,2),
  GCF_023373825.1 = c(3,2,4), # cen B. dorsalis
  GCA_949127995.1 = c(1,2,3,5),
  GCA_963457695.1 = c(1,2,4),
  GCA_964263875.1 = c(5,3,4),
  GCA_943735915.1 = c(2), # cen P. coryleti
  GCA_949752815.1 = c(1,5),
  GCA_951804985.1 = c(1,2,5),
  GCA_947389925.1 = c(1,3,5),
  GCA_964659535.1 = c(2,3),
  GCA_947095585.1 = c(2,3,4,5),
  GCA_963932195.1 = c(1),
  GCA_963855945.1 = c(1,2,4),
  GCA_964036015.1 = c(1,2,5)
)

for (genome_id in names(chroms_to_flip)) {
  
  # Get indices of chromosomes to flip for this genome
  chr_indices <- chroms_to_flip[[genome_id]]
  
  # Get chromosome names for this genome
  genome_chr_names <- names(my_genes[[genome_id]])
  
  for (chr_idx in chr_indices) {
    target_chr <- genome_chr_names[chr_idx]
    
    # Extract chromosome length from your chrom_names file
    header_row <- which(chrom_names[,1] == paste0("chromosome_", genome_id))
    chr_data <- chrom_names[(header_row + 1):nrow(chrom_names), ]
    end_row <- which(chr_data[,1] == "" | is.na(chr_data[,1]))[1]
    chr_data <- chr_data[1:(end_row - 1), ]
    
    chrom_len_individual <- data.frame(
      genome = genome_id,
      chromosome = chr_data[chr_data[, 1] == target_chr, 1],
      length = as.numeric(chr_data[chr_data[, 1] == target_chr, 2]),
      stringsAsFactors = FALSE
    )
    
    L <- chrom_len_individual$length
    pos <- my_genes[[genome_id]][[target_chr]]
    
    # Flip positions
    my_genes[[genome_id]][[target_chr]] <- L - pos + 1
  }
}

#______________________________plot_genomes_____________________________________

plot_genomes_busco(
  genomes = my_genomes,
  genes = my_genes,
  gene_names = my_gene_names,
  ALG_table = ALGs_syngraph_diptera,
  palette_ALG = palette_ALG,
  spacing = 17000000,
  chrom_height = 1,
  filename = "/Users/jg40/Desktop/genomes_conserved.png",
  width = 6000,
  height = 4000,
  genomes_order = genomes_order,
  centromeres = centromeres
)

#_______________________________________________________________________________

