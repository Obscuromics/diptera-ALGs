

#_______________________________________________________________________________
library(ape)
library(zoo)

#________________________________load_data______________________________________

tree_file <- "/Users/jg40/Documents/diptera_centromeres/tree/diptera.supermatrix.phy.treefile.dixa"
tree_diptera <- read.tree(tree_file)

diptera_species <- "/Users/jg40/Documents/diptera_centromeres/tree/20260119_all_species.tsv"
diptera_species <- read.delim(diptera_species, sep = "\t", header = TRUE)

#___________________________overview_species____________________________________

keep_species <- c(
 "Culicoides_brevitarsis",
 "Culicoides_sonorensis",
 "Dixa_nubilipennis",
 "Forcipomyia_palustris",
 "Anopheles_coluzzii",
 "Anopheles_nili"
)

#_____________________________build_tree________________________________________

common_species <- intersect(tree_diptera$tip.label, keep_species)
tree_overview <- keep.tip(tree_diptera, common_species)
overview_tips <- tree_overview$tip.label
tree_depth <- max(node.depth.edgelength(tree_overview))
plot(
  tree_overview,
  cex = 2.5,
  edge.width = 1.9,
  align.tip.label = TRUE,
  no.margin = TRUE,
  x.lim = c(0, tree_depth * 3.0)
)

tree_overview_accession <- data.frame(
  species   = overview_tips,
  accession = diptera_species$accession[
    match(overview_tips, diptera_species$species)
  ],
  stringsAsFactors = FALSE
)

#_____________________generate_raw_data_for_plots_______________________________

file_path_chrom_names <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/chromosome_names_max12_long_filtered_renamed.tsv"
dir_path_busco <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/trial"

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
tree_overview_accession$accession_short <- sub("\\.\\d+$", "", tree_overview_accession$accession)

chrom_counts_summary <- aggregate(chromosome ~ accession_short, data = chrom_counts, FUN = length)
colnames(chrom_counts_summary)[2] <- "n_chromosomes"

tree_overview_accession$n_chromosomes <- chrom_counts_summary$n_chromosomes[
  match(tree_overview_accession$accession_short, chrom_counts_summary$accession_short)
]
tree_overview_accession$accession_short <- NULL
tree_overview_accession <- tree_overview_accession[nrow(tree_overview_accession):1, ]

#________________________chrom_order_list_______________________________________

chrom_order <- list(
  GCF_943734685.1 = c(1,3,2),
  GCF_943737925.1 = c(3,1,2),
  GCA_000.2 = c(4,2,3,1), # D. nubilipennis
  GCF_036172545.1 = c(3,1,2),
  GCA_047716325.1 = c(3,1,2),
  GCA_976917515.1 = c(3,1,2)
)

genomes_order <- tree_overview_accession$accession
tsv_files_busco <- paste0(genomes_order, "_reduced.tsv")

#______________________________load_ALG_data____________________________________

file_path_ALG <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/ALGs_syngraph_diptera.tsv"

palette_ALG <- c("d1" = "#169e73ff", "d2" = "#e59d38ff", "d3" = "#1573afff",
                 "d4" = "#f0e354ff", "d5" = "#60b5e1ff", "d6" = "black")

#palette_ALG <- c("d1" = "#F2F2F2", "d2" = "#F2F2F2", "d3" = "#F2F2F2", 
                 #"d4" = "#F2F2F2", "d5" = "#F2F2F2", "d6" = "#D00000")

ALGs_syngraph_diptera <- read.table(
  file_path_ALG,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)
colnames(ALGs_syngraph_diptera) <- c("busco", "ALG")

#______________________________call_functions___________________________________

genomic_data <- prepare_genome_gene_objects(
  file_path_chrom_names = file_path_chrom_names,
  dir_path_busco = dir_path_busco,
  tsv_files_busco = tsv_files_busco,
  chrom_order = chrom_order
)
my_genomes     <- genomic_data$my_genomes
my_genes       <- genomic_data$my_genes
my_gene_names  <- genomic_data$my_gene_names

#______________________________flip_chromosomes_________________________________

file_path_chrom_names <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/chromosome_names_max12_long_filtered_renamed.tsv"
chrom_names <- read.table(
  file_path_chrom_names,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

chroms_to_flip <- list(
  GCA_976917515.1 = c(3), # F. palustris
  GCA_047716325.1 = c(), # C. sonorensis
  GCF_036172545.1 = c(1), # C. brevitarsis
  GCA_000.2 = c(4), # Dixidae
  GCF_943737925.1 = c(2), # An. nili
  GCF_943734685.1 = c(2) # An. coluzzii
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

#___________________________generate_genome_plot________________________________

plot_genomes_busco(
  genomes = my_genomes,
  genes = my_genes,
  gene_names = my_gene_names,
  ALG_table = ALGs_syngraph_diptera,
  palette_ALG = palette_ALG,
  spacing = 17000000,
  chrom_height = 1,
  filename = "/Users/jg40/Desktop/genomes_ceratopogonidae_trial.png",
  width = 6000, height = 4000,
  genomes_order = genomes_order
)

#_________________________ribbon_plot_X_ALG_6___________________________________

file_path_chrom_names <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/chromosome_names_max12_long_filtered_renamed.tsv"
dir_path_busco <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/trial"

chrom_names <- read.table(
  file_path_chrom_names,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

chrom_order_overview <- list(
  GCA_976917515.1 = c(3,1,2),
  GCA_047716325.1 = c(3,1,2),
  GCF_036172545.1 = c(3,1,2),
  GCA_000.2 = c(4,2,3,1), # Dixidae
  GCF_943737925.1 = c(3,1,2),
  GCF_943734685.1 = c(1,3,2)
)

#________________________

df_genome_chrom_len_3 <- do.call(rbind, lapply(names(chrom_order_overview), function(genome) {
  
  chrom_idx <- 1
  idx <- chrom_order_overview[[genome]][chrom_idx]
  header_row <- which(chrom_names[[1]] == paste0("chromosome_", genome))
  block <- chrom_names[(header_row + 1):nrow(chrom_names), ]
  
  stop_row <- which(
    is.na(block[[1]]) | block[[1]] == "" |
      grepl("^chromosome_", block[[1]])
  )[1]
  if (!is.na(stop_row)) {
    block <- block[1:(stop_row - 1), ]
  }
  
  data.frame(
    genome = genome,
    chromosome = block[idx, 1],
    length = block[idx, 2],
    stringsAsFactors = FALSE
  )
}))

#________________________

# filter each chromosome 
df_chr_list <- lapply(1:nrow(df_genome_chrom_len_3), function(i) {
  genome <- df_genome_chrom_len_3$genome[i]
  chromosome <- df_genome_chrom_len_3$chromosome[i]
  df_chr <- data.frame(
    busco_name = my_gene_names[[genome]][[chromosome]],
    busco_position = my_genes[[genome]][[chromosome]],
    stringsAsFactors = FALSE
  )
  df_chr <- merge(
    df_chr,
    ALGs_syngraph_diptera,
    by.x = "busco_name",
    by.y = "busco",
    all.x = TRUE
  )
  colnames(df_chr)[colnames(df_chr) == "ALG"] <- "ALG_assignment"
  
  df_chr
})
names(df_chr_list) <- paste0(df_genome_chrom_len_3$genome, "_", df_genome_chrom_len_3$chromosome)

#________________________

# Extract the accession (first two parts separated by "_")
df_accessions <- sapply(strsplit(names(df_chr_list), "_"), function(x) paste(x[1:2], collapse = "_"))
accession_to_species <- setNames(
  tree_overview_accession$species,
  tree_overview_accession$accession
)
y_labels <- accession_to_species[df_accessions]

# generate plot
png("ribbon_plot_ceratopogonidae_2.png", width = 2600, height = 1800, res = 300)
par(mar = c(5, 15, 4, 2))

n_chr <- length(df_chr_list)
y_positions <- seq(n_chr, 1)
x_range <- range(sapply(df_chr_list, function(df) df$busco_position), na.rm = TRUE)

plot(
  x_range,
  c(0, n_chr + 1),
  type = "n",
  xlab = "",
  ylab = "",
  yaxt = "n",
  main = ,
  cex.main = 1.5,
  cex.lab = 1.5,
  bty = "n"
)

busco_coords <- list()

for (i in seq_along(df_chr_list)) {
  df <- df_chr_list[[i]]
  y <- y_positions[i]
  
  abline(h = y, col = "black", lty = 2, lwd = 0.5)
  
  for (j in seq_len(nrow(df))) {
    busco <- df$busco_name[j]
    x <- df$busco_position[j]
    ALG <- df$ALG_assignment[j]
    
    if (is.na(ALG)) next
    col_busco <- palette_ALG[ALG]
    points(x, y, pch = 16, col = col_busco, cex = 0.5)
    
    busco_coords[[busco]] <- rbind(
      busco_coords[[busco]],
      data.frame(x = x, y = y, ALG = ALG, stringsAsFactors = FALSE)
    )
  }
}

# Connect same BUSCOs across chromosomes
for (busco in names(busco_coords)) {
  coords <- busco_coords[[busco]]
  
  if (nrow(coords) > 1) {
    for (k in 1:(nrow(coords) - 1)) {
      segments(
        x0 = coords$x[k],  y0 = coords$y[k],
        x1 = coords$x[k+1], y1 = coords$y[k+1],
        col = palette_ALG[coords$ALG[k]],
        lwd = 0.2
      )
    }
  }
}

y_labels_clean <- gsub("_", " ", y_labels)
axis(2, at = y_positions, labels = FALSE)

text(
  x = par("usr")[1] - 0.02 * diff(par("usr")[1:2]),  # position slightly left of axis
  y = y_positions,
  labels = y_labels_clean,
  xpd = TRUE,
  adj = 1,
  cex = 1.1,
  font = 3   # italic
)
dev.off()





















