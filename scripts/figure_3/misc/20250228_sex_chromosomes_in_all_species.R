# sex chromosomes identified by me
# coloured with ALGs
################################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(stringr)
library(ape)
library(ggtree)
library(gsheet)
################################################################################
# read_buscos_2 <- function(file_name, prefix){
#   chr_label <- paste0('chr', prefix)
#   df <- read.csv(file_name, sep = '\t', comment.char = '#', header = FALSE)
#   colnames(df) <- c('busco', chr_label, paste0(prefix, 'start'), paste0(prefix, 'end'))
#   return(df)
# }
################################################################################

####################### WHICH CHROMOSOMES ARE THE X? ###########
# tables/sex_chromosome_turnover_types.tsv
x_chromosome_files <- dir('data/diptera_chromosome_files', pattern = "x.txt", full.names = T)
x_chromosomes <- as.vector(unlist(sapply(x_chromosome_files, read.table)))

all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?gid=1940964825#gid=1940964825", format='csv'),
                            stringsAsFactors = F, header = T, check.names = F)
row.names(all_genome_data) <- all_genome_data[, 'chromosome']

# x_chromosomes[which(!(x_chromosomes %in% all_genome_data[, 'chromosome']))]
# some of the annotated X chromosomes are not in chromosomes!!!
chomosomal_x <- x_chromosomes[x_chromosomes %in% all_genome_data[, 'chromosome']]

print("Number of X chromosomes in assemblies (excl scaffolds):")
print(length(chomosomal_x))

all_genome_data[, 'X'] <- FALSE
all_genome_data[chomosomal_x, 'X'] <- TRUE

################################################################

######################## ORDER OF TIPS IN THE TREE #############
phylo_tree_file <- 'data/diptera.supermatrix.phy.treefile'

tree <- read.tree(phylo_tree_file)
tree_show <- ggtree(tree)
data <- tree_show$data
sorted_data_tips_desc <- data %>% filter(isTip == TRUE) %>% arrange(-y)
################################################################


######################### Chromosomes 2 ALGs ########################
ALG_data <- read.table('data/diptera.no_plecia.mindist.m165_n1_n2.tsv', col.names = c('busco', 'chrQ'))
row.names(ALG_data) <- ALG_data[, 'busco']

all_genome_data[, paste0('d', 1:6)] <- 0

all_busco_data_files <- dir('data/busco_tables', full.names = T)

# head(all_busco_data_files)
# head(all_genome_data[, 'species'])
all(paste0("data/busco_tables/", all_genome_data[, 'species'], ".syngraph.buscos.tsv") %in% all_busco_data_files) 
### ALL FILES ARE WHERE THEY SHOULD BE!!

for (sp in unique(all_genome_data[, 'species'])){
  sp_tab <- read.table(paste0("data/busco_tables/", sp, ".syngraph.buscos.tsv"), col.names = c('busco', 'chr', 'from', 'to'))
  sp_tab[, 'ALG'] <- ALG_data[sp_tab[, 'busco'], 'chrQ']
  sp_tab <- sp_tab[!is.na(sp_tab[, 'ALG']), ]
  ALGs_on_chr <- table(sp_tab$chr, sp_tab$ALG) 
  # subseting to only those that are listed
  ALGs_on_chr <- ALGs_on_chr[row.names(ALGs_on_chr) %in% all_genome_data[, 'chromosome'], ]
  all_genome_data[row.names(ALGs_on_chr), colnames(ALGs_on_chr)] <- ALGs_on_chr
}

write.table(all_genome_data, 'tables/chromosomes_vs_ALGs.tsv', quote = F, sep = '\t', row.names = F, col.names = T)
# all_genome_data <- read.table('tables/chromosomes_vs_ALGs.tsv', sep = '\t', header = T)

####################################################################


########### PLOTTING SEX CHR ###############

source('scripts/20250620_colour_pal.R')

plot_sex_chrom <- function(pdf_file,
                           mode = c("relative", "absolute"),
                           family_gap = 0.0,
                           pdf_height = 26,
                           pdf_width  = 3,
                           multiX_species = character(),
                           show_species_labels = FALSE,
                           show_family_labels  = FALSE,
                           label_cex   = 0.3,
                           label_offset = 0.02,
                           subset_species_file = NULL,
                           subset_species      = NULL,
                           between_chr_gap_size = 0.1,
                           draw_borders = FALSE,
                           border_col   = "black",
                           multiX_layout = c("gapped","columns","stacked"),
                           column_inner_margin = 0.02,
                           repeat_x_axis_in_each_column = TRUE,
                           draw_column_separators = FALSE,
                           stacked_inner_gap = 0.05) {
  mode <- match.arg(mode)
  multiX_layout <- match.arg(multiX_layout)
  
  ## ---------------- CONFIG ----------------
  pal   <- alg_pal
  d_cols <- paste0("d", 1:6)
  
  ## ---------------- INPUT SUBSET (optional) ----------------
  species_to_plot <- sorted_data_tips_desc$label
  if (!is.null(subset_species_file)) {
    subset_from_file <- tryCatch(
      scan(subset_species_file, what = character(), quiet = TRUE),
      error = function(e) character()
    )
    species_to_plot <- intersect(species_to_plot, subset_from_file)
  }
  if (!is.null(subset_species)) {
    species_to_plot <- intersect(species_to_plot, subset_species)
  }
  sdd <- subset(sorted_data_tips_desc, label %in% species_to_plot)
  if (nrow(sdd) == 0) stop("No species to plot after applying subset filters.")
  
  ## ---------------- HELPERS ----------------
  sum_d_cols <- function(df) {
    if (nrow(df) == 0) return(rep(0, length(d_cols)))
    m <- as.matrix(df[, d_cols, drop = FALSE])
    mode(m) <- "numeric"
    m[is.na(m)] <- 0
    colSums(m)
  }
  norm <- function(x) tolower(trimws(x))
  multiX_set <- norm(multiX_species)
  labels_set <- norm(sdd$label)
  not_found  <- setdiff(multiX_set, labels_set)
  if (length(not_found) > 0) {
    message("multiX_species not found in plot labels (will be treated as merged): ",
            paste(not_found, collapse = ", "))
  }
  
  ## ---------------- FAMILY INFO ----------------
  if (family_gap > 0 || show_family_labels) {
    diptera_taxa <- read.csv("data/diptera/diptera_taxa.csv", header = TRUE)
    sdd$family <- diptera_taxa$family[ match(sdd$label, diptera_taxa$species) ]
    sdd$family[is.na(sdd$family)] <- "Unknown"
  } else {
    sdd$family <- NA_character_
  }
  
  ## ---------------- ABSOLUTE MODE PREPASS ----------------
  global_max <- 1
  if (mode == "absolute") {
    bar_totals <- c()
    for (i in 1:nrow(sdd)) {
      species <- as.character(sdd[i, "label"])
      sp_tab <- all_genome_data[all_genome_data[, "species"] == species, , drop = FALSE]
      if (!any(sp_tab[, "X"])) next
      x_rows <- sp_tab[sp_tab[, "X"], , drop = FALSE]
      spec_key <- norm(species)
      if (spec_key %in% multiX_set) {
        if (nrow(x_rows) > 0) {
          m <- as.matrix(x_rows[, d_cols, drop = FALSE]); mode(m) <- "numeric"; m[is.na(m)] <- 0
          bar_totals <- c(bar_totals, rowSums(m))
        }
      } else {
        merged_vals <- sum_d_cols(x_rows)
        bar_totals <- c(bar_totals, sum(merged_vals))
      }
    }
    global_max <- if (length(bar_totals)) max(bar_totals) else 1
  }
  
  ## ---------------- VERTICAL LAYOUT ----------------
  plot_species_order <- rev(sdd$label)
  sdd <- sdd[match(plot_species_order, sdd$label), , drop = FALSE]
  number_of_species <- nrow(sdd)
  
  n_sub_rows <- integer(number_of_species)
  chrom_lists <- vector("list", number_of_species)
  
  for (i in seq_len(number_of_species)) {
    species <- as.character(sdd$label[i])
    sp_tab  <- all_genome_data[all_genome_data[, "species"] == species, , drop = FALSE]
    x_rows  <- if (nrow(sp_tab)) sp_tab[sp_tab[, "X"], , drop = FALSE] else sp_tab[FALSE, , drop = FALSE]
    spec_key <- norm(species)
    
    if (nrow(x_rows) == 0) {
      chrom_lists[[i]] <- x_rows[FALSE, , drop = FALSE]
      n_sub_rows[i] <- 1
    } else if (multiX_layout == "stacked" && (spec_key %in% multiX_set)) {
      chrom_lists[[i]] <- x_rows
      n_sub_rows[i] <- nrow(x_rows)
    } else {
      merged_vals <- sum_d_cols(x_rows)
      one <- x_rows[1, , drop = FALSE]
      one[1, d_cols] <- merged_vals
      chrom_lists[[i]] <- one
      n_sub_rows[i] <- 1
    }
  }
  
  if (multiX_layout != "stacked") {
    if (family_gap > 0) {
      fam_vec_plot <- sdd$family
      family_change_idx <- if (number_of_species > 1) which(fam_vec_plot[-1] != fam_vec_plot[-number_of_species]) else integer(0)
      usable_height <- 1 - length(family_change_idx) * family_gap
      if (usable_height <= 0) stop("family_gap too large for the number of family transitions.")
      species_step <- usable_height / number_of_species
      y_mids <- numeric(number_of_species)
      cursor <- species_step / 2
      for (i in 1:number_of_species) {
        y_mids[i] <- cursor
        cursor <- cursor + species_step
        if (i %in% family_change_idx) cursor <- cursor + family_gap
      }
      species_y <- setNames(y_mids, sdd$label)
    } else {
      species_step <- (1 / number_of_species)
      y_mids <- (1:number_of_species / number_of_species) - (species_step / 2)
      species_y <- setNames(y_mids, sdd$label)
    }
  } else {
    fam_vec_plot <- sdd$family
    family_change_idx <- if (number_of_species > 1) which(fam_vec_plot[-1] != fam_vec_plot[-number_of_species]) else integer(0)
    
    usable_height <- 1 - length(family_change_idx) * family_gap
    if (usable_height <= 0) stop("family_gap too large for the number of family transitions.")
    
    total_units <- sum(n_sub_rows)
    if (total_units == 0) stop("No X chromosomes to plot.")
    unit_h <- usable_height / total_units
    
    bar_h <- unit_h * (1 - stacked_inner_gap)
    gap_h <- unit_h * stacked_inner_gap
    
    species_block_mid <- numeric(number_of_species)
    
    subrow_bounds <- vector("list", number_of_species)
    cursor <- 0
    for (i in 1:number_of_species) {
      k <- n_sub_rows[i]
      if (k == 0) { subrow_bounds[[i]] <- matrix(numeric(0), ncol = 2); next }
      starts <- cursor + (0:(k-1)) * unit_h
      y_bots <- starts + gap_h/2
      y_tops <- y_bots + bar_h
      subrow_bounds[[i]] <- cbind(y_bots, y_tops)
      
      species_block_mid[i] <- cursor + (k * unit_h) / 2
      cursor <- cursor + k * unit_h
      if (i %in% family_change_idx) cursor <- cursor + family_gap
    }
    species_y <- setNames(species_block_mid, sdd$label)
  }
  
  ## ---------------- COLUMNS LAYOUT COUNT (for columns only) ----------------
  n_cols <- 1
  if (multiX_layout == "columns") {
    max_chr <- 1
    for (i in 1:number_of_species) {
      max_chr <- max(max_chr, nrow(chrom_lists[[i]]))
    }
    n_cols <- max_chr
  }
  
  ## ---------------- PLOT ----------------
  xlim <- if (multiX_layout == "columns") c(0, n_cols) else c(0, 1)
  pdf(pdf_file, height = pdf_height, width = pdf_width)
  plot(NULL, xlim = xlim, ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  
  if (multiX_layout == "columns" && draw_column_separators && n_cols > 1) {
    abline(v = 0:n_cols, col = "grey90", lwd = 1)
  }
  
  ## ---------------- DRAW ----------------
  for (i in 1:number_of_species) {
    species <- as.character(sdd$label[i])
    chrom_to_plot <- chrom_lists[[i]]
    
    if (nrow(chrom_to_plot) == 0) next
    
    if (multiX_layout == "gapped" || multiX_layout == "columns") {
      if (multiX_layout == "gapped") {
        chr_number <- nrow(chrom_to_plot)
        y_mid <- species_y[species]
        species_step <- (1 / number_of_species)
      }
      
      if (multiX_layout == "gapped") {
        chr_number <- nrow(chrom_to_plot)
        y_mid <- species_y[species]
        if (!exists("species_step")) {
          species_step <- 1 / number_of_species
        }
        y_bot <- y_mid - (species_step / 2)
        y_top <- y_mid + (species_step / 2)
        
        gap_x <- between_chr_gap_size
        total_bar_width <- max(1 - max(chr_number - 1, 0) * gap_x, 0)
        if (total_bar_width == 0 && chr_number > 0) {
          gap_x <- (1 - 0.02) / max(chr_number - 1, 1)
          total_bar_width <- 1 - (chr_number - 1) * gap_x
        }
        per_bar_width <- if (chr_number > 0) total_bar_width / chr_number else 0
        
        for (k in 1:chr_number) {
          x_left  <- (k - 1) * (per_bar_width + gap_x)
          x_right <- x_left + per_bar_width
          vals <- as.numeric(chrom_to_plot[k, d_cols, drop = TRUE]); vals[is.na(vals)] <- 0
          s <- sum(vals); if (s <= 0) next
          if (mode == "relative") {
            edges01 <- c(0, cumsum(vals) / s)
            edges   <- x_left + (x_right - x_left) * edges01
            for (j in 1:6) rect(edges[j], y_bot, edges[j + 1], y_top, col = pal[j],
                                border = if (draw_borders) border_col else NA)
          } else {
            cum_abs <- c(0, cumsum(vals)) / global_max
            edges_abs <- x_left + pmin(cum_abs, 1) * (x_right - x_left)
            for (j in 1:6) rect(edges_abs[j], y_bot, edges_abs[j + 1], y_top, col = pal[j],
                                border = if (draw_borders) border_col else NA)
          }
        }
      } else if (multiX_layout == "columns") {
        chr_number <- nrow(chrom_to_plot)
        y_mid <- species_y[species]
        if (!exists("species_step")) species_step <- 1 / number_of_species
        y_bot <- y_mid - (species_step / 2)
        y_top <- y_mid + (species_step / 2)
        
        for (k in 1:chr_number) {
          x_left  <- (k - 1) + column_inner_margin
          x_right <-  k      - column_inner_margin
          vals <- as.numeric(chrom_to_plot[k, d_cols, drop = TRUE]); vals[is.na(vals)] <- 0
          s <- sum(vals); if (s <= 0) next
          if (mode == "relative") {
            edges01 <- c(0, cumsum(vals) / s)
            edges   <- x_left + (x_right - x_left) * edges01
            for (j in 1:6) rect(edges[j], y_bot, edges[j + 1], y_top, col = pal[j],
                                border = if (draw_borders) border_col else NA)
          } else {
            cum_abs <- c(0, cumsum(vals)) / global_max
            edges_abs <- x_left + pmin(cum_abs, 1) * (x_right - x_left)
            for (j in 1:6) rect(edges_abs[j], y_bot, edges_abs[j + 1], y_top, col = pal[j],
                                border = if (draw_borders) border_col else NA)
          }
        }
      }
      
    } else if (multiX_layout == "stacked") {
      bounds <- subrow_bounds[[i]]
      for (k in seq_len(nrow(bounds))) {
        y_bot <- bounds[k, 1]; y_top <- bounds[k, 2]
        
        if (k > nrow(chrom_to_plot)) next
        
        vals <- as.numeric(chrom_to_plot[k, d_cols, drop = TRUE]); vals[is.na(vals)] <- 0
        s <- sum(vals); if (s <= 0) next
        
        if (mode == "relative") {
          edges01 <- c(0, cumsum(vals) / s)
          for (j in 1:6) rect(edges01[j], y_bot, edges01[j + 1], y_top, col = pal[j],
                              border = if (draw_borders) border_col else NA)
        } else {
          cum_abs <- c(0, cumsum(vals)) / global_max
          for (j in 1:6) rect(cum_abs[j], y_bot, cum_abs[j + 1], y_top, col = pal[j],
                              border = if (draw_borders) border_col else NA)
        }
      }
    }
    
    ## ---------- OPTIONAL LABELS ----------
    if (show_species_labels) {
      text(x = xlim[1] - label_offset, y = species_y[species], labels = species,
           xpd = TRUE, adj = 1, cex = label_cex)
    }
    if (show_family_labels && !is.na(sdd$family[i])) {
      text(x = xlim[2] + label_offset, y = species_y[species], labels = sdd$family[i],
           xpd = TRUE, adj = 0, cex = label_cex)
    }
  }
  
  ## ---------------- AXIS ----------------
  if (multiX_layout == "columns") {
    usr <- par("usr"); y0 <- usr[3]; tick_len <- strheight("M") * 0.35
    draw_col_axis <- function(x0, x1, at01, labels) {
      segments(x0, y0, x1, y0, xpd = NA)
      xs <- x0 + (x1 - x0) * at01
      segments(xs, y0, xs, y0 - tick_len, xpd = NA)
      text(xs, y0 - tick_len*2.5, labels = labels, xpd = NA, cex = 0.7)
    }
    if (mode == "relative") {
      at01 <- seq(0, 1, by = 0.2); labs <- paste0(seq(0, 100, by = 20), "%")
    } else {
      pretty_breaks <- pretty(c(0, global_max), n = 5)
      at01 <- pretty_breaks / max(global_max, .Machine$double.eps); labs <- pretty_breaks
    }
    for (col in 1:n_cols) {
      x0 <- (col - 1) + column_inner_margin
      x1 <-  col      - column_inner_margin
      draw_col_axis(x0, x1, at01, labs)
    }
    
  } else {
    if (mode == "relative") {
      axis(1, at = seq(0, 1, by = 0.2),
           labels = paste0(seq(0, 100, by = 20), "%"),
           cex.axis = 0.7, line = 0)
    } else {
      pretty_breaks <- pretty(c(0, global_max), n = 5)
      axis(1, at = pretty_breaks / global_max,
           labels = pretty_breaks,
           cex.axis = 0.7, line = 0)
    }
  }
  
  dev.off()
}

# Full dataset, relative
plot_sex_chrom(
  pdf_file = "sex_chrom_base_relative.pdf",
  mode = "relative",
  #family_gap = 0.001,
  #subset_species_file = "tables/sex_chromosome_turnover_types.tsv",
  #pdf_height = 30, pdf_width = 8,
  show_species_labels = FALSE,
  show_family_labels  = FALSE,
  #label_cex = 0.3,
  #between_chr_gap_size = 0.1,
  multiX_species = c("Dolichopus_virgultorum", "Sicus_ferrugineus"),
  draw_borders = TRUE,
  border_col   = "black"
)

# Full dataset, absolute
plot_sex_chrom(
  pdf_file = "sex_chrom_base_absolute.pdf",
  mode = "absolute",
  #family_gap = 0.001,
  #subset_species_file = "tables/sex_chromosome_turnover_types.tsv",
  #pdf_height = 30, pdf_width = 8,
  show_species_labels = FALSE,
  show_family_labels  = FALSE,
  #label_cex = 0.3,
  #between_chr_gap_size = 0.1,
  multiX_species = c("Dolichopus_virgultorum", "Sicus_ferrugineus"),
  draw_borders = TRUE,
  border_col   = "black"
)


# Subset dataset, relative
plot_sex_chrom(
  pdf_file = "sex_chrom_subset_relative.pdf",
  mode = "relative",
  family_gap = 0.01,
  subset_species_file = "tables/sex_chromosome_turnover_types.tsv",
  pdf_height = 8, pdf_width = 8,
  show_species_labels = FALSE,
  show_family_labels  = FALSE,
  #label_cex = 0.3,
  #between_chr_gap_size = 0.1,
  multiX_species = c("Dolichopus_virgultorum", "Sicus_ferrugineus"),
  draw_borders = TRUE,
  border_col   = "black"
)

# Subset dataset, absolute
plot_sex_chrom(
  pdf_file = "sex_chrom_subset_absolute.pdf",
  mode = "absolute",
  family_gap = 0.01,
  subset_species_file = "tables/sex_chromosome_turnover_types.tsv",
  pdf_height = 8, pdf_width = 8,
  show_species_labels = FALSE,
  show_family_labels  = FALSE,
  #label_cex = 0.3,
  #between_chr_gap_size = 0.1,
  multiX_species = c("Dolichopus_virgultorum", "Sicus_ferrugineus"),
  draw_borders = TRUE,
  border_col   = "black"
)

# Full dataset, relative newest version
plot_sex_chrom(
  pdf_file = "sex_chrom_base_relative_v3.pdf",
  mode = "relative",
  #family_gap = 0.001,
  #subset_species_file = "data/diptera/diptera_taxa_selected.txt",
  pdf_height = 26, pdf_width = 3,
  show_species_labels = TRUE,
  show_family_labels  = FALSE,
  #label_cex = 0.3,
  #between_chr_gap_size = 0.1,
  multiX_species = c("Dolichopus_virgultorum", "Sicus_ferrugineus"),
  draw_borders = TRUE,
  border_col   = "black",
  multiX_layout = "stacked"
  #repeat_x_axis_in_each_column = TRUE,
  #stacked_inner_gap = 0.08
)

############################################

# # ALG_data <- ref_df
# # ref_df <- read_buscos_2("C:/Users/julia/Documents/Kamil/busco_fulltable/all_tol/ALG_to_muller/n3.tsv", 'R')
# # ref_chroms <- read.table("C:/Users/julia/Documents/Kamil/busco_fulltable/all_tol/ALG_to_muller/n3_chrominf.tsv", sep = '\t', header = TRUE)
# # ref_chroms <- ref_chroms %>% arrange(order)
# # new_colnames <- c("busco", "chrQ", "Qstart", "Qend")
# # colnames(ALG_data) <- new_colnames
# head(ALG_data)
# ################################################################################
# # X chromosomes identified by me
# # file_path <- "C:/Users/julia/Documents/Kamil/txt/x_list_3.txt"
# # x_list_3 <- read.table(file_path, header = FALSE, sep = "\t")
# # x_list_3 <- as.character(x_list_3[, 1])
# ################################################################################
# # list with 288 genomes in order of the tree
# # most code from 20250222_small_dot_length_for_tree.R
# # small dot status
# file_path <- "C:/Users/julia/Documents/Kamil/txt/small_dot_status.csv"
# small_dot_defin <- read.csv(file_path, header = FALSE, sep = ",")
# small_dot_defin_cp <- small_dot_defin
# colnames(small_dot_defin_cp) <- c("label", "status","position")
# head(small_dot_defin_cp)

# # more info for merging
# names_acc_fam_2 <- read.table("C:/Users/julia/Documents/Kamil/busco_fulltable/names_acc_fam_2.csv", sep = ",", header = TRUE)
# head(names_acc_fam_2)

# small_dot_defin_cp <- small_dot_defin_cp[, "label", drop = FALSE]
# small_dot_defin_cp$Accession <- names_acc_fam_2$Accession[match(small_dot_defin_cp$label, names_acc_fam_2$label)]
# small_dot_defin_cp$Family <- names_acc_fam_2$Family[match(small_dot_defin_cp$label, names_acc_fam_2$label)]
# small_dot_defin_cp$busco_files <- paste0(small_dot_defin_cp$Accession, ".tsv")
# small_dot_defin_cp$chrom_files <- paste0(small_dot_defin_cp$Accession, "_chrominf.tsv")
# head(small_dot_defin_cp)
# ################################################################################
# directory <- "C:/Users/julia/Documents/Kamil/busco_fulltable/all_tol"
# busco_files <- small_dot_defin_cp$busco_files
# chrom_files <- small_dot_defin_cp$chrom_files

# busco_list <- file.path(directory, busco_files)
# chrom_list <- file.path(directory, chrom_files)

# processed_Q_list <- list()
# max_ends <- list()

# temp_ref_chroms <- ref_chroms
# temp_ref_chroms
# temp_ref_df <- ref_df
# temp_ref_df

# for (file in busco_list){
#   i <- match(file, busco_list)
#   query_df <- read_buscos_2(file, 'Q')
#   query_chroms <- read.table(chrom_list[i], sep = '\t', header = TRUE)
#   processed_Q <- make_alignments_table(temp_ref_df, temp_ref_chroms, query_df, query_chroms)
#   alignments <- processed_Q$alignments
#   processed_Q_list <- append(processed_Q_list, processed_Q)
#   max_ends <- append(max_ends, max(alignments$Rend))
#   max_ends <- append(max_ends, max(alignments$Qend))
#   temp_ref_df <- query_df
#   colnames(temp_ref_df) <- c('busco', 'chrR', 'Rstart', 'Rend')
#   temp_ref_chroms <- query_chroms
# }

# # liste mit chromosomen namen erhalten
# chr_order_Q_list <- list()
# main_counter <- 1
# for (k in 1:289) {
#   chr_order_Q <- processed_Q_list[[main_counter + 1]]
#   chr_order_Q_list[[length(chr_order_Q_list) + 1]] <- chr_order_Q
#   main_counter <- main_counter + 5
# }
# chr_order_Q_list
# chr_order_Q <- processed_Q_list[[1443]]
# chr_order_Q_list[[length(chr_order_Q_list) + 1]] <- chr_order_Q
# head(chr_order_Q_list)
# length(chr_order_Q_list)
# ################################################################################
# result_df <- data.frame(
#   M1 = numeric(0),
#   M2 = numeric(0),
#   M3 = numeric(0),
#   M4 = numeric(0),
#   M5 = numeric(0),
#   M6 = numeric(0)
# )

# for (i in 1:289) {
#   percent_row <- rep(0, 6)
#   names(percent_row) <- c("M1", "M2", "M3", "M4", "M5", "M6")
  
#   cat("Processing index", i, "\n")
  
#   # if (any(x_list_3 %in% chr_order_Q_list[[i+1]])) {
#   #   cat("Hit found for index", i, "\n")
#   #   hit_value <- chr_order_Q_list[[i+1]][chr_order_Q_list[[i+1]] %in% x_list_3]
#   #   busco_file_path <- file.path(directory, busco_files[i])
    
#   #   busco_data <- read.table(busco_file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
#   #   colnames(busco_data) <- new_colnames
    
#   #   busco_data_ALG <- merge(busco_data, ALG_data[, c("busco", "chrQ")], by = "busco", all.x = TRUE)
#   #   busco_data_ALG <- busco_data_ALG[!is.na(busco_data_ALG$chrQ.y), ]
    
#   #   sub_busco_data_ALG <- busco_data_ALG[busco_data_ALG$chrQ.x %in% hit_value, ]
    
#   #   freq_table <- table(sub_busco_data_ALG$chrQ.y)
#   #   print(freq_table)
#   #   percent_table <- prop.table(freq_table) * 100
    
#   #   percent_row[names(percent_table)] <- percent_table[names(percent_table)]
    
#   #   missing_cols <- setdiff(names(result_df), names(percent_row))
#   #   percent_row[missing_cols] <- 0
    
#   # } else {
#   #   cat("No hit for index", i, "\n")
#   #   percent_row <- rep(0, 6)
#   #   names(percent_row) <- c("M1", "M2", "M3", "M4", "M5", "M6")
#   # }
  
#   result_df <- rbind(result_df, as.list(percent_row))
# }  
# ################################################################################
# #col_list <- c("M1" = "yellow4", "M2" = "springgreen3", "M3" = "#F3cdC7", 
#               #"M4" = "steelblue2", "M5" = "#FF7B9C", "M6" = "black", "M7" = "white")
# col_list <- c("M1" = "#1573afff", "M2" = "#e59d38ff", "M3" = "#f0e354ff", 
#               "M4" = "#169e73ff", "M5" = "#60b5e1ff", "M6" = "black", "M7" = "white")
# result_df$M7 <- ifelse(rowSums(result_df == 0) == 6, 100, 0)

# long_df <- result_df %>%
#   mutate(row_id = 1:n()) %>%
#   pivot_longer(cols = starts_with("M"), names_to = "category", values_to = "percentage") %>%
#   filter(percentage > 0)

# long_df$label <- NA
# for(i in 1:nrow(long_df)) {
#   long_df$label[i] <- small_dot_defin_cp$label[long_df$row_id[i]]
# }
# print(long_df, n = 50)
# ################################################################################
# p <- ggplot(long_df, aes(x = factor(row_id), y = percentage, fill = category)) +
#   geom_bar(stat = "identity", width = 1, show.legend = TRUE) +
#   scale_fill_manual(values = col_list) +
#   scale_x_discrete(labels = setNames(long_df$label, long_df$row_id)) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(size = 2, angle = 90, hjust = 1),
#     panel.background = element_rect(fill = "white"),
#     plot.background = element_rect(fill = "white"),
#     axis.title = element_blank(),  # Entfernt die Achsentitel
#     legend.position = "none"  # Entfernt die Legende
#   )
# p
# ggsave("barcode_sex_chrom_very_black.png", plot = p, width = 16, height = 2, dpi = 600)
# ################################################################################