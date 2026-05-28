#____________________function_to_prepare_genomic_data___________________________

prepare_genome_gene_objects <- function(
    file_path_chrom_names,
    dir_path_busco,
    tsv_files_busco,
    chrom_order = NULL
) {
  
  chromosome_name_len <- read.table(
    file_path_chrom_names,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  load_selected_tsv_files <- function(dir_path, tsv_files) {
    result <- list()
    for (file_name in tsv_files) {
      file_path <- file.path(dir_path, file_name)
      var_name <- make.names(gsub("\\.tsv$", "", file_name))
      df <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      colnames(df) <- c("busco", "chromosome", "position")
      result[[var_name]] <- df
    }
    result
  }
  
  busco_tables <- load_selected_tsv_files(
    dir_path = dir_path_busco,
    tsv_files = tsv_files_busco
  )
  
  genome_start_rows <- grep("^chromosome_", chromosome_name_len[[1]])
  genomes_chromosomes <- list()
  
  for (i in seq_along(genome_start_rows)) {
    start <- genome_start_rows[i] + 1
    end <- if (i < length(genome_start_rows))
      genome_start_rows[i + 1] - 1
    else
      nrow(chromosome_name_len)
    
    df <- chromosome_name_len[start:end, ]
    colnames(df) <- c("chromosome", "length")
    
    genome_id <- gsub("^chromosome_", "", chromosome_name_len[genome_start_rows[i], 1])
    genomes_chromosomes[[genome_id]] <- df
  }
  
  busco_genomes <- sub("_reduced$", "", names(busco_tables))
  keep_genomes <- intersect(names(genomes_chromosomes), busco_genomes)
  genomes_chromosomes <- genomes_chromosomes[keep_genomes]
  
  my_genomes <- lapply(genomes_chromosomes, function(df) {
    lengths <- as.numeric(df$length)
    names(lengths) <- df$chromosome
    lengths
  })
  names(my_genomes) <- names(genomes_chromosomes)
  
  my_genes <- list()
  my_gene_names <- list()
  
  for (var in names(busco_tables)) {
    df <- busco_tables[[var]]
    genome_name <- sub("_reduced$", "", var)
    
    chrom_split <- split(df, df$chromosome)
    
    my_genes[[genome_name]] <- lapply(chrom_split, function(x) {
      as.numeric(x$position)
    })
    
    my_gene_names[[genome_name]] <- lapply(chrom_split, function(x) {
      x$busco
    })
  }
  
  for (genome_name in intersect(names(my_genomes), names(my_genes))) {
    common_chroms <- intersect(
      names(my_genomes[[genome_name]]),
      names(my_genes[[genome_name]])
    )
    
    my_genomes[[genome_name]]    <- my_genomes[[genome_name]][common_chroms]
    my_genes[[genome_name]]      <- my_genes[[genome_name]][common_chroms]
    my_gene_names[[genome_name]] <- my_gene_names[[genome_name]][common_chroms]
  }
  
  if (!is.null(chrom_order)) {
    for (genome in intersect(names(chrom_order), names(my_genomes))) {
      order_idx <- chrom_order[[genome]]
      chrom_names <- names(my_genomes[[genome]])
      
      my_genomes[[genome]]    <- my_genomes[[genome]][order_idx]
      my_genes[[genome]]      <- my_genes[[genome]][order_idx]
      my_gene_names[[genome]] <- my_gene_names[[genome]][order_idx]
      
      names(my_genomes[[genome]])    <- chrom_names[order_idx]
      names(my_genes[[genome]])      <- chrom_names[order_idx]
      names(my_gene_names[[genome]]) <- chrom_names[order_idx]
    }
  }
  
  list(
    my_genomes = my_genomes,
    my_genes = my_genes,
    my_gene_names = my_gene_names
  )
}

#_________________BUSCO-only plotting function__________________________________

plot_genomes_busco <- function(genomes, genes, gene_names,
                               ALG_table, palette_ALG,
                               chrom_height = 1, spacing = 1,
                               filename = NULL,
                               width = 1200, height = 800,
                               genomes_order = NULL) {
  
  if (!is.null(genomes_order)) {
    genomes <- genomes[genomes_order]
    genes <- genes[genomes_order]
    gene_names <- gene_names[genomes_order]
  }
  
  for (g in names(genomes)) {
    genomes[[g]][is.na(genomes[[g]])] <- 1
  }
  
  num_cols <- 4
  max_col_width <- numeric(num_cols)
  dummy_width <- spacing / 2
  
  for (i in 1:num_cols) {
    col_lengths <- sapply(genomes, function(g) if(length(g) >= i) g[i] else 0)
    max_col_width[i] <- max(col_lengths, na.rm = TRUE)
    if (is.na(max_col_width[i]) || max_col_width[i] == 0)
      max_col_width[i] <- dummy_width
  }
  
  plot_chromosome <- function(x_start, row_y, length, gene_positions, gene_buscos) {
    x_end <- x_start + length
    
    rect(x_start, row_y, x_end, row_y + chrom_height,
         col = "white", border = "white")
    
    if (!is.null(gene_positions) && length(gene_positions) > 0) {
      gene_colors <- palette_ALG[
        ALG_table$ALG[match(gene_buscos, ALG_table$busco)]
      ]
      
      segments(
        x0 = x_start + gene_positions,
        y0 = row_y,
        x1 = x_start + gene_positions,
        y1 = row_y + chrom_height,
        col = gene_colors,
        lwd = 0.18
      )
    }
    
    return(x_end)
  }
  
  total_width <- sum(max_col_width) + spacing * (num_cols - 1)
  total_height <- length(genomes) * (chrom_height + 1.5)
  
  if (!is.null(filename))
    png(filename, width = width, height = height, res = 600)
  
  plot(0, 0, type = "n",
       xlim = c(0, total_width),
       ylim = c(0, total_height),
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", bty = "n")
  
  for (g_idx in seq_along(genomes)) {
    
    genome <- genomes[[g_idx]]
    genome_genes <- genes[[g_idx]]
    genome_buscos <- gene_names[[g_idx]]
    
    row_y <- (length(genomes) - g_idx) * (chrom_height + 1.5)
    
    x_start <- 0
    
    for (i in 1:num_cols) {
      
      chrom_length <- if(i <= length(genome)) genome[i] else 0
      gene_positions <- if(i <= length(genome_genes)) genome_genes[[i]] else NULL
      gene_buscos <- if(i <= length(genome_buscos)) genome_buscos[[i]] else NULL
      
      if (chrom_length > 0) {
        plot_chromosome(x_start, row_y, chrom_length,
                        gene_positions, gene_buscos)
      }
      
      x_start <- x_start + max_col_width[i] + spacing
    }
  }
  
  if (!is.null(filename)) dev.off()
  
  invisible(TRUE)
}
