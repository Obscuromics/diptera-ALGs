# Figure 2C - color given genomes with ALGs
################################################################################
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyr)
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
################################################################################
#root <- "/Users/ab66/Documents/sanger_work/diptera/diptera-ALGs/"
root <- paste0(getwd(), "/")
################################################################################
# color palette
pal <- c("M1" = "#1573afff", "M2" = "#e59d38ff", "M3" = "#f0e354ff", 
         "M4" = "#169e73ff", "M5" = "#60b5e1ff", "M6" = "black", "100" = "white")
################################################################################
n3 <- read.table(paste0(root, "data/minus26_89_n3.tsv"), header = TRUE, sep = "\t")
ALG_6 <- read.table(paste0(root, "data/small_dot_ALG.txt"), header = TRUE, sep = "\t")
ALG_data <- rbind(n3, ALG_6)
colnames(ALG_data) <- c("busco", "chrQ", "Qstart", "Qend")

names_acc_fam_2 <- read.table(paste0(root, "data/names_acc_fam_2.csv"), sep = ",", header = TRUE)
colnames(names_acc_fam_2) <- c("names", "num", "acc", "fam")

directory <- paste0(root, "data/BUSCO_data/")
new_colnames <- c("busco", "chrQ", "Qstart", "Qend")
################################################################################
# Drosophilidae
# Empididae
# Bibionidae
# Ptychopteridae
# Chironomidae
# Culicinae

keep_small_d = names_acc_fam_2
keep_small_d$busco_files <- paste0(keep_small_d$acc, ".tsv")
keep_small_d$chrom_files <- paste0(keep_small_d$acc, "_chrominf.tsv")

spec_keep_dot <- c("Drosophila_melanogaster",
                   "Empis_livida",
                   "Dilophus_febrilis",
                   "Bibio_marci",
                   "Ptychoptera_contaminata",
                   "Ptychoptera_albimana",
                   "Chironomus_riparius",
                   "Aedes_albopictus",
                   "Culex_pipiens",
                   "Topomyia_yanbarensis"
                   #"Tipula_confusa",
                   #"Bombylius_discolor",
                   #"Bombylius_major",
                   #"Anthrax_anthrax",
                   #"Villa_cingulata"
                   )

keep_filtered <- keep_small_d[keep_small_d$names %in% spec_keep_dot, ]
keep_filtered <- keep_filtered[match(rev(spec_keep_dot), keep_filtered$names),]

busco_files <- keep_filtered$busco_files
chrom_files <- keep_filtered$chrom_files

# plot
plot_list <- list()

for (j in 1:length(busco_files)) {
  
  busco_file_path <- file.path(directory, busco_files[j])
  chrom_file_path <- file.path(directory, chrom_files[j])
  species <- keep_filtered$names[j]
  family <- keep_filtered$fam[j]
  
  # Load the busco data
  busco_data <- read.table(busco_file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(busco_data) <- new_colnames
  
  # Match ALG identity of BUSCOs to genome of interest
  busco_data_ALG <- merge(busco_data, ALG_data[, c("busco", "chrQ")], by = "busco", all.x = TRUE)
  busco_data_ALG <- busco_data_ALG[!is.na(busco_data_ALG$chrQ.y), ] # Remove unassigned BUSCOs
  
  # Load the chromosome data
  chrom_data <- read.table(chrom_file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  chromosomes <- unique(busco_data_ALG$chrQ.x)
  chromosomes <- sort(chromosomes)
  chromosomes <- chromosomes[order(chrom_data$order)]
  chrom_data <- chrom_data[order(chrom_data$order), ]
  
  # Sort data by order of chromosomes
  chr_order <- chrom_data$chr
  busco_data_ALG_sorted <- busco_data_ALG[order(factor(busco_data_ALG$chrQ.x, levels = chr_order)), ]
  
  # Insert 80 empty BUSCOs between chromosomes
  busco_data_ALG_sorted <- busco_data_ALG %>%
    mutate(chrQ.x = factor(chrQ.x, levels = chr_order)) %>%
    arrange(chrQ.x, Qstart)
  # remove start and end column
  busco_data_ALG_sorted <- busco_data_ALG_sorted %>%
    select(-Qstart, -Qend)
  
  # invert chromosomes
  filtered_chr <- chrom_data$chr[chrom_data$invert == TRUE]
  unique_chrs <- unique(busco_data_ALG_sorted$chrQ.x)
  flipped_data <- lapply(unique_chrs, function(chr) {
    chr_data <- busco_data_ALG_sorted[busco_data_ALG_sorted$chrQ.x == chr, ]
    if (chr %in% filtered_chr) {
      chr_data <- chr_data[nrow(chr_data):1, ]
    }
    return(chr_data)
  })
  final_data <- do.call(rbind, flipped_data)
  busco_data_ALG_sorted <- final_data
  
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
  # remove everything but chrQ.y
  busco_data_ALG_final <- busco_data_ALG_final %>%
    select(-busco, -chrQ.x)
  
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
  df_aggregated$fill_color <- factor(df_aggregated$value, levels = names(pal))
  
  # Generate the plot
  p <- ggplot(df_aggregated, aes(x = factor(bin), y = count, fill = value)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = pal) +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      panel.grid = element_blank()
    ) +
    labs(title = paste0(species, " (", family, ")"))
  
  # Append the plot to the plot_list
  plot_list[[j]] <- p
}

plot_list <- rev(plot_list)

plt_all <- grid.arrange(grobs = plot_list, ncol = 1)
ggsave(paste0(root, "figures/genomes_painted_by_ALGs.png"), plot = plt_all, dpi = 600)
ggsave(paste0(root, "figures/genomes_painted_by_ALGs.svg"), plot = plt_all, dpi = 600)
################################################################################

# Figure 2D - number of BUSCOs per ALG
agg <- aggregate(ALG_data$busco, by = list(ALG_data$chrQ), FUN = length)
colnames(agg) <- c("ALG", "n_markers")
agg$ALG_label <- sub("M", "", agg$ALG)

plt_n_markers <- ggplot(agg, aes(x = ALG_label, y = n_markers, fill = ALG)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal) +
  labs(x = "ALG", y = "# of BUSCOs") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(paste0(root, "figures/n_buscos_per_ALGs.png"), plot = plt_n_markers, dpi = 600, width = 10, height = 6)
ggsave(paste0(root, "figures/n_buscos_per_ALGs.svg"), plot = plt_n_markers, dpi = 600, width = 10, height = 6)
################################################################################