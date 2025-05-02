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
################################################################################
read_buscos_2 <- function(file_name, prefix){
  chr_label <- paste0('chr', prefix)
  df <- read.csv(file_name, sep = '\t', comment.char = '#', header = FALSE)
  colnames(df) <- c('busco', chr_label, paste0(prefix, 'start'), paste0(prefix, 'end'))
  return(df)
}
################################################################################
# include the small dot chromosome
n3 <- read.table("C:/Users/julia/Documents/Kamil/muller/minus26_89_n3.tsv", header = TRUE, sep = "\t")
ALG_6 <- read.table("C:/Users/julia/Documents/Kamil/muller/small_dot_ALG.txt", header = TRUE, sep = "\t")
ref_df <- rbind(n3, ALG_6)
ALG_data <- ref_df
ref_df <- ref_df <- read_buscos_2("C:/Users/julia/Documents/Kamil/busco_fulltable/all_tol/ALG_to_muller/n3.tsv", 'R')
ref_chroms <- read.table("C:/Users/julia/Documents/Kamil/busco_fulltable/all_tol/ALG_to_muller/n3_chrominf.tsv", sep = '\t', header = TRUE)
ref_chroms <- ref_chroms %>% arrange(order)
new_colnames <- c("busco", "chrQ", "Qstart", "Qend")
colnames(ALG_data) <- new_colnames
head(ALG_data)
################################################################################
# X chromosomes identified by me
file_path <- "C:/Users/julia/Documents/Kamil/txt/x_list_3.txt"
x_list_3 <- read.table(file_path, header = FALSE, sep = "\t")
x_list_3 <- as.character(x_list_3[, 1])
################################################################################
# list with 288 genomes in order of the tree
# most code from 20250222_small_dot_length_for_tree.R
# small dot status
file_path <- "C:/Users/julia/Documents/Kamil/txt/small_dot_status.csv"
small_dot_defin <- read.csv(file_path, header = FALSE, sep = ",")
small_dot_defin_cp <- small_dot_defin
colnames(small_dot_defin_cp) <- c("label", "status","position")
head(small_dot_defin_cp)

# more info for merging
names_acc_fam_2 <- read.table("C:/Users/julia/Documents/Kamil/busco_fulltable/names_acc_fam_2.csv", sep = ",", header = TRUE)
head(names_acc_fam_2)

small_dot_defin_cp <- small_dot_defin_cp[, "label", drop = FALSE]
small_dot_defin_cp$Accession <- names_acc_fam_2$Accession[match(small_dot_defin_cp$label, names_acc_fam_2$label)]
small_dot_defin_cp$Family <- names_acc_fam_2$Family[match(small_dot_defin_cp$label, names_acc_fam_2$label)]
small_dot_defin_cp$busco_files <- paste0(small_dot_defin_cp$Accession, ".tsv")
small_dot_defin_cp$chrom_files <- paste0(small_dot_defin_cp$Accession, "_chrominf.tsv")
head(small_dot_defin_cp)
################################################################################
directory <- "C:/Users/julia/Documents/Kamil/busco_fulltable/all_tol"
busco_files <- small_dot_defin_cp$busco_files
chrom_files <- small_dot_defin_cp$chrom_files

busco_list <- file.path(directory, busco_files)
chrom_list <- file.path(directory, chrom_files)

processed_Q_list <- list()
max_ends <- list()

temp_ref_chroms <- ref_chroms
temp_ref_chroms
temp_ref_df <- ref_df
temp_ref_df

for (file in busco_list){
  i <- match(file, busco_list)
  query_df <- read_buscos_2(file, 'Q')
  query_chroms <- read.table(chrom_list[i], sep = '\t', header = TRUE)
  processed_Q <- make_alignments_table(temp_ref_df, temp_ref_chroms, query_df, query_chroms)
  alignments <- processed_Q$alignments
  processed_Q_list <- append(processed_Q_list, processed_Q)
  max_ends <- append(max_ends, max(alignments$Rend))
  max_ends <- append(max_ends, max(alignments$Qend))
  temp_ref_df <- query_df
  colnames(temp_ref_df) <- c('busco', 'chrR', 'Rstart', 'Rend')
  temp_ref_chroms <- query_chroms
}

# liste mit chromosomen namen erhalten
chr_order_Q_list <- list()
main_counter <- 1
for (k in 1:289) {
  chr_order_Q <- processed_Q_list[[main_counter + 1]]
  chr_order_Q_list[[length(chr_order_Q_list) + 1]] <- chr_order_Q
  main_counter <- main_counter + 5
}
chr_order_Q_list
chr_order_Q <- processed_Q_list[[1443]]
chr_order_Q_list[[length(chr_order_Q_list) + 1]] <- chr_order_Q
head(chr_order_Q_list)
length(chr_order_Q_list)
################################################################################
result_df <- data.frame(
  M1 = numeric(0),
  M2 = numeric(0),
  M3 = numeric(0),
  M4 = numeric(0),
  M5 = numeric(0),
  M6 = numeric(0)
)

for (i in 1:289) {
  percent_row <- rep(0, 6)
  names(percent_row) <- c("M1", "M2", "M3", "M4", "M5", "M6")
  
  cat("Processing index", i, "\n")
  
  if (any(x_list_3 %in% chr_order_Q_list[[i+1]])) {
    cat("Hit found for index", i, "\n")
    hit_value <- chr_order_Q_list[[i+1]][chr_order_Q_list[[i+1]] %in% x_list_3]
    busco_file_path <- file.path(directory, busco_files[i])
    
    busco_data <- read.table(busco_file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    colnames(busco_data) <- new_colnames
    
    busco_data_ALG <- merge(busco_data, ALG_data[, c("busco", "chrQ")], by = "busco", all.x = TRUE)
    busco_data_ALG <- busco_data_ALG[!is.na(busco_data_ALG$chrQ.y), ]
    
    sub_busco_data_ALG <- busco_data_ALG[busco_data_ALG$chrQ.x %in% hit_value, ]
    
    freq_table <- table(sub_busco_data_ALG$chrQ.y)
    print(freq_table)
    percent_table <- prop.table(freq_table) * 100
    
    percent_row[names(percent_table)] <- percent_table[names(percent_table)]
    
    missing_cols <- setdiff(names(result_df), names(percent_row))
    percent_row[missing_cols] <- 0
    
  } else {
    cat("No hit for index", i, "\n")
    percent_row <- rep(0, 6)
    names(percent_row) <- c("M1", "M2", "M3", "M4", "M5", "M6")
  }
  
  result_df <- rbind(result_df, as.list(percent_row))
}  
################################################################################
#col_list <- c("M1" = "yellow4", "M2" = "springgreen3", "M3" = "#F3cdC7", 
              #"M4" = "steelblue2", "M5" = "#FF7B9C", "M6" = "black", "M7" = "white")
col_list <- c("M1" = "#1573afff", "M2" = "#e59d38ff", "M3" = "#f0e354ff", 
              "M4" = "#169e73ff", "M5" = "#60b5e1ff", "M6" = "black", "M7" = "white")
result_df$M7 <- ifelse(rowSums(result_df == 0) == 6, 100, 0)

long_df <- result_df %>%
  mutate(row_id = 1:n()) %>%
  pivot_longer(cols = starts_with("M"), names_to = "category", values_to = "percentage") %>%
  filter(percentage > 0)

long_df$label <- NA
for(i in 1:nrow(long_df)) {
  long_df$label[i] <- small_dot_defin_cp$label[long_df$row_id[i]]
}
print(long_df, n = 50)
################################################################################
p <- ggplot(long_df, aes(x = factor(row_id), y = percentage, fill = category)) +
  geom_bar(stat = "identity", width = 1, show.legend = TRUE) +
  scale_fill_manual(values = col_list) +
  scale_x_discrete(labels = setNames(long_df$label, long_df$row_id)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 2, angle = 90, hjust = 1),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.title = element_blank(),  # Entfernt die Achsentitel
    legend.position = "none"  # Entfernt die Legende
  )
p
ggsave("barcode_sex_chrom_very_black.png", plot = p, width = 16, height = 2, dpi = 600)
################################################################################