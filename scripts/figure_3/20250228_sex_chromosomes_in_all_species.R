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

###################################################


########### PLOTTING SEX CHR ###############

source('scripts/20250620_colour_pal.R')

number_of_species <- nrow(sorted_data_tips_desc)
species_step <- (1 / number_of_species)
species_y <- (1:number_of_species / number_of_species) - (species_step / 2)
names(species_y) <- rev(sorted_data_tips_desc$label)
between_chr_gap_size <- 0.1

pdf('figures/sex_chrom_base.pdf', height = 26, width = 3)

plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = F, xlab = '', ylab = '')

for ( i in 1:nrow(sorted_data_tips_desc)){
  species <- as.character(sorted_data_tips_desc[i, 'label'])
  sp_tab <- all_genome_data[all_genome_data[, 'species'] == species, ]
  
  if(any(sp_tab[, 'X'])){
    chrom_to_plot <- sp_tab[sp_tab[, 'X'], ]
    chr_number <- nrow(chrom_to_plot)
    total_ch_space <- 1 - ((chr_number - 1) * between_chr_gap_size)
    per_ch_space <- total_ch_space / chr_number
    y_bot <- species_y[species] - (species_step / 2)
    y_top <- species_y[species] + (species_step / 2)
    for (sex_ch in 1:chr_number){
      bar_subset <- unlist(chrom_to_plot[sex_ch, paste0('d', 1:6)])
      # making sure the order in the plot will be consisent

      columns_to_plot <- c(0, cumsum(bar_subset)) / sum(bar_subset)

      for(i in 1:length(bar_subset)){
          rect(columns_to_plot[i], y_bot,
               columns_to_plot[i + 1], y_top, 
              col = pal[i], bty = 'n')
      }
    }
  }

}

dev.off()

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