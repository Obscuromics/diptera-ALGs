# detailed illustration of sex chromosome turnover events identified
# with sex chromosomes from NCBI
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
read_buscos_2 <- function(file_name, prefix){
  chr_label <- paste0('chr', prefix)
  df <- read.csv(file_name, sep = '\t', comment.char = '#', header = FALSE)
  colnames(df) <- c('busco', chr_label, paste0(prefix, 'start'), paste0(prefix, 'end'))
  return(df)
}
################################################################################
# X chromosomes identified on NCBI
x_list_2 <- c("CM026972.3",
              "CM015046.2",
              "CM017796.1",
              "CM017781.1",
              "CM024333.1",
              "CM025456.2",
              "CM028102.1",
              "CM028117.1",
              "CM028107.1",
              "CM028123.1",
              "CM032179.1",
              "STHB01000014.1",
              "STHB01000015.1",
              "CM042351.1",
              "CM042371.1",
              "CM055406.1",
              "CM061088.1",
              "CM074802.1",
              "CM081796.1",
              "NC_004354.4",
              "NW_007931106.1",
              "NW_007931105.1",
              "NW_007931107.1",
              "NW_001844851.1",
              "NW_007931121.1",
              "NC_051848.1",
              "NC_046673.1",
              "NC_046674.1",
              "NC_047626.1",
              "NC_046672.1",
              "NC_045954.1",
              "NC_047627.2",
              "NC_046683.1",
              "NC_046608.1",
              "NC_050201.1",
              "NC_050155.1",
              "NC_050737.1",
              "NC_050613.1",
              "NC_053034.1",
              "NC_053021.2",
              "NC_052526.2",
              "NC_052525.2",
              "NC_053519.1",
              "NC_054081.1",
              "NC_057931.1",
              "NC_057932.1",
              "NC_060949.1",
              "NC_084885.1",
              "NC_083459.1",
              "NC_066611.1",
              "NC_071503.1",
              "NC_071672.1",
              "NC_089024.1",
              "NC_069143.1",
              "NC_064876.1",
              "NC_064669.1",
              "NC_064870.1",
              "NC_071290.1",
              "NC_071325.1",
              "NC_064600.1",
              "NC_064873.1",
              "NC_069142.1",
              "NC_080707.1",
              "NC_064597.1",
              "NC_071287.1",
              "NC_071293.1",
              "NC_079150.1",
              "NC_079138.1",
              "NC_088785.1",
              "LR989930.1",
              "LR994575.1",
              "LR999968.1",
              "LR999961.1",
              "HG993129.1",
              "OU026158.1",
              "OU026150.1",
              "OU343119.1",
              "OU343164.1",
              "OU426991.1",
              "OU612061.1",
              "OU612048.1",
              "OU612099.1",
              "OU696534.1",
              "OU696699.1",
              "OU744281.1",
              "OU744321.1",
              "OU744341.1",
              "OU862873.1",
              "OV049927.1",
              "OV050032.1",
              "OV277351.1",
              "OV656869.1",
              "OV743721.1",
              "OV839570.1",
              "OV884006.1",
              "OV884061.1",
              "OV884022.1",
              "OW026363.1",
              "OW026370.1",
              "OW026524.1",
              "OW052039.1",
              "OW052046.1",
              "OW052223.1",
              "OW052186.1",
              "OW121745.1",
              "OW121785.1",
              "OW386171.1",
              "OW387028.1",
              "OW388085.2",
              "OW569398.1",
              "OW569408.1",
              "OW584243.1",
              "OW799238.1",
              "OW818034.2",
              "OX016543.1",
              "OX030884.1",
              "OX030885.1",
              "OX030953.1",
              "OX031011.1",
              "OX101757.1",
              "OX122887.1",
              "OX246757.1",
              "OX276339.1",
              "OX297857.1",
              "OX337249.1",
              "OX346249.1",
              "OX352775.1",
              "OX352766.1",
              "OX371226.1",
              "OX371265.1",
              "OX375760.1",
              "OX376376.1",
              "OX376701.1",
              "OX377615.1",
              #"OX377623.1", # Polietes domitor
              "OX392465.1",
              "OX403620.1",
              "OX411885.1",
              "OX421363.1",
              "OX421847.1",
              "OX421898.1",
              "OX422145.1",
              "OX439144.1",
              "OX439486.1",
              "OX443565.1",
              "OX451220.1",
              "OX451348.1",
              "OX454422.1",
              "OX454529.1",
              "OX454324.1",
              "OX454336.1",
              "OX456987.1",
              "OX457082.1",
              "OX457090.1",
              "OX463759.1",
              "OX465090.1",
              "OX493269.1",
              "OX579681.1",
              "OX596078.1",
              "OX596014.1",
              "OX596017.1",
              "OX596006.1",
              "OX596032.1",
              "OX608055.1",
              "OX637537.1",
              "OX638135.1",
              "OX638314.1",
              "OX638373.1",
              "OX638386.1",
              "OX639998.1",
              "OX645065.1",
              "OX940882.1",
              "OY015315.1",
              "OY101445.1",
              "OY101452.1",
              "OY282645.1",
              "OY282582.1",
              "OY283154.1",
              "OY284473.1",
              "OY288111.1",
              "OY390712.1",
              "OY540807.1",
              "OY720156.1",
              "OY720433.1",
              "OY720446.1",
              "OY720149.1",
              "OY727110.1",
              "OY744478.1",
              "OY757144.1",
              "OY757192.1",
              "OY757240.1",
              "OY770265.1",
              "OY776281.1",
              "OY783203.1",
              "OY804262.1",
              "OY829330.1",
              "OY979699.1",
              "OY987260.1",
              "OY987210.1",
              "OY992534.1",
              "OY992541.1",
              "OZ001351.1",
              "OZ007500.1",
              "OZ007565.1",
              "OZ010650.1",
              "OZ010657.1",
              "OZ012640.1",
              "OZ014535.1",
              "OZ017743.1",
              "OZ018401.1",
              "OZ020121.1",
              "OZ020511.1",
              "OZ020530.1",
              "OZ020605.1",
              "OZ021690.1",
              "OZ022267.1",
              "OZ023276.1",
              "OZ024819.1",
              "OZ035878.1",
              "OZ057401.1",
              "OZ124302.1")
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

# retain species from bottom to top

turnover_spec <-c("Panorpa_germanica",
                  "Bradysia_coprophila",
                  "Machimus_atricapillus",
                  "Microdon_myrmicae",
                  "Thecophora_atra",
                  "Sicus_ferrugineus",
                  "Clusia_tigrina",
                  "Teleopsis_dalmanni",
                  "Ornithomya_chloropus",
                  "Ornithomya_fringillina")

turnover_spec_2 <-c("Anopheles")

turnover_fam <-c("Bombyliidae",
                 "Drosophilidae")

remove_spec <- c("Phortica_okadai",
                 "Hirtodrosophila_cameraria",
                 "Drosophila_bifasciata",
                 "Drosophila_madeirensis",
                 "Drosophila_subobscura")
################################################################################
filtered_data <- small_dot_defin_cp[
  small_dot_defin_cp$label %in% turnover_spec | 
    grepl(paste0("^", paste(turnover_spec_2, collapse="|")), small_dot_defin_cp$label) |
    small_dot_defin_cp$Family %in% turnover_fam, 
]

# Ergebnis anzeigen
head(filtered_data)

filtered_data <- small_dot_defin_cp[
  (small_dot_defin_cp$label %in% turnover_spec | 
     grepl(paste0("^", paste(turnover_spec_2, collapse="|")), small_dot_defin_cp$label) |
     small_dot_defin_cp$Family %in% turnover_fam) & 
    !(small_dot_defin_cp$label %in% remove_spec), 
]

# Ergebnis anzeigen
head(filtered_data)


################################################################################
directory <- "C:/Users/julia/Documents/Kamil/busco_fulltable/all_tol"
busco_files <- filtered_data$busco_files
chrom_files <- filtered_data$chrom_files

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
for (k in 1:57) {
  chr_order_Q <- processed_Q_list[[main_counter + 1]]
  chr_order_Q_list[[length(chr_order_Q_list) + 1]] <- chr_order_Q
  main_counter <- main_counter + 5
}
chr_order_Q_list
chr_order_Q <- processed_Q_list[[283]]
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

for (i in 1:57) {
  percent_row <- rep(0, 6)  # Initialisiere mit Nullen
  names(percent_row) <- c("M1", "M2", "M3", "M4", "M5", "M6")  # Benenne die Spalten
  
  cat("Processing index", i, "\n")
  
  if (any(x_list_2 %in% chr_order_Q_list[[i+1]])) {
    cat("Hit found for index", i, "\n")
    hit_value <- chr_order_Q_list[[i+1]][chr_order_Q_list[[i+1]] %in% x_list_2]
    busco_file_path <- file.path(directory, busco_files[i])
    
    busco_data <- read.table(busco_file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    colnames(busco_data) <- new_colnames
    
    busco_data_ALG <- merge(busco_data, ALG_data[, c("busco", "chrQ")], by = "busco", all.x = TRUE)
    busco_data_ALG <- busco_data_ALG[!is.na(busco_data_ALG$chrQ.y), ]
    
    sub_busco_data_ALG <- busco_data_ALG[busco_data_ALG$chrQ.x %in% hit_value, ]
    
    freq_table <- table(sub_busco_data_ALG$chrQ.y)
    print(freq_table)
    percent_table <- prop.table(freq_table) * 100
    
    # Update `percent_row` mit Werten aus `percent_table`
    percent_row[names(percent_table)] <- percent_table
    
  } else {
    cat("No hit for index", i, "\n")
    percent_row <- rep(0, 6)  # Wenn kein Treffer, alle Werte auf 0 setzen
  }
  
  # Füge die Zeile zu `result_df` hinzu
  result_df <- rbind(result_df, as.list(percent_row))
}
colnames(result_df) <- c("M1", "M2", "M3", "M4", "M5", "M6")
head(result_df, 40)
################################################################################
col_list <- c("M1" = "yellow4", "M2" = "springgreen3", "M3" = "#F3cdC7", 
              "M4" = "steelblue2", "M5" = "#FF7B9C", "M6" = "black", "M7" = "white")
col_list <- c("M1" = "#1573afff", "M2" = "#e59d38ff", "M3" = "#f0e354ff", 
              "M4" = "#169e73ff", "M5" = "#60b5e1ff", "M6" = "black", "M7" = "white")
result_df$M7 <- ifelse(rowSums(result_df == 0) == 6, 100, 0)

long_df <- result_df %>%
  mutate(row_id = 1:n()) %>%
  pivot_longer(cols = starts_with("M"), names_to = "category", values_to = "percentage") %>%
  filter(percentage > 0)

long_df$label <- NA
for(i in 1:nrow(long_df)) {
  long_df$label[i] <- filtered_data$label[long_df$row_id[i]]
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
ggsave("barcode_sex_chrom_zoomed.png", plot = p, width = 16, height = 8, dpi = 600)
################################################################################