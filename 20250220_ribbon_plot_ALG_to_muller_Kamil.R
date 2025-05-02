# make synteny ribbon plot for internal nodes that lead
# from ALGs to Muller elements
################################################################################
library(dplyr)
library(ggplot2)
library(gridExtra)
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
read_buscos_2 <- function(file_name, prefix){
  chr_label <- paste0('chr', prefix)
  df <- read.csv(file_name, sep = '\t', comment.char = '#', header = FALSE)
  colnames(df) <- c('busco', chr_label, paste0(prefix, 'start'), paste0(prefix, 'end'))
  return(df)
}
################################################################################
# initiate reference
ref_df <- ref_df <- read_buscos_2("C:/Users/julia/Documents/Kamil/busco_fulltable/all_tol/ALG_to_muller/n3.tsv", 'R')
ref_chroms <- read.table("C:/Users/julia/Documents/Kamil/busco_fulltable/all_tol/ALG_to_muller/n3_chrominf.tsv", sep = '\t', header = TRUE)
ref_chroms <- ref_chroms %>% arrange(order)

# these are the 5 internal nodes, where rearrangement events from ALGs
# to Muller elements happen
busco_files <- c("n13_final.tsv",
                 "n21_final.tsv",
                 "n35_final.tsv",
                 "n51_final.tsv",
                 "n101_final.tsv")

chrom_files <- c("n13_chrominf.tsv",
                 "n21_chrominf.tsv",
                 "n35_chrominf.tsv",
                 "n51_chrominf.tsv",
                 "n101_chrominf.tsv")
################################################################################
directory <- 'C:/Users/julia/Documents/Kamil/busco_fulltable/all_tol/ALG_to_Muller'
#col_list <- c( "steelblue2","springgreen3", "yellow4", "#F3D8C7", "#FF7B9C")
col_list <- c("#169e73ff", "#e59d38ff", "#1573afff", "#f0e354ff", "#60b5e1ff")
################################################################################
# part of ribbon plot code from Charlotte
chrom_to_col <- list()
counter <- 1
for (chr in ref_chroms$chr){
  chrom_to_col[chr] <- col_list[[counter]]
  counter <- counter + 1
}

## transfer colors from chroms to buscos
buscos_to_col <- list()
for (busco_id in ref_df$busco){
  chr_id <- ref_df[ref_df$busco == busco_id,]$chrR
  
  if (length(chr_id) == 1){
    if (chr_id %in% names(chrom_to_col)){
      buscos_to_col[[busco_id]] <- chrom_to_col[[chr_id]]
    }else{
      buscos_to_col[[busco_id]] <- "grey"
    }
  }else{
    buscos_to_col[[busco_id]] <- "grey"
  }
}
################################################################################
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
  colnames(temp_ref_df) <- c('busco', 'chrR', 'Rstart', 'Rend')#, 'Rstrand')
  temp_ref_chroms <- query_chroms
}
################################################################################
# chralottes plot
max_end <- max(unlist(max_ends))
plot_length <- max_end
gap <- 5
show_outline = TRUE

main_counter <- 1
y_increment <- 14
y_offset <- -(4.8*y_increment)

pdf(paste0('trial_2', '.pdf'),width = 15, height = 100)
print('[+] Generating plot')
plot(0,cex = 0, xlim = c(1, plot_length),
     ylim = c(-gap * 150, max(y_offset + y_increment, gap * 10)),
     xlab = "", ylab = "", bty = "n", yaxt="n", xaxt="n")

main_counter <- 1
y_increment <- 14
y_offset <- -(4.8*y_increment)

for (k in busco_list){
  alignments <- processed_Q_list[[main_counter]]
  chr_order_R <- processed_Q_list[[main_counter+1]]
  chr_order_Q <- processed_Q_list[[main_counter+2]]
  offset_list_R <- processed_Q_list[[main_counter+3]]
  offset_list_Q <- processed_Q_list[[main_counter+4]]
  
  if (max(alignments$Qend) != max_end){
    if (max(alignments$Rend) != max_end){
      adjustment_length_R <- (max_end - max(alignments$Rend)) / 2
      adjustment_length_Q <- (max_end - max(alignments$Qend)) / 2
    }else{
      adjustment_length_R <- 0
      adjustment_length_Q <- (max_end - max(alignments$Qend)) / 2
    }
  }else{
    adjustment_length_Q <- 0
    adjustment_length_R <- (max_end - max(alignments$Rend)) / 2
  }
  
  # plot alignments
  counter <- 1
  for (i in chr_order_R){
    temp <- alignments[alignments$chrR == i,]
    plot_one_ref_chr(temp, buscos_to_col, 
                     adjustment_length_R, adjustment_length_Q, y_offset, alpha=0.6)
    counter <- counter + 1
  }
  
  # plot chromosomes:
  counter <- 1  # plot chromosomes for query:
  for (i in chr_order_Q){
    temp <- alignments[alignments$chrQ == i,]
    Qfirst <- min(temp$Qstart)
    Qlast <- max(temp$Qend)
    
    line_color <- ifelse(i %in% x_list_2, "blue", "black") ###############
    
    segments(Qfirst + adjustment_length_Q, 1 - gap - y_offset, 
             Qlast + adjustment_length_Q, 1 - gap - y_offset, 
             lwd = 5, col = line_color)                  ###################
  }
  
  counter <- 1   # plot chromosomes for ref:
  for (i in chr_order_R){
    temp <- alignments[alignments$chrR == i,]
    Rfirst <- min(temp$Rstart)
    Rlast <- max(temp$Rend)
    segments(Rfirst+adjustment_length_R, gap-y_offset, 
             Rlast+adjustment_length_R, gap-y_offset, lwd = 5)
  }
  
  # add text labels:
  if (main_counter == 1){
    counter <- 1 # text labels for ref:
    for (i in chr_order_R){
      temp <- alignments[alignments$chrR == i,]
      Rfirst <- min(temp$Rstart)
      Rlast <- max(temp$Rend)
      offset <- offset_list_R[[counter]]
      text(x = ((Rlast+Rfirst+1)/2)+adjustment_length_R, y = 4.8*14+8, 
           label = i,
           srt = 0, cex = 0.5)
      counter <- counter + 1
    }
  }
  
  counter <- 1  # text labels for query
  for (i in chr_order_Q){
    temp <- alignments[alignments$chrQ == i,]
    Qfirst <- min(temp$Qstart)
    Qlast <- max(temp$Qend)
  }
  main_counter <- main_counter + 5
  y_offset <- y_offset + y_increment
}
dev.off()
################################################################################
# make function to resort the original nXY files
process_busco_file <- function(file_name, n3) {
  
  file_path <- file.path("C:/Users/julia/Documents/Kamil/busco_fulltable/all_tol/ALG_to_muller", file_name)
  df <- read.table(file_path, sep = "\t", header = FALSE)
  df_merged <- merge(df, n3[, c("busco", "chrR")], by.x = "V1", by.y = "busco", all.x = TRUE)
  colnames(df_merged) <- c("busco", "chromosome", "start", "end", "ALG")
  
  df_sorted <- df_merged %>%
    arrange(chromosome, ALG, busco)
  
  df_final <- df_sorted %>%
    group_by(chromosome) %>%
    mutate(start = seq(0, by = 30000, length.out = n()),
           end = start + 30000) %>%
    ungroup() %>%
    select(-ALG)
  
  output_file <- paste0(sub(".tsv", "", file_name), "_final.tsv")
  write.table(df_final, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  return(output_file)
}

#busco_files <- c("n13.tsv", "n21.tsv", "n35.tsv", "n51.tsv", "n101.tsv")
#output_files <- lapply(busco_files, process_busco_file, n3 = n3)
################################################################################
