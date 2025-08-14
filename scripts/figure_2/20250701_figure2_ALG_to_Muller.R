# make synteny ribbon plot for internal nodes that lead
# from ALGs to Muller elements
################################################################################
library(dplyr)
library(ggplot2)
library(gridExtra)
################################################################################
root <- paste0(getwd(), "/")
################################################################################
# read busco files
read_buscos_2 <- function(file_name, prefix, buscos_to_alg){
  chr_label <- paste0('chr', prefix)
  df <- read.csv(file_name, sep = '\t', comment.char = '#', header = FALSE,
                 na.strings = c("", "NA"))
  colnames(df) <- c("busco", "chr")
  df <- left_join(df, buscos_to_alg)
  df <- df %>% na.omit() %>% arrange(chr)
  
  df_out <- NULL
  
  # add start and end coordinates
  for(d in unique(df$chr)){
    d_df <- df[df$chr == d,]
    d_df <- d_df %>% arrange(desc(ALG))
    d_df$start <- seq(from = 0, by = 30000, length.out = nrow(d_df))
    d_df$end <- seq(from = 30000, by = 30000, length.out = nrow(d_df))
    df_out <- rbind(df_out, d_df[,c('busco', 'chr', 'start', 'end')])
  }
  
  colnames(df_out) <- c('busco', chr_label, paste0(prefix, 'start'), paste0(prefix, 'end'))
  return(df_out)
}

################################################################################
# these are the 5 internal nodes, where rearrangement events from ALGs
# to Muller elements happen
busco_files <- c("diptera.no_plecia.mindist.m165_n1_n2.tsv",
                 "n13_asgn.tsv",
                 "n21_asgn.tsv",
                 "n36_asgn.tsv",
                 "n59_asgn.tsv",
                 "n105_asgn.tsv",
                 "n271_asgn.tsv")

chrom_files <- c("diptera.no_plecia.mindist.m165_n1_n2_chrominf.tsv",
                 "n13_chrominf.tsv",
                 "n21_chrominf.tsv",
                 "n36_chrominf.tsv",
                 "n59_chrominf.tsv",
                 "n105_chrominf.tsv",
                 "n271_chrominf.tsv")
################################################################################
directory <- paste0(root, "data/ALG_to_Muller")
pal <- c("d1" = "#169e73ff", "d2" = "#e59d38ff", "d3" = "#1573afff",
         "d4" = "#f0e354ff", "d5" = "#60b5e1ff", "d6" = "black", "100" = "white")
################################################################################
# read ALGs
buscos_to_alg <- read.table(paste0(root, "data/diptera.no_plecia.mindist.m165_n1_n2.tsv"), header = TRUE, sep = "\t")[1:2]
colnames(buscos_to_alg) <- c("busco", "ALG")

# add color code
col <- data.frame("ALG" = names(pal), "colour" = pal)
buscos_to_col <- left_join(buscos_to_alg, col, by = "ALG")
################################################################################
busco_list <- file.path(directory, busco_files)
chrom_list <- file.path(directory, chrom_files)
minimum_buscos = 1

# load synteny_plotter functions
devtools::source_url("https://github.com/Obscuromics/synteny_plotter/blob/dev/scripts/helper_functions.R?raw=TRUE")

# initiate reference
ref_df <- read_buscos_2(busco_list[1], 'R', buscos_to_alg)
ref_chroms <- read.table(chrom_list[1], sep = '\t', header = TRUE)
ref_chroms <- ref_chroms %>% arrange(order)

# get chromosome sizes
agg <- ref_df %>% group_by(chrR) %>% summarise(length = max(Rend))
ref_chroms <- left_join(ref_chroms, agg, by = c("chr" = "chrR"))

chr_offset <- max(ref_chroms$length) / 2
processed_Q_list <- list()
max_ends <- list()

# save ref as temp_ref
temp_ref_chroms <- ref_chroms
temp_ref_df <- ref_df

for (file in busco_list[-1]){
  i <- match(file, busco_list)
  query_df <- read_buscos_2(file, 'Q', buscos_to_alg)
  query_chroms <- read.table(chrom_list[i], sep = '\t', header = TRUE)
  
  # get sizes
  agg <- query_df %>% group_by(chrQ) %>% summarise(length = max(Qend))
  query_chroms <- left_join(query_chroms, agg, by = c("chr" = "chrQ"))
  
  processed_Q <- make_alignment_table(temp_ref_df, temp_ref_chroms, query_df, query_chroms, chr_offset)
  alignments <- processed_Q$alignments
  processed_Q_list <- append(processed_Q_list, processed_Q)
  max_ends <- append(max_ends, max(alignments$Rend))
  max_ends <- append(max_ends, max(alignments$Qend))
  temp_ref_df <- query_df
  colnames(temp_ref_df) <- c('busco', 'chrR', 'Rstart', 'Rend')
  temp_ref_chroms <- query_chroms
}

################################################################################
max_end <- max(unlist(max_ends))
plot_length <- max_end
gap <- 5
alpha = 0.6
show_outline = TRUE

pdf(paste0(root, "figures/ALG_to_Muller_synteny_plot", '.pdf'))
print('[+] Generating plot')
plot(0,cex = 0, xlim = c(1, plot_length), 
     #ylim = c(((gap+1)*-1*length(busco_list)*2),((gap+1)*length(busco_list)*2)),
     ylim = c(((gap+1)*-1*2*2*2),((gap+1)*2*2*2)),
     #ylim = c(-40, 40),
     xlab = "", ylab = "", bty = "n", yaxt="n", xaxt="n")

main_counter <- 1
y_offset <- -10
y_increment <- 9.5

for (file in busco_list[-1]){
  print(file)
  j <- match(file, busco_list)
  query_chroms <- read.table(chrom_list[j], sep = '\t', header = TRUE)
  
  alignments <- processed_Q_list[[main_counter]]
  chr_order_R <- processed_Q_list[[main_counter+1]]
  chr_order_Q <- processed_Q_list[[main_counter+2]]
  offset_list_R <- processed_Q_list[[main_counter+3]]
  offset_list_Q <- processed_Q_list[[main_counter+4]]
  
  if (max(alignments$Qend) != max_end){ # i.e. if this is the longest chr_set
    if (max(alignments$Rend) != max_end){
      adjustment_length_R <- (max_end - max(alignments$Rend)) / 2 
      adjustment_length_Q <- (max_end - max(alignments$Qend)) / 2 
    }
    else{
      adjustment_length_R <- 0
      adjustment_length_Q <- (max_end - max(alignments$Qend)) / 2 
    }
  }
  else{
    adjustment_length_Q <- 0
    adjustment_length_R <- (max_end - max(alignments$Rend)) / 2 
  }
  
  # plot alignments
  counter <- 1
  for (i in chr_order_R$chr){
    temp <- alignments[alignments$chrR == i,]
    plot_one_ref_chr(temp, adjustment_length_R, adjustment_length_Q, y_offset, buscos_to_col, alpha)
    counter <- counter + 1
  }
  
  # plotting query chromosomes
  if (main_counter == (length(processed_Q_list) - 4)){
    counter <- 0
    offset <- 0
    for (i in chr_order_Q$chr){
      chr_length <- chr_order_Q[chr_order_Q$chr == i,]$length
      Qfirst <- offset
      Qlast <- chr_order_Q[chr_order_Q$chr == i,]$length + offset
      
      if (counter != 0){ # only need to offset start/end if this is not the first chr
        Qfirst <- offset  # allows for accumulative chr positions
        Qlast <- chr_length + offset # allows for accumulative chr positions
      }
      
      offset <- offset + chr_length + chr_offset # accumulative offset
      counter <- counter + 1
      
      segments(Qfirst+adjustment_length_Q, 1-gap-y_offset, 
               Qlast+adjustment_length_Q, 1-gap-y_offset, lwd = 10)
      
      #text(x = ((Qlast+Qfirst+1)/2)+adjustment_length_Q, y = 1-gap-y_offset, 
      #     label = query_chroms[query_chroms$chr == i,]$annot,
      #     srt = 0, cex = 0.5, col = "grey")
      
    }
  }
  
  # plotting reference chromosomes
  counter <- 0
  offset <- 0
  for (i in chr_order_R$chr){
    chr_length <- chr_order_R[chr_order_R$chr == i,]$length
    Rfirst <- offset
    Rlast <- chr_order_R[chr_order_R$chr == i,]$length + offset
    
    if (counter != 0){ # only need to offset start/end if this is not the first chr
      Rfirst <- offset  # allows for accumulative chr positions
      Rlast <- chr_length + offset # allows for accumulative chr positions
    }
    
    offset <- offset + chr_length + chr_offset # accumulative offset
    counter <- counter + 1
    
    segments(Rfirst+adjustment_length_R, gap-y_offset, 
             Rlast+adjustment_length_R, gap-y_offset, lwd = 10)
    
    #text(x = ((Rlast+Rfirst+1)/2)+adjustment_length_R, y = gap-y_offset, 
    #     label = ref_chroms[ref_chroms$chr == i,]$annot,
    #     srt = 0, cex = 0.5, col = "grey")
  }
  
  main_counter <- main_counter + 5
  y_offset <- y_offset + y_increment
  ref_chroms <- query_chroms
}

dev.off()
################################################################################