# make synteny ribbon plot for internal nodes that lead
# from ALGs to Muller elements
################################################################################
library(dplyr)
library(ggplot2)
library(gridExtra)
library(devtools)
################################################################################
root <- getwd()
################################################################################

## To do

## Match nodes between main inference and alg6 inference (DONE)
## Decompose alg6 reconstruction tables to get node tables
## Create new assgn files that contain alg6 markers and give them a chr name
## Also add this chr name to chr info files

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

pal <- c("d1" = "#169e73ff", "d2" = "#e59d38ff", "d3" = "#1573afff",
         "d4" = "#f0e354ff", "d5" = "#60b5e1ff", "d6" = "black", "100" = "white")

# these are the internal nodes, where rearrangement events from ALGs 
# to Muller elements happen
# ALG6 reconstruction corresponding nodes in comments
nodes_to_plot <- c(
  "n13", #n5
  "n21", #n9
  "n36", 
  "n59"  #n18
)

## This is the ALGs + node assignment files in the drive
## Make sure ALG file is actually a tsv - one version wasnt
busco_list <- c("tables/ALGs_syngraph_diptera.tsv", 
  paste0("data/synteny_plot_tables/",nodes_to_plot, "_asgn.tsv"))

## These are tables for the synteny plotter chrom order that look like
#chr order invert  annot
#n59_1 1 F 1
#n59_2 2 F 2
#n59_3 3 F 3
#n59_4 4 F 4
# for the relevant nodes
chrom_list <- c("data/synteny_plot_tables/n1_n2_chrominf.tsv",
                 paste0("data/synteny_plot_tables/",nodes_to_plot, "_chrominf.tsv"))

# This is also in the syngraph drive dir, its unclear to me how this differs from the ALG assignment file
# read ALGs
buscos_to_alg <- read.table(
  file.path("data/syngraph/diptera.no_plecia.mindist.m165_n1_n2.tsv"), 
  header = F, col.names = c('busco', 'ALG'))

# add color code
col <- data.frame("ALG" = names(pal), "colour" = pal)
busco2colour <- left_join(buscos_to_alg, col, by = "ALG")
minimum_buscos = 1

# load synteny_plotter functions
devtools::source_url("https://github.com/Obscuromics/synteny_plotter/blob/dev/scripts/helper_functions.R?raw=TRUE")

# initiate reference
ref_df <- read_buscos_2(busco_list[1], 'R', buscos_to_alg)
ref_chroms <- read.table(chrom_list[1], sep = '\t', header = TRUE) %>% 
  arrange(order)
agg <- ref_df %>% group_by(chrR) %>% summarise(length = max(Rend))
ref_chroms <- left_join(ref_chroms, agg, by = c("chr" = "chrR"))

# crhomosome offset
chr_offset <- 10000000

### generate alignments ###
processed_Q_list <- list()
max_ends <- list()
plot_size <- list()

for (file in busco_list[-1]){
  print(file)
  i <- match(file, busco_list)
  
  # read ref files
  ref_df <- read_buscos_2(busco_list[i-1], 'R', buscos_to_alg)
  ref_chroms <- read.table(chrom_list[i-1], sep = '\t', header = TRUE) %>%
    arrange(order)
  
  # get chromosome sizes for the ref
  agg <- ref_df %>% group_by(chrR) %>% summarise(length = max(Rend))
  ref_chroms <- left_join(ref_chroms, agg, by = c("chr" = "chrR")) %>%
    mutate(Rend = cumsum(as.numeric(length)))
  ref_chroms$Rstart <- c(0, ref_chroms$Rend[1:nrow(ref_chroms) - 1])
  
  # read query file
  query_df <- read_buscos_2(busco_list[i], 'Q', buscos_to_alg)
  query_chroms <- read.table(chrom_list[i], sep = '\t', header = TRUE) %>%
    arrange(order)
  
  # get chromosome sizes for query
  agg <- query_df %>% group_by(chrQ) %>% summarise(length = max(Qend))
  query_chroms <- left_join(query_chroms, agg, by = c("chr" = "chrQ")) %>%
    mutate(Qend = cumsum(as.numeric(length)))
  query_chroms$Qstart <- c(0, query_chroms$Qend[1:nrow(query_chroms) - 1])
  
  # generate alignments
  alignments <- make_alignment_table_upd(
    ref_df, ref_chroms, query_df, query_chroms, chr_offset)
  processed_Q_list <- append(processed_Q_list, 
                             list(alignments, ref_chroms, query_chroms))
  
  max_ends <- append(max_ends, max(ref_chroms$Rend+chr_offset))
  max_ends <- append(max_ends, max(query_chroms$Qend+chr_offset))
  
  plot_size <- append(plot_size, (sum(ref_chroms$length) + chr_offset*(nrow(ref_chroms)-1)))
  plot_size <- append(plot_size, (sum(query_chroms$length) + chr_offset*(nrow(query_chroms)-1)))
}

################################################################################
# plotting
max_end <- max(unlist(max_ends))
plot_length <- max(unlist(plot_size))
gap <- 2
alpha = 0.6
show_outline = TRUE

pdf(file.path(root, "figures/ALG_to_Muller_synteny_plot.pdf"))
print('[+] Generating plot')
plot(0,cex = 0, xlim = c(1, plot_length*1.1), 
     ylim = c(((gap+1)*-1*2*2*2),((gap+1)*2*2*2)),
     xlab = "", ylab = "", bty = "n", yaxt="n", xaxt="n")

main_counter <- 1
y_offset <- -10
y_increment <- 5

while(main_counter <= length(processed_Q_list)){
  
  alignments <- processed_Q_list[[main_counter]]
  ref_chroms <- processed_Q_list[[main_counter+1]]
  query_chroms <- processed_Q_list[[main_counter+2]]
  
  # calculate R and Q adjustment lengths
  if (max(query_chroms$Qend) != max_end){ # i.e. if this is the longest chr_set
    if (max(ref_chroms$Rend) != max_end){
      adjustment_length_R <- (max_end - max(ref_chroms$Rend)) / 2 
      adjustment_length_Q <- (max_end - max(query_chroms$Qend)) / 2 
    }else{
      adjustment_length_R <- 0
      adjustment_length_Q <- (max_end - max(query_chroms$Qend)) / 2 
    }
  }else{
    adjustment_length_Q <- 0
    adjustment_length_R <- (max_end - max(ref_chroms$Rend)) / 2 
  }
  
  if(main_counter != (length(processed_Q_list) - 2)){
    
    # plot alignments
    for (i in ref_chroms$chr){
      if(i %in% alignments$chrR){
        temp <- alignments[alignments$chrR == i,]
        y1 <- gap-y_offset-y_increment
        y2 <- gap-y_offset
        plot_one_ref_chr(temp, adjustment_length_R, adjustment_length_Q,
                         y1, y2, busco2colour, alpha)
      }
    }
    
    # plotting reference chromosomes
    counter <- 0
    offset <- 0
    for (i in ref_chroms$chr){
      Rfirst <- ref_chroms[ref_chroms$chr == i,]$Rstart
      Rlast <- ref_chroms[ref_chroms$chr == i,]$Rend
      
      if(counter != 0){
        Rfirst <- Rfirst + offset
        Rlast <- Rlast + offset
      }
      
      segments(Rfirst+adjustment_length_R, gap-y_offset, 
               Rlast+adjustment_length_R, gap-y_offset, lwd = 10)
      
      #text(x = ((Rlast+Rfirst+1)/2)+adjustment_length_R, y = gap-y_offset, 
       #    label = ref_chroms[ref_chroms$chr == i,]$annot,
        #   srt = 0, cex = 0.5, col = "grey")
      
      offset <- offset + chr_offset # accumulative offset
      counter <- counter + 1
    }
  } else {
    
    # plot alignments
    for (i in ref_chroms$chr){
      temp <- alignments[alignments$chrR == i,]
      y1 <- gap-y_offset-y_increment
      y2 <- gap-y_offset
      plot_one_ref_chr(temp, adjustment_length_R, adjustment_length_Q, 
                       y1, y2, busco2colour, alpha)
    }
    
    # plotting reference chromosomes
    counter <- 0
    offset <- 0
    for (i in ref_chroms$chr){
      Rfirst <- ref_chroms[ref_chroms$chr == i,]$Rstart
      Rlast <- ref_chroms[ref_chroms$chr == i,]$Rend
      
      if(counter != 0){
        Rfirst <- Rfirst + offset
        Rlast <- Rlast + offset
      }
      
      segments(Rfirst+adjustment_length_R, gap-y_offset, 
               Rlast+adjustment_length_R, gap-y_offset, lwd = 10)
      
      #text(x = ((Rlast+Rfirst+1)/2)+adjustment_length_R, y = gap-y_offset, 
       #    label = ref_chroms[ref_chroms$chr == i,]$annot,
        #   srt = 0, cex = 0.5, col = "grey")
      
      offset <- offset + chr_offset # accumulative offset
      counter <- counter + 1
    }
    
    # plotting query chromosomes
    counter <- 0
    offset <- 0
    for (i in query_chroms$chr){
      Qfirst <- query_chroms[query_chroms$chr == i,]$Qstart
      Qlast <- query_chroms[query_chroms$chr == i,]$Qend
      
      if(counter != 0){
        Qfirst <- Qfirst + offset
        Qlast <- Qlast + offset
      }
      
      segments(Qfirst+adjustment_length_Q, gap-y_offset-y_increment, 
               Qlast+adjustment_length_Q, gap-y_offset-y_increment, lwd = 10)
      
      #text(x = ((Qlast+Qfirst+1)/2)+adjustment_length_Q, y = gap-y_offset-y_increment, 
       #    label = query_chroms[query_chroms$chr == i,]$annot,
        #   srt = 0, cex = 0.5, col = "grey")
      
      offset <- offset + chr_offset # accumulative offset
      counter <- counter + 1
    }
  }
  
  main_counter <- main_counter + 3
  y_offset <- y_offset + y_increment
}

dev.off()
################################################################################