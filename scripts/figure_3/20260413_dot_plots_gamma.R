################################################################################
# Genome synteny dotplot (Genome 1 vs Genome 2)
################################################################################
library(dplyr)
library(tidyverse)
library(ggplot2)
setwd("/Users/jg40/Documents")

################################################################################
# FILE PATHS
################################################################################
# genome_1 rearranged
# genome_2 conserved

# Drosophilidae pair
genome_1 <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/GCF_000001215.4_reduced.tsv" # D. mel
genome_2 <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/GCA_034638295.1_reduced.tsv" # P. okadai

# Tachinidae pair
genome_1 <- "diptera_centromeres/20260119_plot_genomes_chromosomes_genes/GCA_956483585.1_reduced.tsv" # G. viridis
genome_2 <- "diptera_centromeres/20260119_plot_genomes_chromosomes_genes/GCA_932526305.1_reduced.tsv" # Epicampocera succincta

# Syrphidae pair
genome_1 <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/GCA_917880715.2_reduced.tsv" # C. berberina
genome_2 <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/GCA_936431705.1_reduced.tsv" # C. pagana

# Asilidae pair
genome_1 <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/GCA_947538895.1_reduced.tsv" # N. cyanurus
genome_2 <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/GCA_964036015.1_reduced.tsv" # Leptarthrus brevirostris

alg_file <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/ALGs_syngraph_brachycera.tsv"

################################################################################
# CHROMOSOME ORDER (genome 1 & 2)
################################################################################

# Drosophilidae pair
# D.mel
chr_order_genome_1 <- c(
  "NC_004354.4",
  "NT_033779.5",
  "NT_037436.4",
  "NC_004353.4"
)
chr_names_genome_1 <- c(
  "X",
  "2",
  "3",
  "4"
)

# P. okadai
chr_order_genome_2 <- c(
  "CM067848.1",
  "CM067850.1",
  "CM067851.1",
  "CM067852.1",
  "CM067849.1",
  "CM067853.1"
)
chr_names_genome_2 <- c(
  "1",
  "3",
  "4",
  "5",
  "2",
  "6"
)
name_map_genome_1 <- setNames(chr_names_genome_1, chr_order_genome_1)
name_map_genome_2 <- setNames(chr_names_genome_2, chr_order_genome_2)
#_______________________________________________________________________________

# Syrphidae pair
# C. berberina
chr_order_genome_1 <- c(
  "OU862871.1",
  "OU862870.1",
  "OU862872.1",
  "OU862873.1"
)
chr_names_genome_1 <- c(
  "2",
  "1",
  "3",
  "X"
)

# C. pagana
chr_order_genome_2 <- c(
  "OW386168.1",
  "OW386170.1",
  "OW386169.1",
  "OW386167.1",
  "OW386166.1",
  "OW386171.1"
)
chr_names_genome_2 <- c(
  "3",
  "5",
  "4",
  "2",
  "1",
  "X"
)

name_map_genome_1 <- setNames(chr_names_genome_1, chr_order_genome_1)
name_map_genome_2 <- setNames(chr_names_genome_2, chr_order_genome_2)
#_______________________________________________________________________________

# Tachinidae pair
# G. viridis
chr_order_genome_1 <- c(
  "OY101441.1",
  "OY101443.1",
  "OY101440.1",
  "OY101444.1",
  "OY101442.1",
  "OY101445.1"
)
chr_names_genome_1 <- c(
  "2",
  "4",
  "1",
  "5",
  "3",
  "X"
)

# Epicampocera succincta
chr_order_genome_2 <- c(
  "OW052038.1",
  "OW052037.1",
  "OW052036.1",
  "OW052035.1",
  "OW052034.1",
  "OW052039.1"
)
chr_names_genome_2 <- c(
  "5",
  "4",
  "3",
  "2",
  "1",
  "X"
)

name_map_genome_1 <- setNames(chr_names_genome_1, chr_order_genome_1)
name_map_genome_2 <- setNames(chr_names_genome_2, chr_order_genome_2)
#_______________________________________________________________________________

# Asilidae pair
# N. cyanurus
chr_order_genome_1 <- c(
  "OX384543.1",
  "OX384547.1",
  "OX384546.1",
  "OX384544.1",
  "OX384541.1",
  "OX384548.1",
  "OX384540.1",
  "OX384542.1",
  "OX384545.1",
  "OX384549.1"
)
chr_names_genome_1 <- c(
  "4",
  "8",
  "7",
  "5",
  "2",
  "9",
  "1",
  "3",
  "6",
  "10"
)

# Leptarthrus brevirostris
chr_order_genome_2 <- c(
  "OZ038374.1",
  "OZ038375.1",
  "OZ038372.1",
  "OZ038373.1",
  "OZ038371.1",
  "OZ038376.1"
)
chr_names_genome_2 <- c(
  "4",
  "5",
  "2",
  "3",
  "1",
  "6"
)

name_map_genome_1 <- setNames(chr_names_genome_1, chr_order_genome_1)
name_map_genome_2 <- setNames(chr_names_genome_2, chr_order_genome_2)

################################################################################
# ALG COLOUR PALETTE
################################################################################

pal <- c("db1a" = "#25d8a0ff", "db1b" = "#168a65ff",
         "db2a" = "#fccf8f", "db2b" = "#f29717ff", 
         "db3a" = "#8d96e5", "db3b" = "#005990ff",
         "db4a" = "#f0e354ff", "db4b" = "#e2d119ff",
         "db5" =  "#60b5e1ff", "db6" = "black")

################################################################################
# LOAD GENOME DATA
################################################################################

g1 <- read.table(genome_1, header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(g1) <- c("marker","chromosome","start")
g1$end <- g1$start
g1$species <- "genome_1"

g2 <- read.table(genome_2, header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(g2) <- c("marker","chromosome","start")
g2$end <- g2$start
g2$species <- "genome_2"

df <- bind_rows(g1,g2)

################################################################################
# EXCLUDE LONG CHROMOSOME NAMES
################################################################################

df <- df %>% filter(nchar(chromosome) <= 12)

################################################################################
# ADD ORDER COLUMN FOR GENOME 1
################################################################################

df <- df %>%
  mutate(chr_order = case_when(
    species == "genome_1" ~ match(chromosome, chr_order_genome_1),
    species == "genome_2" ~ match(chromosome, chr_order_genome_2),
    TRUE ~ NA_integer_
  ))

################################################################################
# CHROMOSOME SIZES
################################################################################

chr_sizes <- df %>%
  group_by(species, chromosome) %>%
  summarise(chromosome_size_b = max(end), .groups="drop")

df <- left_join(df, chr_sizes, by=c("species","chromosome"))

################################################################################
# LOAD ALG TABLE
################################################################################

alg_df <- read.table(
  alg_file,
  header=FALSE,
  stringsAsFactors=FALSE
)

colnames(alg_df) <- c("marker","ALG")

df <- left_join(df, alg_df, by="marker")
df$ALG[is.na(df$ALG)] <- "unassigned"

################################################################################
# LINEARIZE GENOME COORDINATES
################################################################################

df_lin <- df %>% arrange(species, chr_order, chromosome, start)

df_lin <- df %>%
  arrange(species, chr_order, chromosome, start)

df_lin$linear_start <- NA
df_lin$linear_end <- NA

chr_offset <- 0

for (i in 1:nrow(df_lin)) {
  
  if (i==1 || df_lin$species[i] != df_lin$species[i-1]) {
    chr_offset <- 0
    
  } else if (df_lin$chromosome[i] != df_lin$chromosome[i-1]) {
    chr_offset <- chr_offset + as.integer(df_lin$chromosome_size_b[i-1])
  }
  
  df_lin$linear_start[i] <- chr_offset + df_lin$start[i]
  df_lin$linear_end[i]   <- chr_offset + df_lin$end[i]
}

################################################################################
# REMOVE DUPLICATE MARKERS
################################################################################

duplicates <- df_lin |>
  summarise(n=n(), .by=c(marker,species)) |>
  filter(n>1L)

alg_df2 <- df_lin %>%
  select(marker, ALG) %>%
  distinct()

################################################################################
# WIDE FORMAT FOR DOTPLOT
################################################################################

df_wide <- df_lin %>%
  filter(!marker %in% duplicates$marker) %>%
  select(marker, species, chromosome, linear_start, linear_end) %>%
  pivot_wider(
    id_cols = marker,
    names_from = species,
    values_from = c(chromosome, linear_start, linear_end),
    names_sep = "_"
  ) %>%
  left_join(alg_df2, by="marker") %>%
  relocate(ALG,.after=marker)

################################################################################
# SPECIES NAMES
################################################################################

sp_x <- "genome_1"
sp_y <- "genome_2"

################################################################################
# CHROMOSOME BOUNDARIES
################################################################################

chr_info_x <- df_lin %>%
  filter(species==sp_x) %>%
  distinct(chromosome, chromosome_size_b) %>%
  mutate(cum_end=cumsum(chromosome_size_b))

chr_info_y <- df_lin %>%
  filter(species==sp_y) %>%
  distinct(chromosome, chromosome_size_b) %>%
  mutate(cum_end=cumsum(chromosome_size_b))

################################################################################
# CHROMOSOME LABELS
################################################################################

chr_labels_x <- df_lin %>%
  filter(species==sp_x) %>%
  group_by(chromosome) %>%
  summarise(
    mid = mean(c(min(linear_start), max(linear_end))),
    .groups="drop"
  ) %>%
  mutate(label = name_map_genome_1[chromosome])

chr_labels_y <- df_lin %>%
  filter(species==sp_y) %>%
  group_by(chromosome) %>%
  summarise(
    mid = mean(c(min(linear_start), max(linear_end))),
    .groups="drop"
  ) %>%
  mutate(label = name_map_genome_2[chromosome])

################################################################################
# COLUMN NAMES FOR PLOTTING
################################################################################

start_x <- paste0("linear_start_", sp_x)
start_y <- paste0("linear_start_", sp_y)

buffer <- 1e5

################################################################################
# PLOT
################################################################################

p <- ggplot(df_wide) +
  
  geom_point(aes_string(
    x = paste0(start_x, "-buffer"),
    y = paste0(start_y, "-buffer"),
    fill = "ALG",
    colour = "ALG"
  ), size = 2.5) +
  
  geom_vline(data = chr_info_x,
             aes(xintercept = cum_end),
             color = "black",
             linewidth = 0.5) +
  
  geom_hline(data = chr_info_y,
             aes(yintercept = cum_end),
             color = "black",
             linewidth = 0.5) +
  
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  
  scale_x_continuous(
    expand = c(0, 0),
    breaks = chr_labels_x$mid,
    labels = chr_labels_x$label,
    name = "Drosophila melanogaster"
  ) +
  
  scale_y_continuous(
    expand = c(0, 0),
    breaks = chr_labels_y$mid,
    labels = chr_labels_y$label,
    name = "Phortica okadai"
  ) +
  
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    
    axis.text.x = element_text(angle = 0, size = 18),
    axis.text.y = element_text(size = 18),
    
    axis.title.x.bottom = element_text(size = 18, face = "italic"),
    axis.title.y.left   = element_text(size = 18, face = "italic"),
    
    legend.position = "none"
  )

p

################################################################################
# SAVE PLOT (same layout as your other figures)
################################################################################

ggsave("Dotplot_2x_asilidae.png",
       plot = p,
       width = 6,
       height = 5,
       units = "in",
       dpi = 300)

################################################################################
# invert chromosomes
################################################################################

#______genome_1_________________________________________________________________

# Syrphidae: C. berberina
chr_to_flip <- "OU862872.1" # 3 of C. berberina
chr_to_flip <- "OU862870.1" # 1 of C. berberina
chr_to_flip <- "OU862871.1" # 2 of C. berberina

# Tachinidae: G. viridis
chr_to_flip <- "OY101440.1" # 1 of G. viridis
chr_to_flip <- "OY101442.1" # 3 of G. viridis
chr_to_flip <- "OY101441.1" # 2 of G. viridis

# Asilidae: N. cyanurus
chr_to_flip <- "OX384542.1" # 3 of N. cyanurus

# get chromosome bounds from genome_1 coordinates in df_wide
chr_bounds <- df_wide %>%
  filter(chromosome_genome_1 == chr_to_flip) %>%
  summarise(
    chr_min = min(linear_start_genome_1, na.rm = TRUE),
    chr_max = max(linear_start_genome_1, na.rm = TRUE)
  )

min_pos <- chr_bounds$chr_min
max_pos <- chr_bounds$chr_max

# reverse positions
df_wide <- df_wide %>%
  mutate(
    linear_start_genome_1 = if_else(
      chromosome_genome_1 == chr_to_flip,
      (max_pos + min_pos) - linear_start_genome_1,
      linear_start_genome_1
    ),
    linear_end_genome_1 = if_else(
      chromosome_genome_1 == chr_to_flip,
      (max_pos + min_pos) - linear_end_genome_1,
      linear_end_genome_1
    )
  )

#______genome_2_________________________________________________________________

# Syrphidae: C. pagana
chr_to_flip <- "OW386170.1" # 5 of C. pagana
chr_to_flip <- "OW386168.1" # 3 of C. pagana
chr_to_flip <- "OW386169.1" # 4 of C. pagana
chr_to_flip <- "OW386167.1" # 2 of C. pagana
chr_to_flip <- "OW386166.1" # 1 of C. pagana

# Drosophilidae: P. okadai
chr_to_flip <- "CM067848.1" # 1 of P. okadai

# Tachinidae: E. cussincta
chr_to_flip <- "OW052034.1" # 1 of E. succincta
chr_to_flip <- "OW052038.1" # 5 of E. succincta

# Asilidae: L. brevirostris
chr_to_flip <-"OZ038374.1" # 4 of L. brevirostris
chr_to_flip <-"OZ038373.1" # 3 of L. brevirostris
chr_to_flip <-"OZ038371.1" # 1 of L. brevirostris

# get chromosome bounds from genome_1 coordinates in df_wide
chr_bounds <- df_wide %>%
  filter(chromosome_genome_2 == chr_to_flip) %>%
  summarise(
    chr_min = min(linear_start_genome_2, na.rm = TRUE),
    chr_max = max(linear_start_genome_2, na.rm = TRUE)
  )

min_pos <- chr_bounds$chr_min
max_pos <- chr_bounds$chr_max

# reverse positions
df_wide <- df_wide %>%
  mutate(
    linear_start_genome_2 = if_else(
      chromosome_genome_2 == chr_to_flip,
      (max_pos + min_pos) - linear_start_genome_2,
      linear_start_genome_2
    ),
    linear_end_genome_2 = if_else(
      chromosome_genome_2 == chr_to_flip,
      (max_pos + min_pos) - linear_end_genome_2,
      linear_end_genome_2
    )
  )

#_______________________________________________________________________________
