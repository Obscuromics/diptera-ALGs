# plot genome of each species colored in ALGs
# one plot per family
################################################################################
library(tidyr)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(grid)
library(ggtree)
################################################################################
# tree with 289 species (288 + outgroup)
tree <- read.tree("C:/Users/julia/Documents/Kamil/trees/tree_minus26.treefile")
num_tips <- Ntip(tree)
tree2 <- ggtree(tree)
tree_data <- tree2$data

# data frame with species names, accessions and family name in order of the big, original tree
names_acc_fam_2 <- read.table("C:/Users/julia/Documents/Kamil/busco_fulltable/names_acc_fam_2.csv", sep = ",", header = TRUE)
colnames(names_acc_fam_2) <- c("names", "num", "acc", "fam")
head(names_acc_fam_2)

# extract names on current tips
tip_names <- tree$tip.label
head(tip_names)

# get the right order of the species on the current tree
tree_data_filtered <- tree_data %>%
  filter(isTip == TRUE) %>%
  arrange(angle) # angle lists species in correct order

# extract labels
label_df <- tree_data_filtered %>%
  select(label) %>%
  mutate(counter = seq_along(label))

# merge with info about accession and family
combined_df <- inner_join(label_df, names_acc_fam_2, by = c("label" = "acc")) 
subset_names_acc_fam_2 <- combined_df

# generate two columns to get files that are necessary for the analysis
subset_names_acc_fam_2$busco_files <- paste0(subset_names_acc_fam_2$label, ".tsv")
subset_names_acc_fam_2$chrom_files <- paste0(subset_names_acc_fam_2$label, "_chrominf.tsv")
head(subset_names_acc_fam_2)

# extract names in the correct order for small dot purposes
names_column <- subset_names_acc_fam_2$names
names_df <- data.frame(names = names_column)
#write.csv(names_df, "names_column.csv", row.names = FALSE)

diptera_info <- subset_names_acc_fam_2[, c("label", "names", "fam", "counter")]
head(diptera_info)
#write.csv(diptera_info, "diptera_info.csv", row.names = FALSE)
################################################################################
# extract unique family names
unique_fams <- unique(subset_names_acc_fam_2$fam)
length(unique_fams) # 53 --> 52 Diptera families and outgroup family
# generate sub-dataframes for each family, name accordingly
for (i in seq_along(unique_fams)) {
  fam_name <- unique_fams[i]
  sub_df <- subset(subset_names_acc_fam_2, fam == fam_name)
  new_name <- paste0(i, "_", fam_name)  # number_accession
  assign(new_name, sub_df)
}
################################################################################
#  loop over families to generate plot for species within that family
for (i in seq_along(unique_fams)) {
  fam_name <- unique_fams[i]
  sub_df_name <- paste0(i, "_", fam_name)
  if (exists(sub_df_name)) {
    sub_df <- get(sub_df_name)
    busco_files <- sub_df$busco_files
    chrom_files <- sub_df$chrom_files
    all_genome_plots <- list()
    
    for (j in seq_along(busco_files)) {
      busco_file_path <- file.path(directory, busco_files[j])
      busco_data <- read.table(busco_file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
      colnames(busco_data) <- new_colnames
      
      busco_data_ALG <- merge(busco_data, ALG_data[, c("busco", "chrQ")], by = "busco", all.x = TRUE)
      busco_data_ALG <- busco_data_ALG[!is.na(busco_data_ALG$chrQ.y), ]

      chrom_file_path <- file.path(directory, chrom_files[j])
      chrom_data <- read.table(chrom_file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
      
      chromosomes <- unique(busco_data_ALG$chrQ.x)
      chromosomes <- sort(chromosomes)
      reordered_chromosomes <- character(length(chromosomes))
      right_order <- chrom_data$order
      reordered_chromosomes[right_order] <- chromosomes
      chromosomes <- reordered_chromosomes
      len <- length(chromosomes)
      
      chromosome_data_sorted <- busco_data_ALG %>%
        filter(!is.na(chrQ.y)) %>%
        arrange(factor(chrQ.x, levels = chromosomes), Qstart) %>%
        group_by(chrQ.x) %>%
        mutate(group = ceiling(row_number() / 20))
      
      proportions <- chromosome_data_sorted %>%
        group_by(chrQ.x, group, chrQ.y) %>%
        summarise(count = n(), .groups = "drop") %>%
        group_by(chrQ.x, group) %>%
        mutate(prop = count / sum(count))
      
      plots <- proportions %>%
        split(.$chrQ.x) %>%
        lapply(function(chrom_data) {
          ggplot(chrom_data, aes(x = factor(group), y = prop, fill = chrQ.y)) +
            geom_bar(stat = "identity") +
            scale_fill_manual(values = c("M1" = "yellow4", "M2" = "springgreen3", "M3" = "#F3cdC7", "M4" = "steelblue2", "M5" = "#FF7B9C", "M6" = "black")) +
            labs(
              title = paste0(unique(chrom_data$chrQ.x), "_", sub_df$names[j]),
              x = NULL,
              y = NULL,
              fill = "Category"
            ) +
            theme_void() +
            theme(plot.title = element_text(size = 36, face = "bold"),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  legend.position = "none",
                  plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"))
        })
      
      valid_plots <- plots[chromosomes]
      genome_plot <- grid.arrange(grobs = valid_plots, ncol = len)
      
      all_genome_plots[[j]] <- genome_plot
    }
    
    output_file <- paste0(i, "_", fam_name, ".png")
    png(output_file, width = 4500, height = 500 * length(busco_files))
    grid.arrange(grobs = rev(all_genome_plots), ncol = 1)
    dev.off()
    
    cat("Finished plotting for family:", fam_name, "\n")
  }
}
################################################################################
