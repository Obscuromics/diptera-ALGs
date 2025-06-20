# colour linkage groups at internal nodes with ALGs
# /lustre/scratch126/tol/teams/jaron/projects/diptera_alg/03_syngraph_all_11/output/all_26_new/r2.quick/minus26_89/sorted/mega_file_new.tsv
################################################################################
library(dplyr)
library(ggplot2)
library(gridExtra)
################################################################################
# mega_new contains information for all internal nodes
# its called new because also previously missing salmon BUSCOs are included
mega_new <- read.table("C:/Users/julia/Documents/Kamil/muller/mega_file_new.tsv", header = FALSE, sep = "\t", quote = "", fill = TRUE)
split_columns <- do.call(rbind, strsplit(mega_new$V4, " "))
mega_new$V4 <- split_columns[, 1]
mega_new$V5 <- split_columns[, 2]
colnames(mega_new) <- c("busco", "chrR", "Rstart", "Rend", "n")
head(mega_new)
################################################################################
# split mega_new into separate internal nodes
unique_n <- unique(mega_new$n) # 287 nodes
for (val in unique_n) {
  sub_df <- subset(mega_new, n == val)
  var_name <- paste0("df_", val)
  assign(var_name, sub_df)
}
################################################################################
# create plots
unique_n <- unique(mega_new$n)
unique_n <- unique_n[unique_n != "n3"]

for (val in unique_n) {
  sub_df <- subset(mega_new, n == val)
  
  combinations <- list(val = sub_df)
  names(combinations) <- val
  
  plots <- list()
  for (name in names(combinations)) {
    plots[[name]] <- create_stacked_barplot(n3, combinations[[name]], title = name)
  }
  output_file <- paste0("plot_", val, ".png")
  png(output_file, width = 1000, height = 800)
  grid.arrange(grobs = plots, ncol = 1)
  dev.off()
}
################################################################################
################################################################################
# in this part I have investigated whether the inclusion of vanished salmon BUSCOs
# has changed the reference/ ALGs
# answer: no it did not
################################################################################
# reference: ALG node
# old reference, before finding salmon ALGs
n3 <- read.table("C:/Users/julia/Documents/Kamil/muller/minus26_89_n3.tsv", header = TRUE, sep = "\t")
head(n3)
# new reference
head(df_n3) # now 2993 instead of 2965 BUSCOs assigned to this node
################################################################################
# merge old and new reference

# find busco values unique to n3
unique_n3 <- setdiff(n3$busco, df_n3$busco)
# find busco values unique to df_n3
unique_df_n3 <- setdiff(df_n3$busco, n3$busco)

# Subset the rows with unique busco values
rows_unique_to_n3 <- n3[n3$busco %in% unique_n3, ]
rows_unique_to_df_n3 <- df_n3[df_n3$busco %in% unique_df_n3, ]

# Print results
print("Rows unique to n3:")
print(rows_unique_to_n3) # 10 BUSCOs are unique to n3
# assigned to M1, M4 or M5

print("Rows unique to df_n3:")
print(rows_unique_to_df_n3) # 38 BUSCOs are unique to df_n3
# assigned mostly to M3 (salmon ALG, and 4x to M1)
################################################################################
# Find common busco values
common_busco <- intersect(n3$busco, df_n3$busco)

# Merge dataframes by busco to compare chrR values
merged_df <- merge(n3, df_n3, by = "busco", suffixes = c("_n3", "_df_n3"))

# Filter rows where chrR differs
diff_chrR <- merged_df[merged_df$chrR_n3 != merged_df$chrR_df_n3, ]

# Print results
print("Rows where busco appears in both, but chrR differs:")
print(diff_chrR)

# if a BUSCO is assigned in both run, its always assigned to the same group
# n3 can further be used as reference
################################################################################
create_stacked_barplot <- function(n3, n13, title) {
  if (!is.data.frame(n3) || !is.data.frame(n13)) {
    stop("n3 und n13 müssen Datenrahmen sein")
  }
  groups_n13 <- unique(n13$chrR)
  
  list_of_dataframes <- list()
  for (group in groups_n13) {
    genes_in_n13 <- n13$busco[n13$chrR == group]
    chrR_in_n3 <- n3$chrR[match(genes_in_n13, n3$busco)]
    df <- data.frame(Gene = genes_in_n13, ChrR_in_n3 = chrR_in_n3)
    list_of_dataframes[[group]] <- df
  }
  
  list_of_dataframes <- lapply(list_of_dataframes, function(df) {
    df <- df[!is.na(df$ChrR_in_n3), ]
    return(df)
  })
  
  counts_per_group <- lapply(list_of_dataframes, function(df) {
    table(df$ChrR_in_n3)
  })
  
  df_list <- lapply(names(counts_per_group), function(group) {
    data.frame(
      Group = group,
      Category = names(counts_per_group[[group]]),
      Value = as.numeric(counts_per_group[[group]]),
      stringsAsFactors = FALSE
    )
  })
  
  df <- bind_rows(df_list)
  
  #col_list <- c(M1 = "springgreen3", M2 = "steelblue2", M3 = '#F3D8C7', M4 = "yellow4", M5 = "#FF7B9C")
  col_list <- c(M1 = "yellow4", M2 = "springgreen3", M3 = '#F3D8C7', M4 = "steelblue2", M5 = "#FF7B9C")
  
  plot <- ggplot(df, aes(x = Group, y = Value, fill = Category)) +
    geom_bar(stat = "identity", position = "stack", width = 1) +
    scale_fill_manual(values = col_list) +
    ylim(0, 2000) +
    ggtitle(title) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 90, face = "bold")
    )
  
  return(plot)
}
################################################################################
























