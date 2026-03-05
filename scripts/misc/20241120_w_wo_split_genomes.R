# comparison of ALGs inference with and without split genomes
# in the version of the tree that I used here
# both node 4 and node 5 contained valid set of ALGs
# in that version I removed 8 genomes
# see 20240910_remove_tips_tree.R
################################################################################
directory <- 'C:/Users/julia/Documents/Kamil/busco_fulltable/all_tol'

busco_files <- paste0(directory, "/", c("ALGs_n4_72_all_split.tsv",
                                        "ALGs_n5_72_all_split.tsv",
                                        "ALGs_n4_72_all_8.tsv",
                                        "ALGs_n5_72_all_8.tsv",
                                        "ALGs_n1_72_all_split.tsv",
                                        "ALGs_n1_72_all_8.tsv"))

chrom_files <- paste0(directory, "/", c("ALGs_n4_72_all_split_chrominf.tsv",
                                        "ALGs_n5_72_all_split_chrominf.tsv",
                                        "ALGs_n4_72_all_8_chrominf.tsv",
                                        "ALGs_n5_72_all_8_chrominf.tsv",
                                        "ALGs_n1_72_all_split_chrominf.tsv",
                                        "ALGs_n1_72_all_8_chrominf.tsv"))
################################################################################
busco_data <- lapply(busco_files, read.table, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
chrom_data <- lapply(chrom_files, read.table, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

names(busco_data) <- basename(busco_files)
names(chrom_data) <- basename(chrom_files)

new_colnames <- c("busco", "chrQ", "Qstart", "Qend")
busco_data <- lapply(busco_data, function(df) {
  colnames(df) <- new_colnames
  return(df)
})

str(busco_data)
################################################################################
# number of common buscos shared by  n4, n5 of both runs
busco_values_1 <- busco_data[[1]]$busco
busco_values_2 <- busco_data[[2]]$busco
busco_values_3 <- busco_data[[3]]$busco
busco_values_4 <- busco_data[[4]]$busco

# calculate common BUSCOs for the four nodes
common_busco_values <- Reduce(intersect, list(busco_values_1, busco_values_2, busco_values_3, busco_values_4))
print(common_busco_values)
print(length(common_busco_values))
################################################################################
# this code produces a matrix for the 4 nodes
# shows if and how BUSCOs change between nodes

busco_sublist <- busco_data[1:4]
count_inconsistent_chrQ <- function(list1, list2) {
  common_busco <- intersect(list1$busco, list2$busco)
  inconsistent_chrQ <- sapply(common_busco, function(busco_val) {
    chrQ_list1 <- list1$chrQ[list1$busco == busco_val]
    chrQ_list2 <- list2$chrQ[list2$busco == busco_val]
    return(chrQ_list1 != chrQ_list2)
  })
  return(sum(inconsistent_chrQ))
}

inconsistent_counts <- matrix(0, nrow = 4, ncol = 4)

for (i in 1:4) {
  for (j in i:4) {
    inconsistent_counts[i, j] <- count_inconsistent_chrQ(busco_sublist[[i]], busco_sublist[[j]])
    if (i != j) {
      inconsistent_counts[j, i] <- inconsistent_counts[i, j]  # Symmetrisch
    }
  }
}

inconsistent_counts
# this shows the matrix with BUSCOs that changed
# but does not show if BUSCOs were lost
################################################################################
# this code shows names of BUSCOs that change
# its always the same five that change between n4 and n5
# but between the two different n4 and n5 nothing changes

# example for the two n4 nodes
list1 <- busco_data[["ALGs_n4_72_all_split.tsv"]]
list2 <- busco_data[["ALGs_n4_72_all_8.tsv"]]

# finds BUSCOs that are common for two nodes
common_busco <- intersect(list1$busco, list2$busco)

# compare ALG assignment for common BUSCO values
inconsistent_chrQ <- sapply(common_busco, function(busco_val) {
  chrQ_list1 <- list1$chrQ[list1$busco == busco_val]
  chrQ_list2 <- list2$chrQ[list2$busco == busco_val]
  return(chrQ_list1 != chrQ_list2)
})

inconsistent_busco <- common_busco[inconsistent_chrQ]
inconsistent_busco
################################################################################

































