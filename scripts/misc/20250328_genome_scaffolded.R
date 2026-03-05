# how much on everage of each genome is scaffolded into chromosomes?
library(dplyr)
data <- read.delim("C:/Users/julia/Documents/Kamil/txt/all_info_chroms.txt", 
                   header = FALSE, sep = "", check.names = FALSE)
colnames(data) <- c("chromosome", "length", "o", "p", "q", "accession", "identity", "num_BUSCO")
head(data)


result <- data %>%
  group_by(accession) %>%
  summarise(
    total_length = sum(length, na.rm = TRUE),
    chrom_length = sum(length[identity == "chrom"], na.rm = TRUE),
    chrom_percentage = (chrom_length / total_length) * 100
  )

print(result)
mean(result$chrom_percentage)
sd(result$chrom_percentage)



