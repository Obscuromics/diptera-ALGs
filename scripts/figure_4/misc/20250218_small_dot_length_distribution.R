# plot length distribution of potential small dot chromosomes
# dowloaded file from lustre
# from diptera_alg/01_busco_fullt/small_dot_archive/length ...

library(ggplot2)
library(scales)

file_path <- "C:/Users/julia/Documents/Kamil/txt/filtered_potential_small_dots_with_length.txt"
small_d_table <- read.delim(file_path, header = TRUE, sep = "\t")

top_10_longest_chromosomes <- small_d_table[order(-small_d_table$Chromosome_Length), ][1:10, ]
head(top_10_longest_chromosomes)

ggplot(small_d_table, aes(x = Chromosome_Length)) +
  geom_histogram(binwidth = 2500000, fill = "#929292", color = "#929292", alpha = 0.7) +
  labs(x = "chromosome length [Mb]", y = "frequency") +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

ggsave("frequency_chrom_length_ALG_6.png", plot = last_plot(), dpi = 600, width = 8, height = 6)
