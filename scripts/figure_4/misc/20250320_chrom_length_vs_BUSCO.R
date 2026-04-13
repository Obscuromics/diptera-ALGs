# plot length of each chromosome against number of BUSCOs
# color the 212 chromosomes that are syntenic to ALG 6 differently
################################################################################
library(ggplot2)
library(dplyr)
################################################################################
# I collected the data from my files on the cluster
data <- read.delim("C:/Users/julia/Documents/Kamil/txt/all_info_chroms.txt", 
                   header = FALSE, sep = "", check.names = FALSE)
colnames(data) <- c("chromosome", "length", "o", "p", "q", "accession", "identity", "num_BUSCO")

removed_spec <- read.delim("C:/Users/julia/Documents/Kamil/txt/removed.txt", 
                   header = FALSE, sep = "", check.names = FALSE)

data$accession_clean <- sub("\\.\\d+$", "", data$accession)
data_filtered <- data[!(data$accession_clean %in% removed_spec$V1) & data$num_BUSCO > 4, ]
head(data_filtered)
################################################################################
# extract info to integrate into master sheet
busco_summary <- data %>%
  group_by(accession) %>%
  summarise(
    total_num_BUSCO_chrom = sum(num_BUSCO[identity == "chrom"], na.rm = TRUE),
    total_num_BUSCO_scaff = sum(num_BUSCO[identity == "scaff"], na.rm = TRUE),
    total_length = sum(length[identity == "chrom"], na.rm = TRUE)
  )
busco_summary <- as.data.frame(busco_summary)
#write.table(busco_summary, "busco_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
################################################################################
# plot length of chromosomes against BUSCO count
ALG_6_chroms <- read.delim("C:/Users/julia/Documents/Kamil/txt/ALG_6_chroms.txt", sep = "\t", header = TRUE)
head(ALG_6_chroms)

ALG_6_chroms_vec <- ALG_6_chroms[[1]]

filtered_data <- data_filtered %>%
  filter(identity == "chrom") %>%
  mutate(ALG_6_syntenic = ifelse(chromosome %in% ALG_6_chroms_vec, "ALG 6 syntenic chromosomes", "chromosomes"))

lm_model <- lm(num_BUSCO ~ length, data = filtered_data)
r_squared <- summary(lm_model)$r.squared
r <- ggplot(filtered_data, aes(x = length, y = num_BUSCO, color = ALG_6_syntenic)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("chromosomes" = "black", "ALG 6 syntenic chromosomes" = "blue")) +
  labs(x = "Chromosome Length (Mb)", y = "BUSCOs") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  xlim(0, 250000000) +  # to zoom into blue dots
  ylim(0, 100) +  # to zoom into blue dots
  scale_x_continuous(labels = scales::label_number(scale = 1e-6)) +
  guides(color = guide_legend(title = NULL)) #+
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  #annotate("text", x = 400000000, y = 80, label = paste("R² = ", round(r_squared, 2)), size = 3.5, color = "black")
print(r)

ggsave("chrom_length_vs_BUSCOs_zoomed.png", plot = r, dpi = 600, width = 8, height = 6)
################################################################################
# obtain R2
model <- lm(filtered_data$num_BUSCO ~ filtered_data$length)
summary(model)$r.squared
# 0.2688306
# 0.2584953 # zoomed in version
################################################################################
# investigate chromosomes, that are short and BUSCO poor but not related to ALG 6
filtered_data <- data %>%
  filter(identity == "chrom", num_BUSCO <= 100)

filtered_data <- filtered_data %>%
  mutate(is_in_alg_6 = ifelse(chromosome %in% ALG_6_chroms_vec, "match", "no_match"))

no_match_data <- filtered_data %>%
  filter(is_in_alg_6 == "no_match")
no_match_data
# excluded species or fusions
################################################################################

