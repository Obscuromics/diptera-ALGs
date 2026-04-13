# length distribution of small dot across phylogeny
################################################################################
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
################################################################################
# small dot status
file_path <- "C:/Users/julia/Documents/Kamil/txt/small_dot_status.csv"
small_dot_defin <- read.csv(file_path, header = FALSE, sep = ",")
small_dot_defin_cp <- small_dot_defin
colnames(small_dot_defin_cp) <- c("label", "status","position")
head(small_dot_defin_cp)

# lengths of small dots
file_path <- "C:/Users/julia/Documents/Kamil/txt/filtered_potential_small_dots_with_length.txt"
small_d_table <- read.delim(file_path, header = TRUE, sep = "\t")
head(small_d_table)

# more info for merging
names_acc_fam_2 <- read.table("C:/Users/julia/Documents/Kamil/busco_fulltable/names_acc_fam_2.csv", sep = ",", header = TRUE)
head(names_acc_fam_2)

# preparing table with lengths to merge
small_d_table <- small_d_table %>%
  mutate(Accession = str_remove(Genome, "\\..*\\.tsv"))
# add column with label (names)
small_d_table <- small_d_table %>%
  left_join(names_acc_fam_2 %>% select(Accession, label), by = "Accession")
head(small_d_table)

# merge length to status info
small_dot_defin_cp <- small_dot_defin_cp %>%
  left_join(small_d_table %>% select(label, Chromosome_Length), by = "label")
head(small_dot_defin_cp)
na_count <- sum(is.na(small_dot_defin_cp$Chromosome_Length)) # 77, therefore warning message
# manually set to NA, there is something wrong with the assembly
small_dot_defin_cp <- small_dot_defin_cp %>%
  mutate(Chromosome_Length = if_else(label == "Bactrocera_correcta", NA_real_, Chromosome_Length))
################################################################################
# create plot
ggplot(small_dot_defin_cp, aes(x = reorder(label, -position), y = Chromosome_Length)) +
  geom_bar(stat = "identity", fill = "black") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = NULL) +
  scale_y_continuous(labels = label_number(scale = 1e-6), 
                     breaks = c(50e6, 100e6, 150e6, 200e6)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(color = "grey75", size = 0.5),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

ggsave("ALG_6_length_tree.png", plot = last_plot(), width = 3, height = 16, dpi = 600)
################################################################################



