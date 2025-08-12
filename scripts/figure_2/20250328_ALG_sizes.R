# barplot that shows sizes and colours of ALGs

library(ggplot2)

ALG_data <- read.table('data/diptera.no_plecia.mindist.m165_n1_n2.tsv', col.names = c("busco", "chrQ"))

head(ALG_data)

# pal <- c("M1" = "#1573afff", "M2" = "#e59d38ff", "M3" = "#f0e354ff", 
#          "M4" = "#169e73ff", "M5" = "#60b5e1ff", "M6" = "#000000ff")
source('scripts/20250620_colour_pal.R')

ggplot(ALG_data, aes(x = factor(chrQ), fill = chrQ)) +
  geom_bar() +
  scale_fill_manual(values = pal) +
  xlab("ALGs") +
  ylab("BUSCOs") +
  scale_x_discrete(labels = c("d1" = "1", "d2" = "2", "d3" = "3", "d4" = "4", "d5" = "5", "d6" = "6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

ggsave("figures/ALG_sizes_plot.png", plot = last_plot(), dpi = 600, width = 8, height = 6)

