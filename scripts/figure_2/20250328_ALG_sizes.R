# barplot that shows sizes and colours of ALGs

library(ggplot2)

n3 <- read.table("C:/Users/julia/Documents/Kamil/muller/minus26_89_n3.tsv", header = TRUE, sep = "\t")
ALG_6 <- read.table("C:/Users/julia/Documents/Kamil/muller/small_dot_ALG.txt", header = TRUE, sep = "\t")
ALG_data <- rbind(n3, ALG_6)
colnames(ALG_data) <- c("busco", "chrQ", "Qstart", "Qend")
head(ALG_data)

pal <- c("M1" = "#1573afff", "M2" = "#e59d38ff", "M3" = "#f0e354ff", 
         "M4" = "#169e73ff", "M5" = "#60b5e1ff", "M6" = "black")

ggplot(ALG_data, aes(x = factor(chrQ, levels = c("M4", "M2", "M1", "M3", "M5", "M6")), fill = chrQ)) +
  geom_bar() +
  scale_fill_manual(values = pal) +
  xlab("ALGs") +
  ylab("BUSCOs") +
  scale_x_discrete(labels = c("M4" = "1", "M2" = "2", "M1" = "3", "M3" = "4", "M5" = "5", "M6" = "6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

ggsave("ALG_plot.png", plot = last_plot(), dpi = 600, width = 8, height = 6)

