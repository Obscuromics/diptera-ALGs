# plot avergae AII with sd for each ALG across all species
# values obtained from master sheet

averages <- c(0.7088700321, 0.661133294, 0.7278227279, 0.6512500093, 0.7336735769)
std_dev <- c(0.07350246057, 0.09143000156, 0.07693922342, 0.07853134737, 0.06940582159)
ALG <- c("ALG1", "ALG2", "ALG3", "ALG4", "ALG5")
col_list <- c("ALG3" ="#1573afff", "ALG2" = "#e59d38ff", "ALG4" = "#f0e354ff", "ALG1" = "#169e73ff", "ALG5" = "#60b5e1ff")

# Create a data frame
data <- data.frame(ALG, averages, std_dev)

# Plot
library(ggplot2)

ggplot(data, aes(x=ALG, y=averages, color=ALG)) +
  geom_point(size=15.5) +
  geom_errorbar(aes(ymin=averages - std_dev, ymax=averages + std_dev), width=0.2) +
  scale_color_manual(values=col_list) +
  scale_y_continuous(limits=c(0.4, 1.0)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(size=27),
        axis.title.y = element_text(size=25)) +
  labs(y = "AII")

ggsave("AII_average.png", plot = last_plot(), width = 10, height = 5)

