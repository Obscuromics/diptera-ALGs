require('ape')
require('ggtree')
require('ggplot2')
require('dplyr')
require('phytools')
require('gsheet')
require("ggpubr")
require("viridis")
require("patchwork")
require("stringr")
require("tibble")
require('reshape2')
source('scripts/20250620_colour_pal.R')

#load data
tree <- read.tree("data/diptera.supermatrix.phy.treefile") #tree
all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?usp=sharing", format='csv'),
                            stringsAsFactors = F, header = T, check.names = F)
genome_data <- all_genome_data[all_genome_data[, 'TO ADD'] == 'KEEP', ]
chromosome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?gid=1940964825#gid=1940964825", format='csv'),
                            stringsAsFactors = F, header = T, check.names = F)

dir <- 'data/busco_tables' # busco table dir
ALGs <- read.table('tables/ALGs_syngraph_diptera.tsv') # ALG definitions
row.names(ALGs) <- ALGs[, 1]
files <- list.files(path=dir, pattern="syngraph.buscos.tsv", full.names=TRUE, recursive=FALSE)

# label x chromosomes
x_chromosomes <- c(genome_data$X_chrom_NCBI_1, genome_data$X_chrom_NCBI_2)
chromosome_data$is_X <- chromosome_data$chromosome %in% x_chromosomes

## Tree
scale_factor <- 8
outgroup <- "Panorpa_germanica"
tip_labels <- tree$tip.label
rooted_tree <- reroot(tree, node.number = which(tip_labels == outgroup), position =  node.depth.edgelength(tree)[which(tip_labels == outgroup)] / scale_factor)
tree_show <- ggtree(rooted_tree)+ 
              theme_tree2() + 
              scale_x_continuous(labels = abs) + 
              xlab('')

tree_show <- revts(tree_show)
data <- tree_show$data
sorted_data_tips_desc <- data %>% filter(isTip == TRUE) %>% arrange(-y)
colnames(sorted_data_tips_desc)[4] <- 'species'

## get ALG counts per chromosome
alg_per_chromosome <- data.frame()
get_ALGs_per_chromosome <- function(file){
	buscos <- read.table(file, sep = '\t')
	buscos[, 'ALG'] <- ALGs[buscos[, 1], 2]
	buscos <- buscos[!is.na(buscos[, 'ALG']), ]
	species_df <- as.data.frame.matrix(table(buscos[, 2], buscos[, 'ALG']))
	alg_per_chromosome <- rbind(alg_per_chromosome, species_df)
}
alg_per_chromosome <- lapply(files, get_ALGs_per_chromosome)
alg_per_chromosome <- bind_rows(alg_per_chromosome)
alg_per_chromosome <- tibble::rownames_to_column(alg_per_chromosome, "chromosome")

chromosome_data <- merge(chromosome_data, alg_per_chromosome, on='chromosome', all=T)
chromosome_data <- merge(chromosome_data, sorted_data_tips_desc[c('species','y')], on='species')
chromosome_data <- chromosome_data %>% arrange(-y)

## get d6 statistics and labels
### label chromosome with most d6 markers
d_cols <- c('d1','d2','d3','d4','d5','d6')
chromosome_data$alg_sum <- rowSums(chromosome_data[, d_cols])
chromosome_data$d6_prop <- ifelse(chromosome_data$alg_sum > 0, chromosome_data$d6 / chromosome_data$alg_sum, 0)
chromosome_data[d_cols] <- lapply(chromosome_data[d_cols], function(x){ x[is.na(x)] <- 0; x })
chromosome_data <- chromosome_data %>%
  group_by(species) %>%
  mutate(
    d6_max = d6 == max(d6),
    d6_prop_max = d6_prop == max(d6_prop)
  ) %>%
  ungroup()

# check for sensitivity to threshold
results <- c()
for (i in c(1:99)){
	threshold <- i/100
	res <- chromosome_data %>%
  	mutate(dot_d6 = d6_prop > threshold)
  	count <- table(res$dot_d6)[2]
  	results <- append(results, count)
}

results <- unname(results)
png("figures/dot_diag_plot.png")
plot(c(1:99), results)
dev.off()
png("figures/dot_hist_plot.png")
plot(hist(chromosome_data$d6_prop))
dev.off()

threshold <- 0.5
chromosome_data <- chromosome_data %>%
  mutate(dot_d6 = d6_prop > threshold,
         dot_sc = is_X & dot_d6)

# any species with multiple candidates?
dot_count_table <- table(chromosome_data[chromosome_data$dot_d6 == TRUE,]$species)
dot_count_table[dot_count_table>1]
# with threshold 0.2
# Dolichopus_virgultorum    Portevinia_maculata 
#                     2                      2 
chromosome_data[chromosome_data$species == 'Dolichopus_virgultorum',]
chromosome_data[chromosome_data$species == 'Portevinia_maculata',]
# Dolichopus_virgultorum has 2 dots, Portevinia has a 3 marker dot thats still chunky?
# P maculata interesting but no point showing it so could drop it?
# need to figure out how to display both of these

get_wastrels <- function(data_df, columns) {
	wastrels <- unique(chromosome_data$species[which(!chromosome_data$species %in% data_df$species)])
	wastrel_df = data.frame(matrix(nrow = length(wastrels), ncol = length(columns))) 
	colnames(wastrel_df) = columns
	wastrel_df$species <- wastrels
	data_df <- rbind(as.data.frame(data_df), wastrel_df)
	data_df <- merge(data_df, sorted_data_tips_desc[c('species','y')])
	data_df <- data_df %>%
		arrange(-y)
	return(data_df)
}

## dot size distribution
dot_size_df <- chromosome_data %>%
	filter(dot_d6 == TRUE)
size_columns <- c('chromosome','species','chromosome_size_b')
dot_size_df <- dot_size_df[size_columns]
dot_size_df <- get_wastrels(dot_size_df, size_columns)
dot_size_df$chromosome[is.na(dot_size_df$chromosome)] <- 'NA'
dot_size_df$chromosome_size_b[is.na(dot_size_df$chromosome_size_b)] <- 0

dot_size_plot <- ggplot(dot_size_df, aes(x = chromosome_size_b, y = factor(y, levels = rev(unique(y))))) +
			  geom_bar(stat = "identity") +
			  xlab("ALG6* size (Mbp)") +
			  scale_x_continuous(
			    breaks = seq(0, max(dot_size_df$chromosome_size_b, na.rm = TRUE), by = 1e6),
			    labels = scales::label_number(scale = 1e-6)
			  ) +
			  theme_minimal() +
			  theme(panel.grid = element_blank(),
			        #panel.grid.major.x = element_line(color = "gray", size = 0.5),
			        axis.text.y = element_blank(),
			        axis.title.y = element_blank(),
			        axis.line.x.bottom = element_line(size=0.5),
			        axis.text.x = element_blank(),
			        axis.title.x = element_blank(),
			        legend.position = "none")

# painting function
plot_paintings <- function(data_df){
	columns <- c('species', 'd1','d2','d3','d4','d5','d6')
	data_df <- data_df[columns]
	data_df <- get_wastrels(data_df, columns)
	molten_df <- melt(data_df, id.vars=c('species','y'),
	            measure.vars=columns[2:7])
	molten_df$value[is.na(molten_df$value)] <- 0
	molten_df$variable <- factor(molten_df$variable, levels = c('d1','d2','d3','d4','d5','d6'))
	painting_plot <- ggplot(molten_df, aes(x=value, y = factor(y, levels = rev(unique(y))), fill=variable)) +
	         theme_minimal() +
	         geom_bar(position="fill", stat="identity") +
	         scale_fill_manual(values = alg_pal)+
			 theme(panel.grid = element_blank(),
			        #panel.grid.major.x = element_line(color = "gray", size = 0.5),
			        axis.text.y = element_blank(),
			        axis.title.y = element_blank(),
			        axis.line.x.bottom = element_line(linewidth=0.5),
			        axis.text.x = element_blank(),
			        axis.title.x = element_blank(),
			        legend.position = "none")
    return(painting_plot)
}

### ALG6* painting column
dot_data <- chromosome_data %>%
	filter(dot_d6 == TRUE)
alg6_plot <- plot_paintings(dot_data)

#sex chr plot
sex_chr_data <- chromosome_data %>%
	filter(is_X == TRUE) %>%
	filter(dot_d6 == FALSE)

dot_sex_chr_data <- chromosome_data %>%
  mutate(dot_sc = is_X & dot_d6) %>%
  group_by(species) %>%
  mutate(has_dot_sc = any(dot_sc, na.rm = TRUE)) %>%
  ungroup() %>%
  select(species, y, has_dot_sc) %>%
  distinct() %>%
  filter(has_dot_sc)

missing_sex_chr_data <- chromosome_data %>%
  group_by(species) %>%
  mutate(no_X = !any(is_X)) %>%
  ungroup() %>%
  select(species, y, no_X) %>%
  distinct() %>%
  filter(no_X)

#proportion of alg6 on X chromosome
d6_by_species <- chromosome_data %>%
  group_by(species) %>%
  summarise(total_d6 = sum(d6, na.rm = TRUE))

x_data <- chromosome_data %>%
	filter(is_X == TRUE) %>%
	arrange(-y)

x_data_d6sum <- x_data %>%
  left_join(d6_by_species, by = "species")

x_data_d6sum$x_d6_prop <- x_data_d6sum$d6/x_data_d6sum$total_d6

y_lookup <- chromosome_data %>%
  select(species, y) %>%
  distinct()

alg6_x_plot_data <- full_join(
  x_data_d6sum,
  missing_sex_chr_data %>% select(species, no_X),
  by = "species"
) %>%
  select(-matches("^y")) %>%   # remove any messy y columns
  left_join(y_lookup, by = "species") %>%
  arrange(desc(y))
  
alg6_x_plot_data <- alg6_x_plot_data %>%
  mutate(x_d6_prop = tidyr::replace_na(x_d6_prop, 0))

alg6_X_prop_plot <- ggplot(alg6_x_plot_data, aes(x = x_d6_prop, y = factor(y, levels = rev(unique(y))))) +
			  geom_bar(stat = "identity") +
			  xlab("Proportion of ALG6 on X chromosome") +
			  scale_x_continuous(
			    breaks = seq(0, max(alg6_x_plot_data$x_d6_prop, na.rm = TRUE), by = 1e6),
			    labels = scales::label_number(scale = 1e-6)
			  ) +
			  theme_minimal() +
			  theme(panel.grid = element_blank(),
			        #panel.grid.major.x = element_line(color = "gray", size = 0.5),
			        axis.text.y = element_blank(),
			        axis.title.y = element_blank(),
			        axis.line.x.bottom = element_line(size=0.5),
			        axis.text.x = element_blank(),
			        axis.title.x = element_blank(),
			        legend.position = "none")

sex_chr_plot <- plot_paintings(sex_chr_data) +
				  geom_point(
    			  data = dot_sex_chr_data,
				    aes(x = 0.5, y = factor(y, levels = rev(unique(y)))),  # match factor levels exactly
				    shape = 21,   # fixed circle
				    fill = "black",
				    color = "black",
				    size = 0.25
				  ) 
				  #geom_point(
    			  #data = missing_sex_chr_data,
				  #  aes(x = 0.5, y = factor(y, levels = rev(unique(y)))),  # match factor levels exactly
				  #  shape = 4,  # 4 for X, 63 for ?
				  #  fill = "red",
				  #  color = "red",
				  #  size = 0.1
				  #)

### Create compound plot
plt_all <- (tree_show | dot_size_plot | alg6_plot | sex_chr_plot | alg6_X_prop_plot) +
            plot_layout(widths = c(6, 1, 1, 1,1))
ggsave("figures/fig4_core.png", plot=plt_all, dpi=600, width = 6, height = 10)
ggsave("figures/fig4_core.svg", plot=plt_all, dpi=600, width = 6, height = 10)
###


