

#______________________________________________________________________________#
library(dplyr)
#______________________________________________________________________________#


#______________________________________________________________________________#
# centromere positions, obtained from published data

centromeres <- list(
  # B. dorsalis
  GCF_023373825.1 = list(
    c1 = c(68721987,72242337),
    c2 = c(68087510, 78050978),
    c3 = c(44191285, 47030615),
    c4 = c(57308489, 68008738),
    c5 = c(59225403, 66010094),
    c6 = c(34258598, 56875454)
  ),
  
  # B. pandia
  GCA_916048285.2 = list(
    OU696529.1 = c(82880035, 94090860),
    OU696530.1 = c(62064680, 75178345),
    OU696531.1 = c(68676014, 77786097),
    OU696532.1 = NULL,
    OU696533.1 = c(57540376, 68047747),
    OU696534.1 = NULL
  ),
  
  # Pherbina coryleti
  GCA_943735915.1 = list(
    OX030948.1 = c(94460359, 134986254),
    OX030949.1 = c(108362287, 126272581),
    OX030950.1 = c(76357331, 89267023),
    OX030951.1 = c(86086713, 103704457),
    OX030952.1 = c(44539675, 58051709),
    OX030953.1 = c(200853, 2912352)
  ),
  
  # Phyto melanocephala
  GCA_941918925.1 = list(
    OW799233.1 = c(59521213, 78027807),
    OW799234.1 = c(27720101, 42230621),
    OW799235.1 = c(41122629, 55030281),
    OW799236.1 = c(45933060, 62044654),
    OW799237.1 = c(43220832, 57427679),
    OW799238.1 = c(1304798, 4416235)
  ),
  
  # Pollenia amentaria
  GCA_943735925.1 = list(
    OX031006.1 = c(106750087, 125258770),
    OX031007.1 = c(90220578, 105624090),
    OX031008.1 = c(45827093, 54932475),
    OX031009.1 = c(61833755, 78642929),
    OX031010.1 = c(70841411, 92053809),
    OX031011.1 = c(11824998, 25554019)
  ),
  
  # Sarcophaga rosellei
  GCA_930367235.1 = list(
    OV884017.1 = c(67974787, 78686571),
    OV884018.1 = c(28428892, 37337945),
    OV884019.1 = c(47545107, 61558400),
    OV884020.1 = c(26220170, 32625096),
    OV884021.1 = c(39135559, 45341197),
    OV884022.1 = c(6526523, 7036407)
  ),
  
  # Stomorhina lunata
  GCA_933228675.1 = list(
    OW121740.1 = c(63233460, 81743253),
    OW121741.1 = c(48747413, 66764936),
    OW121742.1 = c(81994587, 89102786),
    OW121743.1 = c(27113456, 32916334),
    OW121744.1 = c(71038652, 77041917),
    OW121745.1 = c(100837, 1512542)
  )
)
#______________________________________________________________________________#



#______________________________________________________________________________#
# file paths
dir_path_busco <- "/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes"
file_path_ALG <- file.path(dir_path_busco, "ALGs_syngraph_brachycera.tsv")
#______________________________________________________________________________#



#______________________________________________________________________________#
# load ALG data
ALGs_syngraph_diptera <- read.table(
  file_path_ALG,
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)
colnames(ALGs_syngraph_diptera) <- c("busco", "ALG")
#______________________________________________________________________________#


#______________________________________________________________________________#
# function to split chromosomes at the centromere
get_side <- function(chr, pos, centro_info) {
  if (!chr %in% names(centro_info)) return(NA)
  
  cpos <- centro_info[[chr]]
  if (is.null(cpos)) return(NA)
  
  if (pos < cpos[1]) return("left")
  if (pos > cpos[2]) return("right")
  
  return(NA)
}
#______________________________________________________________________________#



#______________________________________________________________________________#
# function to compute ADS
compute_ads <- function(df, centro_info) {
  
  if (nrow(df) == 0) return(NULL)
  
  chr <- unique(df$chromosome)
  if (!(chr %in% names(centro_info))) return(NULL)
  
  cpos <- centro_info[[chr]]
  if (is.null(cpos)) return(NULL)
  
  # split by centromere START
  centromere_start <- cpos[1]
  
  df$arm <- ifelse(df$position < centromere_start, "left",
                   ifelse(df$position > centromere_start, "right", NA))
  
  df <- df[!is.na(df$arm), ]
  if (nrow(df) == 0) return(NULL)
  
  # two most common ALGs
  alg_counts <- table(df$ALG)
  if (length(alg_counts) < 2) return(NULL)
  
  top_ALGs <- names(sort(alg_counts, decreasing = TRUE))[1:2]
  
  A <- top_ALGs[1]
  B <- top_ALGs[2]
  
  df2 <- df[df$ALG %in% c(A, B), ]
  if (nrow(df2) == 0) return(NULL)
  
  # function to compute fractions
  get_fraction <- function(alg, side) {
    subset <- df2[df2$ALG == alg, ]
    
    total <- nrow(subset)
    if (total == 0) return(0)
    
    sum(subset$arm == side) / total
  }
  
  # fractions
  A_left  <- get_fraction(A, "left")
  A_right <- get_fraction(A, "right")
  B_left  <- get_fraction(B, "left")
  B_right <- get_fraction(B, "right")
  
  # ADS definition
  ADS_left  <- max(A_left, B_left) - min(A_left, B_left)
  ADS_right <- max(A_right, B_right) - min(A_right, B_right)
  
  #--------------------------------------------------
  # CHECK (SYMMETRY)
  #--------------------------------------------------
  if (abs(ADS_left - ADS_right) > 1e-10) {
    warning(
      paste0("ADS mismatch in chromosome ", chr,
             " | left: ", ADS_left,
             " | right: ", ADS_right)
    )
  }
  
  ADS <- max(ADS_left, ADS_right)
  
  data.frame(
    chromosome = chr,
    ALG_A = A,
    ALG_B = B,
    A_left = A_left,
    A_right = A_right,
    B_left = B_left,
    B_right = B_right,
    ADS = ADS
  )
}
#______________________________________________________________________________#



#______________________________________________________________________________#
# loop to generate results
all_results <- list()

for (genome_id in names(centromeres)) {
  
  file <- file.path(dir_path_busco, paste0(genome_id, "_reduced.tsv"))
  if (!file.exists(file)) next
  
  data <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(data) <- c("busco", "chromosome", "position")
  
  centro_info <- centromeres[[genome_id]]
  
  # assign side (for merging with ALG later)
  data$side <- mapply(get_side,
                      data$chromosome,
                      data$position,
                      MoreArgs = list(centro_info = centro_info))
  
  # merge ALG
  data <- merge(data, ALGs_syngraph_diptera, by = "busco")
  
  # remove NA
  data <- data[!is.na(data$side), ]
  
  if (nrow(data) == 0) next
  
  # compute ADS per chromosome
  res <- data %>%
    group_split(chromosome) %>%
    lapply(function(df) compute_ads(df, centro_info)) %>%
    bind_rows()
  
  res$genome <- genome_id
  
  all_results[[genome_id]] <- res
}

ADS_summary <- bind_rows(all_results)
#______________________________________________________________________________#


#______________________________________________________________________________#
# extraction of relevant LG pairs
important_pairs <- data.frame(
  #ALG_pair1 = c("db1b","db3a","db1a","db3b"),
  #ALG_pair2 = c("db2a","db4a","db2b","db4b"),
  ALG_pair1 = c("db1b","db3a","db1a","db3b"),
  ALG_pair2 = c("db2b","db4a","db2a","db4b"),
  stringsAsFactors = FALSE
)

ADS_summary <- ADS_summary %>%
  mutate(
    pair1 = pmin(ALG_A, ALG_B),
    pair2 = pmax(ALG_A, ALG_B),
    pair = paste(pair1, pair2, sep = "_")
  )

important_pairs$pair <- paste(
  pmin(important_pairs$ALG_pair1, important_pairs$ALG_pair2),
  pmax(important_pairs$ALG_pair1, important_pairs$ALG_pair2),
  sep = "_"
)

# filter important pairs
ADS_filtered <- ADS_summary %>%
  filter(pair %in% important_pairs$pair)

ADS_filtered$pair <- factor(ADS_filtered$pair, levels = important_pairs$pair)
#______________________________________________________________________________#


#______________________________________________________________________________#
# plot
xlabels <- c(
  "db1b_db2b" = "y1", # C
  "db3a_db4a" = "y3",  # E
  "db1a_db2a" = "y2",  # D
  "db3b_db4b" = "y4"   # A
)


png("/Users/jg40/Desktop/ADS_plot.png",
    width = 6, height = 5, units = "in", res = 300)
par(mar = c(5, 5, 1, 1)) # bottom, left, top, right

stripchart(ADS ~ pair,
           data = ADS_filtered,
           method = "jitter",
           jitter = 0.2,
           vertical = TRUE,
           pch = 16,
           col = "steelblue",
           cex = 2,
           ylab = "ADS",
           xlab = "ancestral karyotype y",
           xaxt = "n",
           ylim = c(0, 1),
           cex.lab = 1.5,
           cex.axis = 1.5)

# medians
med <- aggregate(ADS ~ pair, data = ADS_filtered, median)

points(1:nrow(med), med$ADS,
       pch = 95, cex = 4, col = "coral")

axis(1,
     at = 1:nrow(med),
     labels = xlabels[as.character(med$pair)],
     las = 1,
     cex.axis = 1.5)
dev.off()
#______________________________________________________________________________#
