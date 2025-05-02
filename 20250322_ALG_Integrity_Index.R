# this script calculates AII
# ATTENTION: the files needed are in a different folder
# in: ALG_painting/BUSCO_data
################################################################################
library(tidyr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(gridExtra)
################################################################################
names_acc_fam_2 <- read.table("C:/Users/julia/Documents/Kamil/busco_fulltable/names_acc_fam_2.csv", sep = ",", header = TRUE)
colnames(names_acc_fam_2) <- c("names", "num", "acc", "fam")
head(names_acc_fam_2)
diptera_info <- names_acc_fam_2
colnames(diptera_info) <- c("label", "num", "acc", "fam")
head(diptera_info)

families_order_tree <- read.table("C:/Users/julia/Documents/Kamil/txt/families_order_tree.tsv",
                                  sep = "\t", header = TRUE, stringsAsFactors = FALSE)
################################################################################
# function to prepare/ aggregate ALG data to calculate ALG integrity index
process_busco_data <- function(directory, busco_files, chrom_files, ALG_data, new_colnames) {
  all_results <- list()
  all_freq_data <- list()
  
  for (i in seq_along(busco_files)) {
    busco_file_path <- file.path(directory, busco_files[i])
    chrom_file_path <- file.path(directory, chrom_files[i])
    
    result_key <- tools::file_path_sans_ext(busco_files[i])
    busco_data <- read.table(busco_file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    colnames(busco_data) <- new_colnames
    
    chrom_data <- read.table(chrom_file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    busco_data_ALG <- merge(busco_data, ALG_data[, c("busco", "chrQ")], by = "busco", all.x = TRUE)
    busco_data_ALG <- busco_data_ALG[!is.na(busco_data_ALG$chrQ.y), ]
    
    unique_chroms <- unique(busco_data_ALG$chrQ.x)
    combined_df_list <- list()
    freq_data_list <- list()
    
    for (chrom in unique_chroms) {
      trial_4 <- subset(busco_data_ALG, chrQ.x == chrom)
      trial_4 <- trial_4[order(trial_4$Qstart), ]
      
      only_M <- trial_4[[5]]
      result <- rep(NA, length(only_M))
      counter <- 1
      
      for (j in seq_along(only_M)[-1]) {
        if (only_M[j] == only_M[j - 1]) {
          counter <- counter + 1
          result[j] <- NA
        } else {
          result[j - 1] <- counter
          counter <- 1
        }
      }
      
      result[length(only_M)] <- counter
      output <- data.frame(only_M, result)
      output_clean <- output[!is.na(output$result), ]
      
      freq_data <- as.data.frame(table(output_clean$only_M, output_clean$result))
      colnames(freq_data) <- c("category", "result", "frequency")
      freq_data$result <- as.numeric(as.character(freq_data$result))
      
      freq_data_list[[chrom]] <- freq_data
      
      category_sum <- output_clean %>%
        group_by(only_M) %>%
        summarise(busco_total = sum(result))
      
      total_sum <- sum(category_sum$busco_total)
      category_sum$busco_percentage <- (category_sum$busco_total / total_sum) * 100
      
      weighted_avg_length <- freq_data %>%
        group_by(category) %>%
        summarise(avg_block_size = sum(frequency * result) / sum(frequency)) %>%
        ungroup()
      
      busco_num <- sum(category_sum$busco_total)
      weighted_avg_length$block_percentage <- (weighted_avg_length$avg_block_size / busco_num) * 100
      
      combined_df <- category_sum %>%
        rename(category = only_M) %>%
        left_join(weighted_avg_length, by = "category")
      
      combined_df_list[[chrom]] <- combined_df
    }
    
    all_results[[result_key]] <- combined_df_list
    all_freq_data[[result_key]] <- freq_data_list
  }
  
  all_aggregated_data <- list()
  
  for (i in seq_along(all_freq_data)) {
    freq_data_cur <- all_freq_data[[i]]
    M1_df <- data.frame()
    M2_df <- data.frame()
    M3_df <- data.frame()
    M4_df <- data.frame()
    M5_df <- data.frame()
    
    for (chromosome in names(freq_data_cur)) {
      freq_data <- freq_data_cur[[chromosome]]
      M1_df <- rbind(M1_df, freq_data[freq_data$category == "M1", ])
      M2_df <- rbind(M2_df, freq_data[freq_data$category == "M2", ])
      M3_df <- rbind(M3_df, freq_data[freq_data$category == "M3", ])
      M4_df <- rbind(M4_df, freq_data[freq_data$category == "M4", ])
      M5_df <- rbind(M5_df, freq_data[freq_data$category == "M5", ])
    }
    
    aggregate_frequencies <- function(df) {
      df %>%
        group_by(result) %>%
        summarise(total_frequency = sum(frequency), .groups = "drop") %>%
        arrange(result)
    }
    
    M1_aggregated <- aggregate_frequencies(M1_df)
    M2_aggregated <- aggregate_frequencies(M2_df)
    M3_aggregated <- aggregate_frequencies(M3_df)
    M4_aggregated <- aggregate_frequencies(M4_df)
    M5_aggregated <- aggregate_frequencies(M5_df)
    
    aggregated_data <- list(
      M1 = M1_aggregated,
      M2 = M2_aggregated,
      M3 = M3_aggregated,
      M4 = M4_aggregated,
      M5 = M5_aggregated
    )
    all_aggregated_data[[i]] <- aggregated_data
  }
  
  for (i in seq_along(all_aggregated_data)) {
    aggregated_data <- all_aggregated_data[[i]]
    
    aggregated_data$M1 <- aggregated_data$M1 %>% filter(total_frequency > 0)
    aggregated_data$M2 <- aggregated_data$M2 %>% filter(total_frequency > 0)
    aggregated_data$M3 <- aggregated_data$M3 %>% filter(total_frequency > 0)
    aggregated_data$M4 <- aggregated_data$M4 %>% filter(total_frequency > 0)
    aggregated_data$M5 <- aggregated_data$M5 %>% filter(total_frequency > 0)
    
    all_aggregated_data[[i]] <- aggregated_data
  }
  
  return(all_aggregated_data)
}
################################################################################
# data and files that are needed to calculate the index
directory <- 'C:/Users/julia/Documents/Kamil/busco_fulltable/all_tol'
new_colnames <- c("busco", "chrQ", "Qstart", "Qend")
# list with busco_files and chrom_files
# for this execute below
# busco_files 
# chrom_files

# ALG_data_ALG
# ALGs of node 3 from tree with 289 genomes ( -m 89 -a quick -r 2)
directory2 <- 'C:/Users/julia/Documents/Kamil/muller'
file_path <- file.path(directory2, "minus26_89_n3.tsv")
ALG_data <- read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(ALG_data) <- new_colnames
ALG_data_ALG <- ALG_data

# match names accessions and family
names_acc_fam_2 <- read.table("C:/Users/julia/Documents/Kamil/busco_fulltable/names_acc_fam_2.csv", sep = ",", header = TRUE)
colnames(names_acc_fam_2) <- c("names", "num", "acc", "fam")
head(names_acc_fam_2)

# below are busco_files and chrom_files
################################################################################
busco_files <-c("GCA_963675445.tsv",
                "GCA_951800035.tsv",
                "GCA_956483635.tsv",
                "GCA_947311025.tsv",
                "GCA_963402855.tsv",
                "GCA_956483585.tsv",
                "GCA_949987645.tsv",
                "GCA_905220375.tsv",
                "GCA_944452675.tsv",
                "GCA_936439885.tsv",
                "GCA_963924685.tsv",
                "GCA_963681545.tsv",
                "GCA_932526305.tsv",
                "GCA_914767995.tsv",
                "GCA_963662145.tsv",
                "GCA_937654795.tsv",
                "GCA_916610165.tsv",
                "GCA_949628195.tsv",
                "GCA_963932375.tsv",
                "GCA_947397855.tsv",
                "GCA_963931875.tsv",
                "GCA_930367215.tsv",
                "GCA_949318255.tsv",
                "GCA_943735925.tsv",
                "GCF_958450345.tsv",
                "GCA_942486065.tsv",
                "GCA_950370525.tsv",
                "GCA_916048285.tsv",
                "GCF_022045245.tsv",
                "GCA_963583895.tsv",
                "GCA_941918925.tsv",
                "GCA_933228675.tsv",
                "GCA_951394005.tsv",
                "GCA_932274085.tsv",
                "GCA_958299815.tsv",
                "GCA_932273835.tsv",
                "GCA_936449025.tsv",
                "GCA_930367235.tsv",
                "GCA_927399465.tsv",
                "GCA_014635995.tsv",
                "GCA_963576795.tsv",
                "GCA_963682025.tsv",
                "GCA_040938175.tsv",
                "GCA_963921195.tsv",
                "GCA_949987735.tsv",
                "GCA_949748255.tsv",
                "GCA_021234595.tsv",
                "GCF_963082655.tsv",
                "GCA_963930755.tsv",
                "GCA_947397865.tsv",
                "GCA_963513945.tsv",
                "GCA_958296145.tsv",
                "GCA_964034835.tsv",
                "GCA_963931995.tsv",
                "GCA_963978525.tsv",
                "GCA_949710015.tsv",
                "GCA_963971445.tsv",
                "GCF_016746395.tsv",
                "GCF_004382195.tsv",
                "GCF_004382145.tsv",
                "GCF_000001215.tsv",
                "GCF_016746365.tsv",
                "GCF_016746245.tsv",
                "GCF_016746235.tsv",
                "GCF_037355615.tsv",
                "GCF_014743375.tsv",
                "GCF_025231255.tsv",
                "GCA_014170255.tsv",
                "GCF_017639315.tsv",
                "GCF_009870125.tsv",
                "GCF_003369915.tsv",
                "GCA_008121275.tsv",
                #"GCA_037356375.tsv",
                #"GCA_008121215.tsv",
                "GCA_963969585.tsv",
                "GCF_008121235.tsv",
                "GCA_040893045.tsv",
                "GCA_009664405.tsv",
                "GCF_023558535.tsv",
                "GCF_009650485.tsv",
                "GCA_023558455.tsv",
                "GCA_023558445.tsv",
                "GCF_023558435.tsv",
                "GCA_963583835.tsv",
                "GCA_951394115.tsv",
                "GCF_004354385.tsv",
                "GCA_963924055.tsv",
                "GCA_958299025.tsv",
                "GCA_949708635.tsv",
                "GCA_030788265.tsv",
                "GCF_011750605.tsv",
                "GCA_949987675.tsv",
                "GCA_034638295.tsv",
                "GCA_020466095.tsv",
                "GCA_002237135.tsv",
                "GCF_024586455.tsv",
                "GCF_016617805.tsv",
                "GCF_023373825.tsv",
                "GCA_027475135.tsv",
                "GCF_028554725.tsv",
                "GCA_031772095.tsv",
                "GCA_030068015.tsv",
                "GCF_028408465.tsv",
                "GCF_027943255.tsv",
                "GCA_963971545.tsv",
                "GCA_951828415.tsv",
                "GCA_951509765.tsv",
                "GCA_949987695.tsv",
                "GCA_964058745.tsv",
                "GCA_040869045.tsv",
                "GCA_960531455.tsv",
                "GCA_964056695.tsv",
                "GCA_963668005.tsv",
                "GCA_963971225.tsv",
                "GCA_963691655.tsv",
                "GCA_958295655.tsv",
                "GCA_964194425.tsv",
                "GCA_963854835.tsv",
                "GCA_951394025.tsv",
                "GCA_949127995.tsv",
                "GCA_963457695.tsv",
                "GCA_920105625.tsv",
                "GCA_949629155.tsv",
                "GCA_943737955.tsv",
                "GCA_937620795.tsv",
                "GCA_963966595.tsv",
                "GCA_949752815.tsv",
                "GCA_922984085.tsv",
                "GCA_963920525.tsv",
                "GCA_914767935.tsv",
                "GCA_943735915.tsv",
                "GCA_958299015.tsv",
                "GCA_951804985.tsv",
                "GCA_947389925.tsv",
                "GCA_017309665.tsv",
                "GCA_963680825.tsv",
                "GCA_963662105.tsv",
                "GCA_949775025.tsv",
                "GCA_951812925.tsv",
                "GCA_947095585.tsv",
                "GCA_951509635.tsv",
                "GCA_920104205.tsv",
                "GCF_945859685.tsv",
                "GCA_905146935.tsv",
                "GCA_958431115.tsv",
                "GCA_951394125.tsv",
                "GCA_929447395.tsv",
                "GCA_947049315.tsv",
                "GCA_932273885.tsv",
                "GCA_949320155.tsv",
                "GCA_963921185.tsv",
                "GCA_946251815.tsv",
                "GCA_963931855.tsv",
                "GCA_911387755.tsv",
                "GCA_910595825.tsv",
                "GCF_945859705.tsv",
                "GCA_948107695.tsv",
                "GCA_943590905.tsv",
                "GCA_963971375.tsv",
                "GCA_916050605.tsv",
                "GCA_951217065.tsv",
                "GCA_949752695.tsv",
                "GCA_914767635.tsv",
                "GCA_947623515.tsv",
                "GCA_963082955.tsv",
                "GCA_949126925.tsv",
                "GCA_946477585.tsv",
                "GCA_951230905.tsv",
                "GCA_948293265.tsv",
                "GCA_916610125.tsv",
                "GCA_936431705.tsv",
                "GCA_955612985.tsv",
                "GCA_949372485.tsv",
                "GCA_963576555.tsv",
                "GCA_949715645.tsv",
                "GCA_949824845.tsv",
                "GCA_932526625.tsv",
                "GCA_949716465.tsv",
                "GCA_945910035.tsv",
                "GCA_963583995.tsv",
                "GCA_905220385.tsv",
                "GCA_963082735.tsv",
                "GCA_946965025.tsv",
                "GCA_951813785.tsv",
                "GCA_917880715.tsv",
                "GCA_949775065.tsv",
                "GCA_963966105.tsv",
                "GCA_949129095.tsv",
                "GCA_907269105.tsv",
                "GCA_928272305.tsv",
                "GCA_963555765.tsv",
                "GCA_916610145.tsv",
                "GCA_907269125.tsv",
                "GCA_964034865.tsv",
                "GCA_905231855.tsv",
                "GCA_930367185.tsv",
                "GCA_955652365.tsv",
                "GCA_944738645.tsv",
                "GCA_905187475.tsv",
                "GCA_958301585.tsv",
                "GCA_951905685.tsv",
                "GCA_963931805.tsv",
                "GCA_963082915.tsv",
                "GCA_947095535.tsv",
                "GCA_964204755.tsv",
                "GCA_963942445.tsv",
                "GCA_964204855.tsv",
                "GCA_964007475.tsv",
                "GCA_963932195.tsv",
                "GCA_949752835.tsv",
                "GCA_963920725.tsv",
                "GCA_963971155.tsv",
                "GCA_951394055.tsv",
                "GCA_939192795.tsv",
                "GCA_932526495.tsv",
                "GCA_963920795.tsv",
                "GCA_951509405.tsv",
                "GCA_959613345.tsv",
                "GCA_933228815.tsv",
                "GCA_963969385.tsv",
                "GCA_949715965.tsv",
                "GCA_947538895.tsv",
                "GCA_963930735.tsv",
                "GCA_963924295.tsv",
                "GCA_964036015.tsv",
                "GCA_963082835.tsv",
                "GCA_963855945.tsv",
                "GCA_949987705.tsv",
                "GCA_963969285.tsv",
                "GCA_954870665.tsv",
                "GCA_963942545.tsv",
                "GCA_963978885.tsv",
                "GCA_947369275.tsv",
                "GCA_964017055.tsv",
                "GCA_949715475.tsv",
                "GCF_905115235.tsv",
                "GCA_951812415.tsv",
                "GCA_949128065.tsv",
                "GCA_963669355.tsv",
                "GCA_963422695.tsv",
                "GCA_964211845.tsv",
                "GCA_964034945.tsv",
                "GCA_963971475.tsv",
                "GCA_958298945.tsv",
                "GCA_947359355.tsv",
                "GCA_030463065.tsv",
                "GCA_028476595.tsv",
                "GCA_029228625.tsv",
                "GCA_014529535.tsv",
                "GCA_958336335.tsv",
                "GCA_910594885.tsv",
                "GCA_963854165.tsv",
                "GCA_963556165.tsv",
                "GCA_958295665.tsv",
                "GCA_951394425.tsv",
                "GCA_947310385.tsv",
                "GCA_932526605.tsv",
                "GCA_963691935.tsv",
                "GCA_963932295.tsv",
                "GCA_963556175.tsv",
                "GCA_963693315.tsv",
                "GCA_963082725.tsv",
                "GCA_963942525.tsv",
                "GCA_961205885.tsv",
                "GCF_943734845.tsv",
                "GCA_016170035.tsv",
                "GCA_016254315.tsv",
                "GCA_016170015.tsv",
                "GCF_943734755.tsv",
                "GCF_943734725.tsv",
                "GCF_943734695.tsv",
                "GCF_013141755.tsv",
                "GCF_943734735.tsv",
                "GCF_943734685.tsv",
                "GCF_016920715.tsv",
                "GCF_017562075.tsv",
                "GCF_943737925.tsv",
                "GCF_943734765.tsv",
                "GCF_943734705.tsv",
                "GCF_013758885.tsv",
                "GCA_943734665.tsv",
                "GCF_943734745.tsv",
                "GCF_943735745.tsv",
                "GCF_943734635.tsv",
                "GCF_030247195.tsv",
                "GCF_030247185.tsv",
                "GCF_943734655.tsv",
                "GCF_029784165.tsv",
                "GCF_029784135.tsv",
                "GCF_035046485.tsv",
                "GCF_002204515.tsv",
                "GCF_024139115.tsv",
                "GCA_040801935.tsv",
                "GCF_016801865.tsv",
                "GCF_015732765.tsv",
                "GCA_964187845.tsv",
                "GCF_029784155.tsv",
                "GCA_963573255.tsv",
                "GCA_917627325.tsv",
                "GCA_026123125.tsv",
                "GCA_018290095.tsv",
                "GCA_033064975.tsv",
                "GCA_033063855.tsv",
                "GCA_018397935.tsv",
                "GCF_036172545.tsv",
                "GCF_024763615.tsv",
                "GCF_024334085.tsv")
################################################################################
chrom_files <-c("GCA_963675445_chrominf.tsv",
                "GCA_951800035_chrominf.tsv",
                "GCA_956483635_chrominf.tsv",
                "GCA_947311025_chrominf.tsv",
                "GCA_963402855_chrominf.tsv",
                "GCA_956483585_chrominf.tsv",
                "GCA_949987645_chrominf.tsv",
                "GCA_905220375_chrominf.tsv",
                "GCA_944452675_chrominf.tsv",
                "GCA_936439885_chrominf.tsv",
                "GCA_963924685_chrominf.tsv",
                "GCA_963681545_chrominf.tsv",
                "GCA_932526305_chrominf.tsv",
                "GCA_914767995_chrominf.tsv",
                "GCA_963662145_chrominf.tsv",
                "GCA_937654795_chrominf.tsv",
                "GCA_916610165_chrominf.tsv",
                "GCA_949628195_chrominf.tsv",
                "GCA_963932375_chrominf.tsv",
                "GCA_947397855_chrominf.tsv",
                "GCA_963931875_chrominf.tsv",
                "GCA_930367215_chrominf.tsv",
                "GCA_949318255_chrominf.tsv",
                "GCA_943735925_chrominf.tsv",
                "GCF_958450345_chrominf.tsv",
                "GCA_942486065_chrominf.tsv",
                "GCA_950370525_chrominf.tsv",
                "GCA_916048285_chrominf.tsv",
                "GCF_022045245_chrominf.tsv",
                "GCA_963583895_chrominf.tsv",
                "GCA_941918925_chrominf.tsv",
                "GCA_933228675_chrominf.tsv",
                "GCA_951394005_chrominf.tsv",
                "GCA_932274085_chrominf.tsv",
                "GCA_958299815_chrominf.tsv",
                "GCA_932273835_chrominf.tsv",
                "GCA_936449025_chrominf.tsv",
                "GCA_930367235_chrominf.tsv",
                "GCA_927399465_chrominf.tsv",
                "GCA_014635995_chrominf.tsv",
                "GCA_963576795_chrominf.tsv",
                "GCA_963682025_chrominf.tsv",
                "GCA_040938175_chrominf.tsv",
                "GCA_963921195_chrominf.tsv",
                "GCA_949987735_chrominf.tsv",
                "GCA_949748255_chrominf.tsv",
                "GCA_021234595_chrominf.tsv",
                "GCF_963082655_chrominf.tsv",
                "GCA_963930755_chrominf.tsv",
                "GCA_947397865_chrominf.tsv",
                "GCA_963513945_chrominf.tsv",
                "GCA_958296145_chrominf.tsv",
                "GCA_964034835_chrominf.tsv",
                "GCA_963931995_chrominf.tsv",
                "GCA_963978525_chrominf.tsv",
                "GCA_949710015_chrominf.tsv",
                "GCA_963971445_chrominf.tsv",
                "GCF_016746395_chrominf.tsv",
                "GCF_004382195_chrominf.tsv",
                "GCF_004382145_chrominf.tsv",
                "GCF_000001215_chrominf.tsv",
                "GCF_016746365_chrominf.tsv",
                "GCF_016746245_chrominf.tsv",
                "GCF_016746235_chrominf.tsv",
                "GCF_037355615_chrominf.tsv",
                "GCF_014743375_chrominf.tsv",
                "GCF_025231255_chrominf.tsv",
                "GCA_014170255_chrominf.tsv",
                "GCF_017639315_chrominf.tsv",
                "GCF_009870125_chrominf.tsv",
                "GCF_003369915_chrominf.tsv",
                "GCA_008121275_chrominf.tsv",
                #"GCA_037356375_chrominf.tsv",
                #"GCA_008121215_chrominf.tsv",
                "GCA_963969585_chrominf.tsv",
                "GCF_008121235_chrominf.tsv",
                "GCA_040893045_chrominf.tsv",
                "GCA_009664405_chrominf.tsv",
                "GCF_023558535_chrominf.tsv",
                "GCF_009650485_chrominf.tsv",
                "GCA_023558455_chrominf.tsv",
                "GCA_023558445_chrominf.tsv",
                "GCF_023558435_chrominf.tsv",
                "GCA_963583835_chrominf.tsv",
                "GCA_951394115_chrominf.tsv",
                "GCF_004354385_chrominf.tsv",
                "GCA_963924055_chrominf.tsv",
                "GCA_958299025_chrominf.tsv",
                "GCA_949708635_chrominf.tsv",
                "GCA_030788265_chrominf.tsv",
                "GCF_011750605_chrominf.tsv",
                "GCA_949987675_chrominf.tsv",
                "GCA_034638295_chrominf.tsv",
                "GCA_020466095_chrominf.tsv",
                "GCA_002237135_chrominf.tsv",
                "GCF_024586455_chrominf.tsv",
                "GCF_016617805_chrominf.tsv",
                "GCF_023373825_chrominf.tsv",
                "GCA_027475135_chrominf.tsv",
                "GCF_028554725_chrominf.tsv",
                "GCA_031772095_chrominf.tsv",
                "GCA_030068015_chrominf.tsv",
                "GCF_028408465_chrominf.tsv",
                "GCF_027943255_chrominf.tsv",
                "GCA_963971545_chrominf.tsv",
                "GCA_951828415_chrominf.tsv",
                "GCA_951509765_chrominf.tsv",
                "GCA_949987695_chrominf.tsv",
                "GCA_964058745_chrominf.tsv",
                "GCA_040869045_chrominf.tsv",
                "GCA_960531455_chrominf.tsv",
                "GCA_964056695_chrominf.tsv",
                "GCA_963668005_chrominf.tsv",
                "GCA_963971225_chrominf.tsv",
                "GCA_963691655_chrominf.tsv",
                "GCA_958295655_chrominf.tsv",
                "GCA_964194425_chrominf.tsv",
                "GCA_963854835_chrominf.tsv",
                "GCA_951394025_chrominf.tsv",
                "GCA_949127995_chrominf.tsv",
                "GCA_963457695_chrominf.tsv",
                "GCA_920105625_chrominf.tsv",
                "GCA_949629155_chrominf.tsv",
                "GCA_943737955_chrominf.tsv",
                "GCA_937620795_chrominf.tsv",
                "GCA_963966595_chrominf.tsv",
                "GCA_949752815_chrominf.tsv",
                "GCA_922984085_chrominf.tsv",
                "GCA_963920525_chrominf.tsv",
                "GCA_914767935_chrominf.tsv",
                "GCA_943735915_chrominf.tsv",
                "GCA_958299015_chrominf.tsv",
                "GCA_951804985_chrominf.tsv",
                "GCA_947389925_chrominf.tsv",
                "GCA_017309665_chrominf.tsv",
                "GCA_963680825_chrominf.tsv",
                "GCA_963662105_chrominf.tsv",
                "GCA_949775025_chrominf.tsv",
                "GCA_951812925_chrominf.tsv",
                "GCA_947095585_chrominf.tsv",
                "GCA_951509635_chrominf.tsv",
                "GCA_920104205_chrominf.tsv",
                "GCF_945859685_chrominf.tsv",
                "GCA_905146935_chrominf.tsv",
                "GCA_958431115_chrominf.tsv",
                "GCA_951394125_chrominf.tsv",
                "GCA_929447395_chrominf.tsv",
                "GCA_947049315_chrominf.tsv",
                "GCA_932273885_chrominf.tsv",
                "GCA_949320155_chrominf.tsv",
                "GCA_963921185_chrominf.tsv",
                "GCA_946251815_chrominf.tsv",
                "GCA_963931855_chrominf.tsv",
                "GCA_911387755_chrominf.tsv",
                "GCA_910595825_chrominf.tsv",
                "GCF_945859705_chrominf.tsv",
                "GCA_948107695_chrominf.tsv",
                "GCA_943590905_chrominf.tsv",
                "GCA_963971375_chrominf.tsv",
                "GCA_916050605_chrominf.tsv",
                "GCA_951217065_chrominf.tsv",
                "GCA_949752695_chrominf.tsv",
                "GCA_914767635_chrominf.tsv",
                "GCA_947623515_chrominf.tsv",
                "GCA_963082955_chrominf.tsv",
                "GCA_949126925_chrominf.tsv",
                "GCA_946477585_chrominf.tsv",
                "GCA_951230905_chrominf.tsv",
                "GCA_948293265_chrominf.tsv",
                "GCA_916610125_chrominf.tsv",
                "GCA_936431705_chrominf.tsv",
                "GCA_955612985_chrominf.tsv",
                "GCA_949372485_chrominf.tsv",
                "GCA_963576555_chrominf.tsv",
                "GCA_949715645_chrominf.tsv",
                "GCA_949824845_chrominf.tsv",
                "GCA_932526625_chrominf.tsv",
                "GCA_949716465_chrominf.tsv",
                "GCA_945910035_chrominf.tsv",
                "GCA_963583995_chrominf.tsv",
                "GCA_905220385_chrominf.tsv",
                "GCA_963082735_chrominf.tsv",
                "GCA_946965025_chrominf.tsv",
                "GCA_951813785_chrominf.tsv",
                "GCA_917880715_chrominf.tsv",
                "GCA_949775065_chrominf.tsv",
                "GCA_963966105_chrominf.tsv",
                "GCA_949129095_chrominf.tsv",
                "GCA_907269105_chrominf.tsv",
                "GCA_928272305_chrominf.tsv",
                "GCA_963555765_chrominf.tsv",
                "GCA_916610145_chrominf.tsv",
                "GCA_907269125_chrominf.tsv",
                "GCA_964034865_chrominf.tsv",
                "GCA_905231855_chrominf.tsv",
                "GCA_930367185_chrominf.tsv",
                "GCA_955652365_chrominf.tsv",
                "GCA_944738645_chrominf.tsv",
                "GCA_905187475_chrominf.tsv",
                "GCA_958301585_chrominf.tsv",
                "GCA_951905685_chrominf.tsv",
                "GCA_963931805_chrominf.tsv",
                "GCA_963082915_chrominf.tsv",
                "GCA_947095535_chrominf.tsv",
                "GCA_964204755_chrominf.tsv",
                "GCA_963942445_chrominf.tsv",
                "GCA_964204855_chrominf.tsv",
                "GCA_964007475_chrominf.tsv",
                "GCA_963932195_chrominf.tsv",
                "GCA_949752835_chrominf.tsv",
                "GCA_963920725_chrominf.tsv",
                "GCA_963971155_chrominf.tsv",
                "GCA_951394055_chrominf.tsv",
                "GCA_939192795_chrominf.tsv",
                "GCA_932526495_chrominf.tsv",
                "GCA_963920795_chrominf.tsv",
                "GCA_951509405_chrominf.tsv",
                "GCA_959613345_chrominf.tsv",
                "GCA_933228815_chrominf.tsv",
                "GCA_963969385_chrominf.tsv",
                "GCA_949715965_chrominf.tsv",
                "GCA_947538895_chrominf.tsv",
                "GCA_963930735_chrominf.tsv",
                "GCA_963924295_chrominf.tsv",
                "GCA_964036015_chrominf.tsv",
                "GCA_963082835_chrominf.tsv",
                "GCA_963855945_chrominf.tsv",
                "GCA_949987705_chrominf.tsv",
                "GCA_963969285_chrominf.tsv",
                "GCA_954870665_chrominf.tsv",
                "GCA_963942545_chrominf.tsv",
                "GCA_963978885_chrominf.tsv",
                "GCA_947369275_chrominf.tsv",
                "GCA_964017055_chrominf.tsv",
                "GCA_949715475_chrominf.tsv",
                "GCF_905115235_chrominf.tsv",
                "GCA_951812415_chrominf.tsv",
                "GCA_949128065_chrominf.tsv",
                "GCA_963669355_chrominf.tsv",
                "GCA_963422695_chrominf.tsv",
                "GCA_964211845_chrominf.tsv",
                "GCA_964034945_chrominf.tsv",
                "GCA_963971475_chrominf.tsv",
                "GCA_958298945_chrominf.tsv",
                "GCA_947359355_chrominf.tsv",
                "GCA_030463065_chrominf.tsv",
                "GCA_028476595_chrominf.tsv",
                "GCA_029228625_chrominf.tsv",
                "GCA_014529535_chrominf.tsv",
                "GCA_958336335_chrominf.tsv",
                "GCA_910594885_chrominf.tsv",
                "GCA_963854165_chrominf.tsv",
                "GCA_963556165_chrominf.tsv",
                "GCA_958295665_chrominf.tsv",
                "GCA_951394425_chrominf.tsv",
                "GCA_947310385_chrominf.tsv",
                "GCA_932526605_chrominf.tsv",
                "GCA_963691935_chrominf.tsv",
                "GCA_963932295_chrominf.tsv",
                "GCA_963556175_chrominf.tsv",
                "GCA_963693315_chrominf.tsv",
                "GCA_963082725_chrominf.tsv",
                "GCA_963942525_chrominf.tsv",
                "GCA_961205885_chrominf.tsv",
                "GCF_943734845_chrominf.tsv",
                "GCA_016170035_chrominf.tsv",
                "GCA_016254315_chrominf.tsv",
                "GCA_016170015_chrominf.tsv",
                "GCF_943734755_chrominf.tsv",
                "GCF_943734725_chrominf.tsv",
                "GCF_943734695_chrominf.tsv",
                "GCF_013141755_chrominf.tsv",
                "GCF_943734735_chrominf.tsv",
                "GCF_943734685_chrominf.tsv",
                "GCF_016920715_chrominf.tsv",
                "GCF_017562075_chrominf.tsv",
                "GCF_943737925_chrominf.tsv",
                "GCF_943734765_chrominf.tsv",
                "GCF_943734705_chrominf.tsv",
                "GCF_013758885_chrominf.tsv",
                "GCA_943734665_chrominf.tsv",
                "GCF_943734745_chrominf.tsv",
                "GCF_943735745_chrominf.tsv",
                "GCF_943734635_chrominf.tsv",
                "GCF_030247195_chrominf.tsv",
                "GCF_030247185_chrominf.tsv",
                "GCF_943734655_chrominf.tsv",
                "GCF_029784165_chrominf.tsv",
                "GCF_029784135_chrominf.tsv",
                "GCF_035046485_chrominf.tsv",
                "GCF_002204515_chrominf.tsv",
                "GCF_024139115_chrominf.tsv",
                "GCA_040801935_chrominf.tsv",
                "GCF_016801865_chrominf.tsv",
                "GCF_015732765_chrominf.tsv",
                "GCA_964187845_chrominf.tsv",
                "GCF_029784155_chrominf.tsv",
                "GCA_963573255_chrominf.tsv",
                "GCA_917627325_chrominf.tsv",
                "GCA_026123125_chrominf.tsv",
                "GCA_018290095_chrominf.tsv",
                "GCA_033064975_chrominf.tsv",
                "GCA_033063855_chrominf.tsv",
                "GCA_018397935_chrominf.tsv",
                "GCF_036172545_chrominf.tsv",
                "GCF_024763615_chrominf.tsv",
                "GCF_024334085_chrominf.tsv")
################################################################################
all_aggregated_data_ALG_2 <- process_busco_data(directory, busco_files, chrom_files, ALG_data_ALG, new_colnames)
################################################################################
chrQ_frequencies <- c(M1 = 799, M2 = 447, M3 = 549, M4 = 668, M5 = 502)
################################################################################
# this is important to correctly map names of species!
# make busco_list to use for mapping below
busco_files_clean <- sub("\\.tsv$", "", busco_files)
busco_list <- setNames(vector("list", length(busco_files_clean)), busco_files_clean)
################################################################################
calculate_conservedness <- function(data, S_total) {
  weighted_sum <- sum(data$total_frequency * log(data$result + 1))
  conservedness_index <- weighted_sum / S_total
  conservedness_index <- 1 - conservedness_index
  return(conservedness_index)
}
################################################################################
results <- list()
for (genome_idx in seq_along(all_aggregated_data_ALG_2)) {
  genome_data <- all_aggregated_data_ALG_2[[genome_idx]]
  genome_results <- list()
  
  for (alg_idx in seq_along(genome_data)) {
    alg_data <- genome_data[[alg_idx]]
    S_total <- chrQ_frequencies[names(genome_data)[alg_idx]]
    conservedness <- calculate_conservedness(alg_data, S_total)
    genome_results[[names(genome_data)[alg_idx]]] <- conservedness
  }
  
  genome_name <- sprintf("%03d", genome_idx)
  results[[genome_name]] <- genome_results
}

# Convert results to a tidy format
tidy_results <- tibble(
  Genome = rep(names(results), sapply(results, length)),
  ALG = unlist(lapply(results, names)),
  Conservedness = unlist(results)
)

mapping <- setNames(names(busco_list), sprintf("%03d", seq_along(busco_list)))

tidy_results <- tidy_results %>%
  mutate(Names = mapping[Genome])
names <- read.table("C:/Users/julia/Documents/Kamil/trees/acc_names_match.txt", sep = ",")
tidy_results <- tidy_results %>%
  left_join(names, by = c("Names" = "V1")) %>%
  rename(Species = V2)

tidy_results$Genome <- factor(tidy_results$Genome, levels = unique(tidy_results$Genome))

species_labels <- tidy_results %>%
  distinct(Genome, Species) %>%
  arrange(Genome) %>%
  pull(Species)

tidy_results <- tidy_results %>%
  left_join(names_acc_fam_2, by = c("Species" = "names")) %>%
  select(Genome, ALG, Conservedness, Names, Species, fam, num)

tidy_results <- tidy_results %>%
  rename(Family = fam, Num = num)

family_summary <- tidy_results %>%
  group_by(Family, ALG) %>%
  summarise(
    mean_conservedness = mean(Conservedness, na.rm = TRUE),
    sd_conservedness = ifelse(n() > 1, sd(Conservedness, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  left_join(
    tidy_results %>% select(Family, Num) %>% distinct(),
    by = "Family"
  ) %>%
  arrange(Num, ALG)

family_summary <- family_summary %>%
  group_by(Family) %>%
  filter(Num == max(Num)) %>%
  ungroup()

unique_families <- family_summary %>%
  select(Family, Num) %>%
  distinct() %>%
  arrange(-Num)
unique_families

family_summary$Family <- factor(family_summary$Family, levels = unique_families$Family)

family_summary$Num <- families_order_tree$y[match(family_summary$Family, families_order_tree$label)]
family_summary_sorted <- family_summary %>%
  arrange(desc(Num))
print(family_summary_sorted)

family_summary_sorted$Family <- factor(family_summary_sorted$Family, 
                                       levels = unique(family_summary_sorted$Family))

#col_list <- c("yellow4", "springgreen3", "#F3cdC7", "steelblue2", "#FF7B9C")
col_list <- c("#1573afff", "#e59d38ff", "#f0e354ff", "#169e73ff", "#60b5e1ff")

plot <- ggplot(family_summary_sorted, aes(x = Family, y = mean_conservedness, color = ALG)) +
  geom_point(size = 3.3, position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = mean_conservedness - sd_conservedness, ymax = mean_conservedness + sd_conservedness),
    width = 0.2, position = position_dodge(width = 0.5)
  ) + # Fehlerbalken
  scale_color_manual(values = col_list, name = "Muller e.") +
  labs(
    x = NULL,
    y = NULL,
    color = "ALG"
  ) +
  ylim(0.4, 1.0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1, size = 13
    ),
    axis.text.y = element_text(size = 17),
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    #panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
plot

ggsave("family_conservedness_ALG_new_colour.png", plot = plot, width = 12, height = 6, dpi = 500)
################################################################################

