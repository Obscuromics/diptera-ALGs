require('gsheet')
require('dplyr')

all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?usp=sharing", format='csv'), stringsAsFactors = F, header = T, check.names = F)
chr_size_and_busco_data <- read.csv(text = gsheet2text('https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?pli=1&gid=1940964825#gid=1940964825', format='csv'), stringsAsFactors = F, header = T, check.names = F)
busco_score_data <- read.csv(text = gsheet2text('https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?pli=1&gid=1742946067#gid=1742946067', format='csv'), stringsAsFactors = F, header = T, check.names = F)
sex_chr_transition_data <- read.csv(text = gsheet2text('https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?gid=295210694#gid=295210694', format='csv'), stringsAsFactors = F, header = T, check.names = F)

# Families where all genomes were excluded from the reconstruction
family_exclusion_tally <- all_genome_data %>%
  select(family, excluded_from_ALG_inference) %>%
  group_by(family, excluded_from_ALG_inference) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n))
family_exclusion_tally  %>% filter(excluded_from_ALG_inference == "yes", proportion == 1)

# Sheet 1
sheet_1_columns = c('taxon_id',  
            'suborder',
            'superfamily',
            'family',
            'TO ADD',  
            'subfamily',
            'species',
            'excluded_from_ALG_inference',
            'exclusion_reason',
            'accession',
            'assembly',
            'source',
            'DOI',
            'X_chrom_NCBI_1',
            'X_chrom_NCBI_2',
            'n_chromosomes_in_fasta',
            'genome_size',
            'alg_6_resembling_chromosome'
)

sheet_1_df <- all_genome_data %>%
  select(any_of(sheet_1_columns)) #%>%
  # filter(TO ADD != EXCLUDE) Could remove no permission entries

write.table(sheet_1_df, "tables/supplementary_table_1.tsv", row.names=FALSE, sep='\t', quote=F)
write.table(chr_size_and_busco_data, "tables/supplementary_table_2.tsv", row.names=FALSE, sep='\t', quote=F)
write.table(busco_score_data, "tables/supplementary_table_3.tsv", row.names=FALSE, sep='\t', quote=F)
write.table(sex_chr_transition_data, "tables/supplementary_table_4.tsv", row.names=FALSE, sep='\t', quote=F)
