require('gsheet')

all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?usp=sharing", format='csv'),
                        stringsAsFactors = F, header = T, check.names = F)



write.table(all_genome_data[, 'accession'], col.names = F, file = 'list_of_accessions.list', row.names = F, quote = F)

command <- "bash ./scripts/figure_1/get_genome_sizes.sh list_of_accessions.list"
#cat(command)
system(command)

file.remove('list_of_accessions.list')