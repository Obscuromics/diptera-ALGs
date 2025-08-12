require('gsheet')

all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?usp=sharing", format='csv'),
                        stringsAsFactors = F, header = T, check.names = F)


print("Keeping only genomes that have TO ADD labeled as KEEP")
genome_in_analysis <- all_genome_data[all_genome_data[, 'TO ADD'] == 'KEEP', ]

if (any(duplicated(genome_in_analysis[, 'taxon_id']))){ # all species are unique at least by name
    stop("There is a species there twice (checked by TaxonID)")
}

# removing outgroup
genome_in_analysis <- genome_in_analysis[genome_in_analysis[, 'species'] != 'Panorpa_germanica', ]

print(paste("Dipteran species:", nrow(genome_in_analysis)))
print(paste("Dipteran families:", length(table(genome_in_analysis[, 'family']))))

# head(genome_in_analysis)
chromsome_numbers <- genome_in_analysis[,'n_chromosomes_in_fasta-y']
min_chrom <- min(chromsome_numbers, na.rm = T) # 3
max_chrom <- max(chromsome_numbers, na.rm = T) # 10

# families with 3 chromsomes
# consideraton: is the 'num_chrom' what we want?

print(paste("Minimum chromosomes:", min_chrom, " in ", genome_in_analysis[which(chromsome_numbers == min_chrom), 'family']))
print(paste("Minimum chromosomes:", max_chrom, " in ", genome_in_analysis[which(chromsome_numbers == max_chrom), 'family']))

# hist(genome_in_analysis[,'genome_size'] / 1e6, breaks = 60) 

passing_QC_genomes <- genome_in_analysis[genome_in_analysis[, "QC"] == "", ]

# barplot(table(chromsome_numbers))
# removing the two with wrong numbers, probably should more than that
print(paste("Dipteran species with 6 chromosomes assembled:", sum(passing_QC_genomes[,'n_chromosomes_in_fasta-y'] == 6) / 338))


smallest_genome_index <- which.min(passing_QC_genomes[,'genome_size']) #, na.rm = T
largest_genome_index <- which.max(passing_QC_genomes[,'genome_size']) #, na.rm = T

print(paste("Smallest genome: ", round(passing_QC_genomes[smallest_genome_index, 'genome_size'] / 1e6, 2), "Mbp of ", genome_in_analysis[smallest_genome_index, 'species'], ';', genome_in_analysis[smallest_genome_index, 'family']))
# 31037916 ### BS
print(paste("Largest genome: ", round(passing_QC_genomes[largest_genome_index, 'genome_size'] / 1e6, 2), "Mbp of ", genome_in_analysis[largest_genome_index, 'species'], ';', genome_in_analysis[largest_genome_index, 'family']))

# head(passing_QC_genomes[order(passing_QC_genomes[,'genome_size']), ])
# tail(passing_QC_genomes[order(passing_QC_genomes[,'genome_size']), ])