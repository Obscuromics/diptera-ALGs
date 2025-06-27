require('gsheet')

all_genome_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?usp=sharing", format='csv'),
                        stringsAsFactors = F, header = T, check.names = F)

# subsetting the genome data to those that were used for the analysis
genome_data <- all_genome_data[all_genome_data[, 'TO ADD'] == 'KEEP', ]

# THIS WILL AHVE TO BE REPLACED WITH NCK ORDR OF SPECIES
species_family_tab <- genome_data[, c('species', 'family')]

number_of_species <- nrow(species_family_tab)
families <- unique(species_family_tab[, 'family'])
number_of_families <- length(families)

species_y <- 1:number_of_species / number_of_species
names(species_y) <- species_family_tab[, 'species']

family_y <- 1:number_of_families / number_of_families
names(family_y) <- families

# code from https://stackoverflow.com/a/32100956/29281675
curveMaker <- function(x1, y1, x2, y2, ...){
    return(curve( plogis( x, scale = 0.08, loc = (x1 + x2) /2 ) * (y2-y1) + y1, 
                   x1, x2, add = TRUE, ...))
}

pdf('figures/figure_1_sp_family_connectors.pdf', width = 4, height = 10)
# blank canvas
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = F, xlab = '', ylab = '')

text(0.05, species_y, species_family_tab[, 'species'], cex = 0.15, pos = 2)
text(0.60, family_y, families, cex = 0.70, pos = 4)

# polygon(c(0.1, 0.5, 0.9, 0.5), c(0.5, 0.9, 0.5, 0.1), col = 'red')
x1 <- 0.05
x2 <- 0.6
for ( i in 1:length(families)){
    fam <- families[i]
    col <- ifelse(i %% 2 == 0, 'gray60', 'gray')
    all_species <- species_family_tab[species_family_tab[, 'family'] == fam, 'species']
    ys1 <- max(species_y[all_species])
    ys2 <- min(species_y[all_species])
    yf <- family_y[fam]
    top_curve <- curveMaker(x1, ys1, x2, yf)
    bot_curve <- curveMaker(x1, ys2, x2, yf)
    polygon(c(top_curve[['x']], rev(bot_curve[['x']])), c(top_curve[['y']], rev(bot_curve[['y']])), col = col)
    # curveMaker(0.05, species_y[species_family_tab[i, 'species']], 0.6, family_y[species_family_tab[i, 'family']])
}

dev.off()

