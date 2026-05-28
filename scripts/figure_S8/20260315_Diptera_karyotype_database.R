
#------------------------------------------------------------------------------#
library(fields)
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# csv file downloaded from the diptera karyotype database
diptera_karyotypes <- read.csv("/Users/jg40/Documents/diptera_centromeres/20260119_plot_genomes_chromosomes_genes/diptera-2026-03-15.csv")
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# families descending from ancestral karyotype gamma
families <- c(
  "Asilidae","Therevidae","Empididae","Dolichopodidae","Hybotidae","Lonchopodidiae",
  "Phoridae","Syrphidae","Pipunculidae","Lauxaniidae","Coelopidae", "Conopidae", "Dryomyzidae", "Sciomyzidae",
  "Anthomyzidae", "Heleomyzidae", "Diopsidae", "Clusiidae","Micropezidae","Pallopteridae","Tephritidae","Ulidiidae",
  "Megamerinidae","Opomyzidae","Sepsidae","Chloropidae","Agromyzidae","Drosophilidae",
  "Hippoboscidae","Glossinidae","Muscidae","Scathophagidae","Anthomyiidae",
  "Sarcophagidae","Tachinidae","Pollenidae","Rhiniphoridae","Calliphoridae","Rhiniidae",
  "Chamaemyiidae", "Oestridae", "Piophilidae", "Platystomatidae", "Psilidae", "Richardiidae" # additional families
)
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# only keep gamma families
filtered <- diptera_karyotypes[diptera_karyotypes$Family %in% families, ]
found_families <- intersect(families, unique(diptera_karyotypes$Family))
found_families

missing_families <- setdiff(families, unique(diptera_karyotypes$Family))
missing_families
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# only keep species with annotated karyotype/ centromeres
filtered2 <- filtered[rowSums(!is.na(filtered[,5:10])) > 0, ]
# for 152 species the karyotype/ centromeres were annotated

# only keep those with haploid karyotype 6
filtered3 <- filtered2[filtered2$Haploid.Number == 6, ]
# for 130 out of 152 the haploid chromosome number is 6

# combine sub and actual (e.g. submetacentric and metacentric)
filtered4 <- filtered3
filtered4$linke6  <- rowSums(filtered4[, c(5,6)], na.rm = TRUE)
filtered4$linke8  <- rowSums(filtered4[, c(7,8)], na.rm = TRUE)
filtered4$linke10 <- rowSums(filtered4[, c(9,10)], na.rm = TRUE)
filtered4 <- filtered4[, -c(5:10)]
colnames(filtered4)[colnames(filtered4) == "linke6"]  <- "metacentric"
colnames(filtered4)[colnames(filtered4) == "linke8"]  <- "telocentric"
colnames(filtered4)[colnames(filtered4) == "linke10"] <- "acrocentric"
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# count species with at least 5 metacentric chromosomes out of 6
filtered5 <- filtered4[filtered4$metacentric >= 5, ]
sort(table(filtered5$Family), decreasing = TRUE)
# for 120 species, out of 130 with HaploidNumber 6, there are 120 with at least 5 metacentric chromosomes
# "Calliphoridae" "Lauxaniidae"   "Sarcophagidae" "Sciomyzidae"   "Tachinidae"    "Tephritidae"   "Ulidiidae"
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# investaige the remaining 32 species
# not-haploid 6, at least one value in column 5:10
non6 <- filtered2[filtered2$Haploid.Number != 6 & rowSums(!is.na(filtered2[,5:10])) > 0, ]
# Haploid = 6, but submetacentric + metacentric <= 4
hap6_low_meta <- filtered2[filtered2$Haploid.Number == 6 & 
                             rowSums(filtered2[, c("Submetacentric","Metacentric")], na.rm = TRUE) <= 4, ]
excluded_species <- rbind(non6, hap6_low_meta)
excluded_species$linke6  <- rowSums(excluded_species[, c("Submetacentric","Metacentric")], na.rm = TRUE)
excluded_species$linke8  <- rowSums(excluded_species[, c("Subtelocentric","Telocentric")], na.rm = TRUE)
excluded_species$linke10 <- rowSums(excluded_species[, c("Subacrocentric","Acrocentric")], na.rm = TRUE)
excluded_species <- excluded_species[, -c(5:10)]
colnames(excluded_species)[colnames(excluded_species) == "linke6"]  <- "metacentric"
colnames(excluded_species)[colnames(excluded_species) == "linke8"]  <- "telocentric"
colnames(excluded_species)[colnames(excluded_species) == "linke10"] <- "acrocentric"
head(excluded_species)
dim(excluded_species)
unique_families <- unique(excluded_species$Family)
unique_families
# "Chloropidae" "Glossinidae" "Lauxaniidae" "Muscidae"    "Phoridae"    "Sciomyzidae" "Tephritidae" "Ulidiidae"  
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# plot the heatmap

# Combine the two datasets
combined_species <- rbind(
  excluded_species[, c("Haploid.Number", "metacentric")],
  filtered5[, c("Haploid.Number", "metacentric")]
)
# remove invalid samples
combined_species <- combined_species[
  combined_species$metacentric <= combined_species$Haploid.Number, 
]

# Recreate the table
haplo_meta <- table(combined_species$Haploid.Number,
                    combined_species$metacentric)

haplo_meta_matrix <- as.matrix(haplo_meta)
m <- haplo_meta_matrix
m_log10 <- log10(m + 1)  # log2(number + 1) to handle zeros

# Create color palette
cols <- colorRampPalette(c("lightyellow", "steelblue", "darkred"))(110)

breaks <- seq(0, max(m_log10, na.rm=TRUE), length.out=111)

# Define legend labels at nice intervals of original numbers
legend_vals <- pretty(c(0, max(m)), n=5)  # original numbers
legend_at <- log10(legend_vals + 1)

# Plot
image.plot(1:ncol(m_log10), 1:nrow(m_log10),
           t(m_log10[nrow(m_log10):1, ]),
           col=cols,
           breaks=breaks,
           axes=FALSE,
           xlab="Metacentric chromosomes",
           ylab="Haploid number",
           legend.lab="Number of species",
           legend.cex = 1.3,
           cex.lab =2,
           legend.line=2.5,
           legend.shrink=1,
           axis.args=list(at=legend_at, labels=legend_vals, cex.axis=1.3))

axis(1, at=1:ncol(m_log10), labels=colnames(m_log10), las=2, cex.axis=1.5)
axis(2, at=1:nrow(m_log10), labels=rev(rownames(m_log10)), las=1, cex.axis=1.5)
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# save the plot
png("/Users/jg40/Desktop/m_log_plot.png",
    width = 6, height = 5, units = "in", res = 300)

par(mar = c(5, 5, 1, 1))  # bottom, left, top, right

# Create color palette
cols <- colorRampPalette(c("lightyellow", "steelblue", "darkred"))(110)

breaks <- seq(0, max(m_log10, na.rm=TRUE), length.out=111)

# Define legend labels at nice intervals of original numbers
legend_vals <- pretty(c(0, max(m)), n=5)  # original numbers
legend_at <- log10(legend_vals + 1)

# Plot
image.plot(1:ncol(m_log10), 1:nrow(m_log10),
           t(m_log10[nrow(m_log10):1, ]),
           col=cols,
           breaks=breaks,
           axes=FALSE,
           xlab="metacentric chromosomes",
           ylab="haploid number",
           legend.lab="Number of species",
           legend.cex = 1.3,
           cex.lab =2,
           legend.line=2.5,
           legend.shrink=1,
           axis.args=list(at=legend_at, labels=legend_vals, cex.axis=1.0))

axis(1, at=1:ncol(m_log10), labels=colnames(m_log10), las=2, cex.axis=1.5)
axis(2, at=1:nrow(m_log10), labels=rev(rownames(m_log10)), las=1, cex.axis=1.5)

dev.off()
#------------------------------------------------------------------------------#
