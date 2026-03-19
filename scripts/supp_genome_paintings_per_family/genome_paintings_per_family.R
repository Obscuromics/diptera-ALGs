# =========================
# Load required libraries
# =========================
library(ape)

# =========================
# File paths
# =========================
setwd("scripts/supp_genome_paintings_per_family")
tree_file <- "diptera.supermatrix.phy.treefile"
species_file <- "diptera_species_340.tsv"

# =========================
# Load data
# =========================
tree <- read.tree(tree_file)

fam_data <- read.delim(
  species_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

# Quick inspection
head(fam_data)

# =========================
# Prepare data structures
# =========================

# Species present in the tree
tree_species <- tree$tip.label

# Group species by family
family_species <- split(fam_data$species, fam_data$family)

# =========================
# Helper functions
# =========================

# Subset tree by tips
get_family_tree <- function(tree, species_vector) {
  keep.tip(tree, species_vector)
}

# Reorder dataframe to match tree tip order
order_df_by_tree <- function(df, tree_species) {
  df_ordered <- df[match(tree_species, df$species), ]
  df_ordered <- df_ordered[, c("species", "family", "accession")]
  df_ordered
}

# Reverse dataframe rows
reverse_df <- function(df) {
  df_rev <- df[nrow(df):1, ]
  rownames(df_rev) <- NULL
  df_rev
}

# Convert vector to dataframe row
vec_to_df <- function(v) {
  data.frame(
    species = v[1],
    family = v[2],
    accession = v[3],
    stringsAsFactors = FALSE
  )
}

# =========================
# Build trees per family
# =========================
family_trees <- lapply(family_species, function(sp) {
  get_family_tree(tree, sp)
})

# =========================
# Build ordered + reversed dataframes per family
# =========================
all_family_dfs <- lapply(names(family_trees), function(fam) {
  
  fam_tree <- family_trees[[fam]]
  
  # Match dataframe to tree order
  fam_df_ordered <- order_df_by_tree(fam_data, fam_tree$tip.label)
  
  # Reverse order (your original logic)
  fam_df_reversed <- reverse_df(fam_df_ordered)
  
  fam_df_reversed
})

names(all_family_dfs) <- names(family_trees)

# =========================
# Quick checks / exploration
# =========================
head(all_family_dfs)
head(all_family_dfs[["Syrphidae"]])

# Plot example family tree
plot(family_trees[["Syrphidae"]], cex = 0.7)

# =========================
# Add additional species
# =========================
new_species <- list(
  c("Culicoides_sonorensis", "Ceratopogonidae", "GCA_047716325.1"),
  c("Forcipomyia_palustris", "Ceratopogonidae", "GCA_976917515.1"),
  c("Atylotus_latistriatus", "Tabanidae", "GCA_977012865.1")
)

# Add each new species to the corresponding family
for (v in new_species) {
  
  new_row <- vec_to_df(v)
  fam <- new_row$family
  
  if (fam %in% names(all_family_dfs)) {
    # Append to existing family
    all_family_dfs[[fam]] <- rbind(all_family_dfs[[fam]], new_row)
  } else {
    # Create new family entry
    all_family_dfs[[fam]] <- new_row
  }
}

# =========================
# Final checks
# =========================
all_family_dfs[["Ceratopogonidae"]]
all_family_dfs[["Tabanidae"]]
