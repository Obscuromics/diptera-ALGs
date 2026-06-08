# 340 dipteran genomes reveal the origin of Muller elements and sex chromosomes in Diptera

This is a repository accompaning this preprint: https://www.biorxiv.org/content/10.64898/2026.06.01.729285v1

## Organisation

This repository was used for the work, therefore the history contains lots of code specific to our cluster environment and the scripts still need some reorganisation and polishing. This repository contains only the code and some small files and tables

### data directly in the repository 

The tree based on 5067 BUSCO genes found in 340 dipteran genomes and a mecopteran _Panorpa germanica_ as an outgroup. Sequences were alinged using MAFFT v7.525, trimAl v1.4.rev15 using `-gt 0.8 -st 0.001 -resoverlap 0.75 -seqoverlap 80` parameters. 4298 alignment passing the trimming step were used to infer the tree using iqtree with autoselection of the substitution model (Q.insect+I+G4 substitution model with gamma rate variation) and 1000 ultrafast bootstrap replicates. The final tree is available in this repository:

```
data/syngraph/iqtree/diptera.supermatrix.phy.treefile
```

The ancestral linkage groups were infered using syngraph (see Methods of the manuscript). Here are the final assignments of BUSCOs to ALGs for common ancestor of all Diptera: `tables/ALGs_syngraph_diptera.tsv` (karyotype α), Brachycera: `tables/ALGs_syngraph_brachycera.tsv` (karyotype β). We include also assignments for the common ancestor of schizophora and relatives (karyotype γ), however, we would advise using Brachyceran ALGs instead, as typical schizophoran chormosomes are metacentric chromosomes with two chromosome arms corresponding to two different brachyceran ALGs. This ALGs are found here: `ALGs_syngraph_schizophora_syrphidae.tsv`. Karyotype γ is moreless corresponding to Muller Elements (the biggest difference is in 'the dot').

### data to fetch from a cloud

TODO
This section will be done after finalisation of all the files during the review process.

## Rerunning scripts

TODO (most of the script will have missing data)

## Basic Genome Stats


## Regenerating figures

 - [drafting slideck](https://docs.google.com/presentation/d/1M6Vxeka0wPZDs_RF7VTUbJoILGp6jp1BJcpDHgwQSoI/edit?usp=sharing)

### Plots

Genric

```bash
mkdir -p figures/paints
# dotplot
Rscript scripts/plot_dotplot.R -s1 Bibio_marci -s2 Ptychoptera_contaminata -o figures/dotplot_test
# paintings
Rscript scripts/ALG_painter_ME_building_blocks.R -a tables/ALGs_syngraph_diptera.tsv -l tables/species_lists/selected_mosquitos.tsv -o figures/paints/mosquito_pick
Rscript scripts/ALG_painter_ME_building_blocks.R -a tables/ALGs_syngraph_diptera.tsv -f Tabanidae -o figures/paints/Tabanidae
# table of chromosome vs ancestry
Rscript scripts/table_of_chromosomes_vs_ancestry.R -s Machimus_atricapillus -a tables/ALGs_syngraph_diptera.tsv
```

#### Figure 1

#### Figure 2

```bash
Rscript scripts/figure_2/20250627_figure_2_tree_with_families.R
Rscript scripts/misc/plot_tree_with_ALGs_at_nodes.R -t data/syngraph/diptera.no_plecia.mindist.m165.newick.txt -a tables/ALGs_syngraph_diptera.tsv -l family -o figures/diptera_syngraph_tree_of_changes -r data/syngraph/diptera.no_plecia.mindist.m165.rearrangements.tsv
```

TODO: Julia needs to add centromere plotting scripts and learn how to ust git



#### Figure 3

Sex chromosome panel

```
scripts/figure_3/20250228_sex_chromosomes_in_all_species.R
```

#### Supps

Suppmenetary figure showing % of ALGs together on chromosomes of contemporary species.

```bash
Rscript scripts/20251106_stability.R -a tables/ALGs_syngraph_schizophora_syrphidae.tsv -o figures/ALG_stability_histograms_ds
Rscript scripts/20251106_stability.R -a tables/ALGs_syngraph_brachycera.tsv -o figures/ALG_stability_histograms_db
Rscript scripts/20251106_stability.R -a tables/ALGs_syngraph_diptera.tsv -o figures/ALG_stability_histograms_d
```

#### Paintings

```bash
Rscript scripts/ALG_painter_ME_building_blocks.R -a tables/ALGs_syngraph_diptera.tsv -l tables/species_lists/selected_mosquitos.tsv -o figures/paints/mosquito_pick # selected culicomorphan genomes
Rscript scripts/ALG_painter_ME_building_blocks.R -a tables/ALGs_syngraph_diptera.tsv -l tables/species_lists/excluded_species -o figures/paints/excluded
```

#### plotting dotplots

```bash
Rscript scripts/plot_dotplot.R -s1 Bibio_marci -s2 Ptychoptera_contaminata -o figures/dotplot_test

Rscript scripts/plot_dotplot.R -s1 Anopheles_ziemanni -s2 Anopheles_coustani -o figures/dotplot_Ano_ziemanni_coustani
Rscript scripts/plot_dotplot.R -s1 Anopheles_moucheti -s2 Anopheles_marshallii -o figures/dotplot_Ano_moucheti_marshallii
Rscript scripts/plot_dotplot.R -s1 Anopheles_pretoriensis -s2 Anopheles_maculipalpis -o figures/dotplot_Ano_pretoriensis_maculipalpis
Rscript scripts/plot_dotplot.R -s1 Anopheles_squamosus -s2 Anopheles_pharoensis -o figures/dotplot_Ano_squamosus_pharoensis
Rscript scripts/plot_dotplot.R -s1 Anopheles_longipalpis -s2 Anopheles_stephensi -o figures/dotplot_Ano_longipalpis_stephensi
Rscript scripts/plot_dotplot.R -s1 Culex_modestus -s2 Aedes_aegypti -o figures/dotplot_Culex_modestus_Aedes_aegypti
 

 
Rscript scripts/plot_dotplot.R -s1 Drosophila_sechellia -s2 Drosophila_melanogaster -o figures/dotplot_Dros_1_sechellia_melanogaster
Rscript scripts/misc/plot_dotplot.R -s1 Drosophila_triauraria -s2 Drosophila_melanogaster -o figures/dotplot_Dros_triauraria_melanogaster
Rscript scripts/plot_dotplot.R -s1 Drosophila_yakuba -s2 Drosophila_melanogaster -o figures/dotplot_Dros_yakuba_melanogaster
Rscript scripts/plot_dotplot.R -s1 Drosophila_ananassae -s2 Drosophila_melanogaster -o figures/dotplot_Dros_ananassae_melanogaster
Rscript scripts/plot_dotplot.R -s1 Drosophila_pseudoobscura -s2 Drosophila_melanogaster -o figures/dotplot_Dros_pseudoobscura_melanogaster
Rscript scripts/plot_dotplot.R -s1 Drosophila_funebris -s2 Drosophila_melanogaster -o figures/dotplot_Dros_funebris_melanogaster
 
# TODO
Rscript scripts/plot_dotplot.R -s1 Tipula_fascipennis -s2 Tipula_helvola -o figures/dotplot_Tipula_1_fascipennis_Tipula_helvola
Rscript scripts/plot_dotplot.R -s1 Tipula_fascipennis -s2 Tipula_lateralis -o figures/dotplot_Tipula_2_fascipennis_Tipula_lateralis
Rscript scripts/plot_dotplot.R -s1 Tipula_fascipennis -s2 Tipula_unca -o figures/dotplot_Tipula_3_fascipennis_Tipula_unca
Rscript scripts/plot_dotplot.R -s1 Tipula_fascipennis -s2 Tipula_scripta -o figures/dotplot_Tipula_4_fascipennis_Tipula_scripta
Rscript scripts/plot_dotplot.R -s1 Tipula_fascipennis -s2 Nephrotoma_flavescens -o figures/dotplot_Tipula_5_fascipennis_Nephrotoma_flavescens


Rscript scripts/plot_dotplot.R -s1 Tipula_fascipennis -s2 Tipula_helvola -o figures/dotplot_Tipula_1_fascipennis_Tipula_helvola
Rscript scripts/plot_dotplot.R -s1 Tipula_fascipennis -s2 Tipula_lateralis -o figures/dotplot_Tipula_2_fascipennis_Tipula_lateralis
Rscript scripts/plot_dotplot.R -s1 Tipula_fascipennis -s2 Tipula_unca -o figures/dotplot_Tipula_3_fascipennis_Tipula_unca
Rscript scripts/plot_dotplot.R -s1 Tipula_fascipennis -s2 Tipula_scripta -o figures/dotplot_Tipula_4_fascipennis_Tipula_scripta
Rscript scripts/plot_dotplot.R -s1 Tipula_fascipennis -s2 Nephrotoma_flavescens -o figures/dotplot_Tipula_5_fascipennis_Nephrotoma_flavescens

Rscript scripts/plot_dotplot.R -s1 Chironomus_riparius -s2 Chironomus_tentans -o figures/_dotplots/dotplot_Chironomus_riparius_1_Chironomus_tentans
Rscript scripts/plot_dotplot.R -s1 Chironomus_riparius -s2 Chironomus_striatipennis -o figures/_dotplots/dotplot_Chironomus_riparius_2_Chironomus_striatipennis
Rscript scripts/plot_dotplot.R -s1 Chironomus_riparius -s2 Polypedilum_vanderplanki -o figures/_dotplots/dotplot_Chironomus_riparius_3_Polypedilum_vanderplanki
Rscript scripts/plot_dotplot.R -s1 Chironomus_riparius -s2 Smittia_aterrima -o figures/_dotplots/dotplot_Chironomus_riparius_4_Smittia_aterrima
Rscript scripts/plot_dotplot.R -s1 Chironomus_riparius -s2 Propsilocerus_akamusi -o figures/_dotplots/dotplot_Chironomus_riparius_5_Smittia_aterrima

Rscript scripts/plot_dotplot.R -s1 Chironomus_riparius -s2 Propsilocerus_akamusi -o figures/_dotplots/dotplot_Chironomus_riparius_5_Smittia_aterrima

Rscript scripts/misc/plot_dotplot.R -s1 Merodon_equestris -s2 Blera_fallax -o figures/_dotplots/dotplot_Merodon_equestris_Blera_fallax
```
