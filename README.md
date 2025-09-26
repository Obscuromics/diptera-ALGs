# diptera-ALGs

## Basic Genome Stats

To get basic summaries about the genomes run following script

```bash
Rscript scripts/20250702_genome_summaries.R 
```

## Writing

- [Julia's thesis](https://docs.google.com/document/d/1bfhqvAD_j3B5yVDeCrSt0J82xKEFOqUfy64Hq547xIE)
- [Drafting space](https://docs.google.com/document/d/1fc5yGJq-1laB7v4YSdKUrmTqiotB-RvVQa6od_DRAJU/edit?tab=t.0)
- [Clean manuscript](https://docs.google.com/document/d/18QSRY8k8Z5xc9gZTiWP54esINj1GxGQF7-9OgIyaJoo/edit?tab=t.0)
- [Supps](https://docs.google.com/document/d/1qPkfJt6ZOQZzHlvTzJtfGthx2J_9wtrOV4KQKjn2PjE/edit?tab=t.0#heading=h.32dn0pko8wct)

## Figures

 - [drafting slideck](https://docs.google.com/presentation/d/1M6Vxeka0wPZDs_RF7VTUbJoILGp6jp1BJcpDHgwQSoI/edit?usp=sharing)


#### Figure 3

Sex chromosome panel

```
scripts/figure_3/20250228_sex_chromosomes_in_all_species.R
```

## Others

- [The task list](https://github.com/orgs/Obscuromics/projects/1)


### plotting dotplots


```bash
Rscript scripts/plot_dotplot.R -s1 Bibio_marci -s2 Ptychoptera_contaminata -o figures/dotplot_test

Rscript scripts/plot_dotplot.R -s1 Anopheles_ziemanni -s2 Anopheles_coustani -o figures/dotplot_Ano_ziemanni_coustani
Rscript scripts/plot_dotplot.R -s1 Anopheles_moucheti -s2 Anopheles_marshallii -o figures/dotplot_Ano_moucheti_marshallii
Rscript scripts/plot_dotplot.R -s1 Anopheles_pretoriensis -s2 Anopheles_maculipalpis -o figures/dotplot_Ano_pretoriensis_maculipalpis
Rscript scripts/plot_dotplot.R -s1 Anopheles_squamosus -s2 Anopheles_pharoensis -o figures/dotplot_Ano_squamosus_pharoensis
Rscript scripts/plot_dotplot.R -s1 Anopheles_longipalpis -s2 Anopheles_stephensi -o figures/dotplot_Ano_longipalpis_stephensi
Rscript scripts/plot_dotplot.R -s1 Culex_modestus -s2 Aedes_aegypti -o figures/dotplot_Culex_modestus_Aedes_aegypti
 

 
Rscript scripts/plot_dotplot.R -s1 Drosophila_sechellia -s2 Drosophila_melanogaster -o figures/dotplot_Dros_1_sechellia_melanogaster
Rscript scripts/plot_dotplot.R -s1 Drosophila_triauraria -s2 Drosophila_melanogaster -o figures/dotplot_Dros_triauraria_melanogaster
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




```