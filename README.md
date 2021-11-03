# ox-stress_histones

**Perturbed fatty-acid metabolism is linked to chromatin hyperacetylation, increased stress-response gene expression and resistance to oxidative stress**

Jarmila Princová, Clàudia Salat-Canela, Petr Daněk, Anna Marešová, Laura de Cubas, Jürg Bähler, José Ayté, Elena Hidalgo, and Martin Převorovský

<br>

**Figure 4**

Figure4.R - R script for generating Figure 4A and 4B

files - list of expression microarray data files

_expression microarray data files:_ <br>
These are fully processed data from https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6761/. <br>
martinp_2436_11_WT_0min_norm_genes.txt <br>
martinp_2436_12_WT_15min_norm_genes.txt <br>
martinp_2436_13_cbf11KO_0min_norm_genes.txt <br>
martinp_2436_14_cbf11KO_15min_norm_genes.txt <br>
martinp_2436_15_cbf12KO_0min_norm_genes.txt <br>
martinp_2436_16_cbf12KO_15min_norm_genes.txt <br>
martinp_2436_22_cbf12KO_60min_norm_genes.txt <br>
martinp_2436_23_cbf11KO_60min_norm_genes.txt <br>
martinp_2436_24_WT_60min_norm_genes.txt

<br>
  
**Figure 5 and Supplementary Figure 4**

Histone_ChIP-seq_analysis_CPM.Rmd - Rmarkdown document with all code needed for getting coverage tracks from raw ChIP-seq data (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11081/) <br>
Histone_ChIP-seq_analysis_CPM.html

profiles_CPM-normalized.Rmd - Rmarkdown document with code needed for generating Figure 5A and Supplementary Figure 4 <br>
profiles_CPM-normalized.html

_gene lists:_ <br>
CESRup_All.txt - all CESR genes upregulated during stress (taken from https://doi.org/10.1091/mbc.e02-08-0499) <br>
cbf11 2xUP (3of4) in YES log.txt - genes upregulated in _cbf11KO_ grown in YES (taken from https://doi.org/10.1371/journal.pone.0137820)

S. pombe _genome annotation:_ <br>
Schizosaccharomyces_pombe_all_chromosomes.gff3.gz (obtained from https://www.pombase.org/data/genome_sequence_and_features/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz)
