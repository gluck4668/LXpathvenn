
library(openxlsx)
gene_KEGG_pathways <- read.xlsx("gene KEGG pathways.xlsx")
pro_KEGG_pathways <- read.xlsx("pro KEGG pathways.xlsx")
metabolite_kegg_pathways <- read.xlsx("metabolite_kegg_pathways.xlsx")

usethis::use_data(gene_KEGG_pathways,overwrite = T)
usethis::use_data(pro_KEGG_pathways,overwrite = T)
usethis::use_data(metabolite_kegg_pathways,overwrite = T)

rm(list=ls())

data(gene_KEGG_pathways)
data(pro_KEGG_pathways)
data(metabolite_kegg_pathways)
