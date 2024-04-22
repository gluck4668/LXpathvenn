
library(openxlsx)
gene_KEGG_pathways <- read.xlsx("gene KEGG pathways.xlsx")
pro_KEGG_pathways <- read.xlsx("pro KEGG pathways.xlsx")
meta_pathways_impala <- read.csv("meta_pathways_impala.csv")
meta_pathways_metaboAnalyst <- read.csv("meta_pathways_metaboAnalyst.csv")

usethis::use_data(gene_KEGG_pathways,overwrite = T)
usethis::use_data(pro_KEGG_pathways,overwrite = T)
usethis::use_data(meta_pathways_impala,overwrite = T)
usethis::use_data(meta_pathways_metaboAnalyst,overwrite = T)

rm(list=ls())

data(gene_KEGG_pathways)
data(pro_KEGG_pathways)
data(meta_pathways_impala) # obtained from http://impala.molgen.mpg.de/
data(meta_pathways_metaboAnalyst) # obtained from https://www.metaboanalyst.ca/
