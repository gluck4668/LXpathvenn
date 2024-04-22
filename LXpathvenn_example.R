
if(!requireNamespace("devtools"))
  install.packages("devtools")
library(devtools)

install_github("gluck4668/LXpathvenn")
library(LXpathvenn)

#---------------------------------
data(gene_KEGG_pathways)
data(pro_KEGG_pathways)
data(meta_pathways_impala) # obtained from http://impala.molgen.mpg.de/
data(meta_pathways_metaboAnalyst) # obtained from https://www.metaboanalyst.ca/
#------------------------------
rm(list=ls())

devtools::load_all()

gene_path_data = "gene KEGG pathways.xlsx"
protein_path_data="pro KEGG pathways.xlsx"
meta_path_data="meta_pathways_impala.csv"
# meta_path_data="meta_pathways_metaboAnalyst.csv"

LXpathvenn(gene_path_data,protein_path_data,meta_path_data)
