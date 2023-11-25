\name{LXpathvenn}
\alias{LXpathvenn}
\title{LXpathvenn}
\usage{
LXpathvenn(gene_path_data,protein_path_data,meta_pathe_data)
}
\description{
LXpathvenn can be used to analyze the jont pathways.
}
\examples{
if(!requireNamespace("devtools"))
  install.packages("devtools")
library(devtools)

install_github("gluck4668/LXpathvenn")
library(LXpathvenn)

#---------------------------------
data(gene_KEGG_pathways)
data(pro_KEGG_pathways)
data(metabolite_kegg_pathways)
#------------------------------
rm(list=ls())

#devtools::load_all()

gene_path_data = "gene KEGG pathways.xlsx"
protein_path_data="pro KEGG pathways.xlsx"
meta_pathe_data="metabolite_kegg_pathways.xlsx"

LXpathvenn(gene_path_data,protein_path_data,meta_pathe_data)
}
