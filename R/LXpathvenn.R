
LXpathvenn <- function(gene_path_data,protein_path_data,meta_pathe_data){

library(openxlsx)
library(dplyr)
library(ggplot2)

gene_path <- read.xlsx(gene_path_data)[,c(2,5,9)]
names(gene_path) <- c("pathways","gene_p","genes")
gene_path$pathways <- trimws(gene_path$pathways)

protein_path <- read.xlsx(protein_path_data)[,c(2,5,8)]
names(protein_path) <-c("pathways","protein_p","proteins")
protein_path$pathways <- trimws(protein_path$pathways)

meta_path <- read.xlsx(meta_pathe_data)[,c(1,6,4)]
names(meta_path) <- c("pathways","metabolite_p","metabolites")
meta_path$pathways <- trimws(meta_path$pathways)


venn_path_all <- inner_join(gene_path,protein_path,"pathways") %>%
                  inner_join(.,meta_path,"pathways")


venn_path_p005 <- filter(venn_path_all,metabolite_p<0.05)

if(!dir.exists("analysis results"))
  dir.create("analysis results")
  dir.file="analysis results"

write.xlsx(venn_path_all,"analysis results/gene_protein_metabolites_kegg_pathways_all.xlsx")

write.xlsx(venn_path_p005,"analysis results/gene_protein_metabolites_kegg_pathways(metabolite_p 0.05).xlsx")

#---------Jiont analysis of the gene and metabolite enriched pathways---------#

venn_path <- venn_path_p005

gene_pro_meta_path <- data.frame(pathways=rep(venn_path$pathways,3),
                                 log2_p=c(-log2(venn_path$gene_p),-log2(venn_path$protein_p),-log2(venn_path$metabolite_p)),
                                 type=c(rep("genes",nrow(venn_path)),rep("proteins",nrow(venn_path)),rep("metabolites",nrow(venn_path)))
                                  )

y_p <- max(gene_pro_meta_path$log2_p)

height_y <- y_p*1.2

nrow_path <- nrow(gene_pro_meta_path)/2

joint_title_size <- case_when(nrow_path>=30 ~12,
                              nrow_path>=20 ~12,
                              TRUE ~14)

joint_x_size <- case_when(nrow_path>=20 ~9,
                          nrow_path>=10 ~10,
                          TRUE ~12)

joint_y_size <- case_when(nrow_path>=20 ~12,
                          nrow_path>=10 ~12,
                          TRUE ~12)

joint_legend_size <- case_when(nrow_path>=30 ~12,
                               nrow_path>=20 ~12,
                               TRUE ~12)


bar_width <- case_when(nrow_path>=30 ~0.9,
                       nrow_path>=20 ~0.8,
                       TRUE ~0.7)

f1 <- ggplot(gene_pro_meta_path, aes(x = pathways, y = log2_p,fill=type))+
  geom_bar(position = "dodge",stat = "identity",width = bar_width)+
  scale_fill_manual(values=c("#008b8b","#f08080","#87ceeb"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(colour = "black", face="bold",size=12)) +
  labs(x="",y = "-log2(Pvalue)",title = 'Gene-protein-metabolite Joint Pathways')+
  scale_y_continuous(expand = c(0, 0),limits = c(0, height_y))

f1

log05 <- -log2(0.05)
log01 <- -log2(0.01)

line1 <- geom_hline(yintercept = c(log05),
                    linewidth = 0.6,
                    color = "blue",
                    lty = "dashed")
line2 <- geom_hline(yintercept = c(log01),
                    linewidth = 0.6,
                    color = "red",
                    lty = "dashed")

y1 <- geom_text(x=nrow(gene_pro_meta_path)/3,y=log05+0.8,label = c("p<0.05"),
                size=4,color="blue",fontface="italic")
y2 <- geom_text(x=nrow(gene_pro_meta_path)/3,y=log01+0.8,label = c("p<0.01"),
                size=4,color="blue",fontface="italic")

f2 <- f1+line1+line2+y1+y2

mytheme<-theme_bw()+
  theme(text=element_text(family = "sans",colour ="black",face="bold",size =joint_title_size),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"),
        plot.margin = unit(c(t=0.5, r=0.5, b=0.5, l=2), "cm")
  )+
  theme(plot.title = element_text(hjust = 0.5))


xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=joint_x_size,angle =45,hjust=1))+
  theme(axis.text.y = element_text(face="bold",color="black",size=joint_y_size))

legend_theme <- theme(
  legend.title = element_blank(),
  legend.text = element_text(size = joint_legend_size, face = "bold"),
  legend.direction = "vertical",
  #legend.position = c(0.5,0.9),
  legend.background = element_blank()
)

f3 <- f2+mytheme+xytheme


f3_name <- paste("Gene_protein_metabolite_Joint_pathways 01.png")
f3_name <-paste0(dir.file,"/", f3_name)

ggsave(f3_name,f3,width=1200, height =1000, dpi=150,units = "px")


f4 <- f3+
  labs(fill="")+
  theme(legend.direction = "horizontal",
        legend.position = c(0.5,0.92),
        legend.text = element_text(size=14,face = "bold") )

f4_name <- paste("Gene_protein_metabolite_Joint_pathways 02.png")
f4_name <-paste0(dir.file,"/", f4_name)

ggsave(f4_name,f4,width=1200, height =1000, dpi=150,units = "px")

print("--------------------------------------------------------------")

print_text <- paste("The results can be found in the folder of",dir.file)

print(print_text)

f4

}






