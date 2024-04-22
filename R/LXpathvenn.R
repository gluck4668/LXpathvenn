
LXpathvenn <- function(gene_path_data,protein_path_data,meta_path_data){

packs <- c("openxlsx","dplyr","ggplot2","stringr","readr")
lib_fun <- function(i){library(i,character.only = T)}
sapply(packs,lib_fun)

#-----gene pathways data-----------
gene_type <- str_extract(gene_path_data,"(?<=[.]).*")
if(tolower(gene_type) =="txt")
   gene_path <- read_table(gene_path_data) else
   gene_path <- eval(str2expression(paste0("read.",gene_type,"(gene_path_data)")))
names(gene_path)
gene_pathways <- gene_path[,c("Description","pvalue","GeneRatio","geneSymbol")]
names(gene_pathways) <- c("pathways","gene_p","GeneRatio","geneSymbol")
gene_pathways$pathways <- trimws(gene_pathways$pathways)

#-------protein pathways data-------
protein_type <- str_extract(protein_path_data,"(?<=[.]).*")
if(tolower(protein_type) =="txt")
  protein_path <- read_table(protein_path_data) else
    protein_path <- eval(str2expression(paste0("read.",protein_type,"(protein_path_data)")))
names(protein_path)
protein_pathways <- protein_path[,c("Pathways","pvalue","Protein/Ratio","ProteinID")]
names(protein_pathways) <- c("pathways","protein_p","ProteinRatio","ProteinID")
protein_pathways$pathways <- trimws(protein_pathways$pathways)

num_fun <- function(i){eval(str2expression(i))} # 转变为小数
protein_pathways$ProteinRatio <- sapply(protein_pathways$ProteinRatio,num_fun)

#------metabolite pathways data------
meta_type <- str_extract(meta_path_data,"(?<=[.]).*")
if(tolower(meta_type) =="txt")
  meta_path <- read_table(meta_path_data) else
    meta_path <- eval(str2expression(paste0("read.",meta_type,"(meta_path_data)")))

if(grepl("overlapping",names(meta_path),ignore.case = T) %>% any())
  {
   meta_pathways <- filter(meta_path,pathway_source=="KEGG")
   meta_pathways$num_all_pathway_metabolites <- str_extract(meta_pathways$num_all_pathway_metabolites,".*(?=[(])") %>% as.numeric()
   meta_pathways$MetaRatio <- meta_pathways$num_overlapping_metabolites/meta_pathways$num_all_pathway_metabolites
   meta_pathways <- meta_pathways[,c("pathway_name","P_metabolites","MetaRatio","overlapping_metabolites")]
   names(meta_pathways) <-  c("pathways","metabolite_p","MetaRatio","Metabolites")
   meta_pathways$pathways <- str_extract(meta_pathways$pathways,".*(?= -)")
  }else
   {
   meta_pathways <- meta_path[,c(1,5,4,2)]
   meta_pathways$MetaRatio <- meta_pathways$Hits/meta_pathways$Total
   meta_pathways <- meta_pathways[,-c(3,4)]
   names(meta_pathways) <-  c("pathways","metabolite_p","MetaRatio")
    }
 meta_pathways$pathways <- trimws(meta_pathways$pathways)

#-----venn analysis----------
venn_path_all <- inner_join(gene_pathways,protein_pathways,"pathways") %>%
                  inner_join(.,meta_pathways,"pathways")


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

#---------ggplot point----------------
mytheme<-theme_bw()+
  theme(text=element_text(family = "sans",colour ="black",face="bold",size =12),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))+
  theme(panel.grid =element_line(colour="#dcdcdc",linewidth=0.2,linetype = "dashed"))


xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=10,angle =0,hjust=1))+
  theme(axis.text.y = element_text(face="bold",color="black",size=10))+
  theme(legend.text=element_text(face="bold",color="black",size=10))


venn_paths_p005 <- data.frame(pathways=rep(venn_path_p005$pathways,3),
                         pvalue=c(venn_path_p005$gene_p,venn_path_p005$protein_p,venn_path_p005$metabolite_p),
                         Ratio=c(venn_path_p005$GeneRatio,venn_path_p005$ProteinRatio,venn_path_p005$MetaRatio),
                         Type=c(rep("gene",nrow(venn_path_p005)),rep("protein",nrow(venn_path_p005)),rep("metabolite",nrow(venn_path_p005)))
                         )


f5 <- ggplot(venn_paths_p005)+
  geom_point(aes(x=Ratio,
                 y=pathways,
                 shape=Type,
                 color=-log2(pvalue),
                 size=Ratio))+
  scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
  labs(x = 'Ratio', y = 'Pathways',title=c('Gene-Protein-Metabolite Joint Pathways'),
       col="-log2(pvalue)",shape="Type",size="Ratio")+ # 修改图例名称
  mytheme+xytheme

f5

f5_name <-c("Gene-Protein-metabolite Joint pathways 03.png")
f5_name <-paste0(dir.file,"/", f5_name)
ggsave(f5_name,f5,width=1450, height =1200, dpi=150,units = "px")

#------------------
line05 <- geom_vline(xintercept = c(log05),
                     linewidth = 0.5,
                     color = "black",
                     lty = "dashed")

txt05 <- geom_text(x=log05+4,y=2.5,label = c("p<0.05"),
                   size=5,color="blue",fontface="italic")

f6 <- ggplot(venn_paths_p005)+
  geom_point(aes(x=-log2(pvalue),
                 y=pathways,
                 shape=Type,
                 color=-log2(pvalue),
                 size=Ratio))+
  scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
  labs(x = '-log2(pvalue)', y = 'Pathways',title=c('Gene-Protein-Metabolite Joint Pathways'),
       col="-log2(pvalue)",shape="Type",size="Ratio")+ # 修改图例名称
  mytheme+xytheme+
  line05+txt05

f6

f6_name <-c("Gene-Protein-metabolite Joint pathways 04.png")
f6_name <-paste0(dir.file,"/", f6_name)
ggsave(f6_name,f6,width=1450, height =1200, dpi=150,units = "px")


#--------------
f4

}






