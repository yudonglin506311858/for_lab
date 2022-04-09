Sys.setenv(LANGUAGE = "en") #显示英文报错信息
gc()
memory.limit(9999999999)
set.seed(123)
rm(list = ls())  
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
library(pheatmap)
library("biomaRt")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("DESeq2")
library(edgeR)
library(limma)
library("Rgraphviz")
setwd("D:/anxiiuli_data/count_炎红")


dir.create("./new")
setwd("./new")

dir.create("./1为阈值")
setwd("./1为阈值")
dir.create("./2为阈值")
setwd("./2为阈值")

dir.create("./3为阈值")
setwd("./3为阈值")




#读入数据并合并

non_immune_ery13<-read.table("13-非炎红_FRAS220042674-1r.count",header=T)
non_immune_ery14<-read.table("14-非炎红_FRAS220042675-1r.count",header=T)
non_immune_ery15<-read.table("15-非炎红_FRAS220042676-1r.count",header=T)
immune_ery13<-read.table("13-炎红_FRAS220042667-1r.count",header=T)
immune_ery14<-read.table("14-炎红_FRAS220042668-1r.count",header=T)
immune_ery15<-read.table("15-炎红_FRAS220042669-1r.count",header=T)


#合并数据
data<-merge(non_immune_ery13,non_immune_ery14,by=c("Geneid","Chr","Start","End","Strand","Length"))
data<-merge(data,non_immune_ery15,by=c("Geneid","Chr","Start","End","Strand","Length"))
data<-merge(data,immune_ery13,by=c("Geneid","Chr","Start","End","Strand","Length"))
data<-merge(data,immune_ery14,by=c("Geneid","Chr","Start","End","Strand","Length"))
data<-merge(data,immune_ery15,by=c("Geneid","Chr","Start","End","Strand","Length"))


WT_count<-data

head(WT_count)
dim(WT_count)
colnames(WT_count)


#删除ensmusg的版本号

#WT_count$Geneid<-gsub("\\.*","",WT_count$Geneid)
#WT_count$Geneid <- gsub("\\.[0-9]*$", "", WT_count$Geneid)
rownames(WT_count)<-WT_count$Geneid
#https://www.nhooo.com/note/qa02b6.html
#https://www.biostars.org/p/178726/
WT_count_all<-WT_count[,c(1,6:12)]
head(WT_count_all)
colnames(WT_count_all)<-c("Geneid","Length","non_immune_ery13","non_immune_ery14","non_immune_ery15","immune_ery13","immune_ery14","immune_ery15"
                          )


cts<-WT_count_all
head(cts)
dim(cts)


cts<-cts[which(rowSums(cts[,c(3:8)]) > 0),]#去掉全为零的行 情况

write.csv(cts,"ensembl_gene_expression_count.csv")



#基因ID转换
library('biomaRt')
library("curl")
library(ensembldb)
library(dplyr)
library(AnnotationHub)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
saveRDS(mart,"mart_human.rds")
mart<-readRDS("mart_human.rds")

gene<-read.csv("ensembl_gene_expression_count_anxiuli.csv")
gene<-as.matrix(gene$X)
head(gene)
colnames(gene)[1]<-"ensembl_gene_id"
#listAttributes(mart)
id_con<-getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),filters = 'ensembl_gene_id', values = gene, mart = mart)
head(id_con)
write.csv(id_con,"human_gene_ensembl_transition.csv")





library(stringr)
#cts$Geneid<-str_sub(cts$Geneid,1,str_locate(cts$Geneid,"\\.")[1]-1)
id_con<-read.csv("D:/anxiiuli_data/count_炎红/mouse_gene_ensembl_transition.csv")
id_con<-id_con[,c(2,3)]
colnames(cts)[1]<-"ensembl_gene_id"
head(cts)
head(id_con)
cts<-cts[which(rowSums(cts[c(3:8)]) > 0),]
cts<-merge(id_con,cts,by=c("ensembl_gene_id"))
dim(cts)
write.csv(cts,file = "ternimalE_count_genesymbol.csv",row.names = F)
cts<-cts[,-1]
cts$external_gene_name<-make.names(cts$external_gene_name, unique = TRUE)
rownames(cts)<-cts$external_gene_name
write.table(cts,"data.txt", sep = "\t")



#读取count文件：
write.csv(cts,file="ternimalE_count.csv")

#读取tpm文件：
head(cts)
kb <- cts$Length / 1000
head(kb)
countdata <- cts[,3:8]
rpk <- countdata / kb
head(rpk)
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
head(tpm)
write.csv(tpm,file="ternimalE_tpm.csv")

#读取fpkm文件：
fpkm <- t(t(rpk)/colSums(countdata) * 10^6) 
head(fpkm)
write.csv(fpkm,file="ternimalE_fpkm.csv")



#PCA
# condition<-factor(c(
#   "non_immune_ery","non_immune_ery","non_immune_ery",
#   "immune_ery","immune_ery","immune_ery"),
#   levels = c("non_immune_ery","immune_ery"))
condition<-factor(c(
  "non_immune_ery","non_immune_ery","non_immune_ery",
  "immune_ery","immune_ery","immune_ery"),
  levels = c("immune_ery","non_immune_ery"))

head(cts)
tmp<-cts[,3:8]
head(tmp)
colData <- data.frame(row.names=colnames(cts[,3:8]), condition)
head(colData,10)
dds_all<- DESeqDataSetFromMatrix(countData = cts[,3:8],colData = colData,design= ~condition)
head(dds_all)
dds_all<- DESeq(dds_all)
vsd_all<-vst(dds_all,blind=FALSE)
head(vsd_all)
dist(t(assay(vsd_all)))
plotPCA(vsd_all,intgroup="condition")
vsd_all$condition<-factor(c("immune_ery13","immune_ery14","immune_ery15",
                            "non_immune_ery13","non_immune_ery14","non_immune_ery15"),levels = c("immune_ery13","immune_ery14","immune_ery15",
                                                                                                 "non_immune_ery13","non_immune_ery14","non_immune_ery15"))

#样本的聚类图
sampleDists <- dist(t(assay(vsd_all)))
library("RColorBrewer")
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_all$condition, vsd_all$type, sep="")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


#做一下差异表达分析
res_all <- results(dds_all)
head(res_all)
res_all <- res_all[order(res_all$padj),]
head(res_all)
diff_gene <- subset(res_all, padj < 0.05 & (log2FoldChange < -1|log2FoldChange > 1))
head(diff_gene)
dim(diff_gene)
head(res_all)
write.csv(diff_gene,file = "diff_gene_all.csv",row.names = T)
diff_gene_up <- subset(res_all, padj < 0.05 & (log2FoldChange > 1))
write.csv(diff_gene_up,file = "diff_gene_up.csv",row.names = T)
diff_gene_down <- subset(res_all, padj < 0.05 & (log2FoldChange < -1))
write.csv(diff_gene_down,file = "diff_gene_down.csv",row.names = T)
dim(diff_gene)
dim(diff_gene_up)
dim(diff_gene_down)
head(diff_gene_up)
head(diff_gene_down)



resdata <- merge(as.data.frame(res_all), as.data.frame(counts(dds_all, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file = "DEG_all.csv",row.names = F)


head(resdata)
resdata = res_all[order(resdata$pvalue),]
summary(resdata)
table(resdata$padj<0.05)#number of true 小于0.05 的基因个数

# grammar
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# analysis
library(DESeq2)

# graphics
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)
library(latex2exp)
diffData <- fread("DEG_all.csv")

colnames(diffData)[1] <- "gene"

diffData[is.na(padj), padj := 1][]
diffData[, p := -log10(padj)][]


diffData[, type := "ns"][]
diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]

labelGene1 <- diffData[order(p, decreasing = T)][type == "up"][1:10]
labelGene2 <- diffData[order(p, decreasing = T)][type == "down"][1:10]
labelGene <- rbind(labelGene1,labelGene2)

#pal_nejm()(8) %>% show_col()
typeColor <- structure(
  c(pal_nejm()(2), "gray80"),
  names = c("up", "down", "ns")
)
options(ggrepel.max.overlaps = Inf)
ggplot(diffData, aes(x = log2FoldChange, y = p)) +
  geom_point(aes(color = type, size = p), show.legend = F) +
  geom_hline(yintercept = -log10(0.05), color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = -1, color = "gray60", linetype = "dashed") +
  geom_text_repel(
    data = labelGene, aes(label = gene),
    size = 3, fontface = 3,
    nudge_x = .5, nudge_y = .5) +
  scale_radius(range = c(.1, 2)) +
  scale_color_manual(values = typeColor) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  labs(
    x = TeX("$log_{2}(Fold\\,Change)$"),
    y = TeX("$-log_{10}(\\textit{P}\\,value)$")) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line())+labs(title= "")
ggsave("炎性红细胞_volcano.pdf",width = 7,height = 7)






fpkm<-read.csv(file="ternimalE_fpkm.csv",row.names = 1)
head(fpkm)

data <-fpkm
data<-data[which(rowSums(data) > 0),]#去掉全为零的行 情况

newdata<-data[c(rownames((diff_gene))),]
newdata<-data[c(rownames((diff_gene_up))),]
newdata<-data[c(rownames((diff_gene_down))),]

newdata<-na.omit(newdata)
# data<-newdata
# data<-data[which(rowSums(data) > 0),]#去掉全为零的行 情况

library(pheatmap)
library(ggplot2)
library(colorRamps)
library(RColorBrewer)
library(viridis)
library(cowplot)

pheatmap::pheatmap(newdata,scale = "row",
                   cluster_rows = T,
                   cluster_cols = F)

pheatmap::pheatmap(newdata,scale = "row",clustering_method = "ward.D2",
                   cluster_rows = T,show_rownames = F,
                   cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap(
  mat   = newdata,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  scale = "row",
  cluster_rows = T,
  cluster_cols = F,
  border_color = NA,
  #annotation_col = annotation_col,
  #annotation_row = annotation_row,
  show_colnames = TRUE,
  show_rownames = TRUE,
  drop_levels   = TRUE,
  fontsize  = 8#,filename = "fig2h.pdf"
)
pheatmap(
  mat   = newdata,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  scale = "row",
  cluster_rows = T,
  cluster_cols = F,
  border_color = NA,
  #annotation_col = annotation_col,
  #annotation_row = annotation_row,
  show_colnames = TRUE,
  show_rownames = F,
  drop_levels   = TRUE,
  fontsize  = 8#,filename = "fig2h.pdf"
)


pheatmap(
  mat   = newdata,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  scale = "row",clustering_method = "average",
  cluster_rows = T,
  cluster_cols = F,
  border_color = NA,
  #annotation_col = annotation_col,
  #annotation_row = annotation_row,
  show_colnames = TRUE,
  show_rownames = F,
  drop_levels   = TRUE,
  fontsize  = 8#,filename = "fig2h.pdf"
)




resname<-c("immune_ery_vs_non_immune_ery")
# 
# library(clusterProfiler)
# library("org.Mm.eg.db")
# library(ggplot2)
# library(DO.db)
# b<-as.data.frame(rownames(diff_gene))
# eg = bitr(b$`rownames(diff_gene)`, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
# gene <- eg[,2]
# head(gene)
# #分子功能(MolecularFunction)
# 
# ego <- enrichGO(
#   
#   gene          = gene,
#   
#   keyType = "ENTREZID",
#   
#   OrgDb         = org.Mm.eg.db,
#   
#   ont           = "MF",
#   
#   pAdjustMethod = "BH",
#   
#   pvalueCutoff  = 0.05,
#   
#   qvalueCutoff  = 0.05,
#   
#   readable      = TRUE)
# 
# barplot(ego, showCategory =50)
# ggsave(paste0("barplot_GO_MF_all_",resname,".pdf"),width=20,height=20)
# 
# dotplot(ego, showCategory =50)
# ggsave(paste0("dotplot_GO_MF_all_",resname,".pdf"),width=20,height=20)
# write.csv(ego,file=paste0("GO_MF_all_",resname,".csv"))
# 
# 
# 
# #生物过程(biologicalprocess)
# 
# ego <- enrichGO(
#   
#   gene          = gene,
#   
#   keyType = "ENTREZID",
#   
#   OrgDb         = org.Mm.eg.db,
#   
#   ont           = "BP",
#   
#   pAdjustMethod = "BH",
#   
#   pvalueCutoff  = 0.05,
#   
#   qvalueCutoff  = 0.05,
#   
#   readable      = TRUE)
# 
# barplot(ego, showCategory =50)
# ggsave(paste0("barplot_GO_BP_all_",resname,".pdf"),width=20,height=20)
# 
# dotplot(ego, showCategory =50)
# ggsave(paste0("dotplot_GO_BP_all_",resname,".pdf"),width=20,height=20)
# write.csv(ego,file=paste0("GO_BP_all_",resname,".csv"))
# 
# 
# #细胞组成(cellularcomponent)
# 
# ego <- enrichGO(
#   
#   gene          = gene,
#   
#   keyType = "ENTREZID",
#   
#   OrgDb         = org.Mm.eg.db,
#   
#   ont           = "CC",
#   
#   pAdjustMethod = "BH",
#   
#   pvalueCutoff  = 0.05,
#   
#   qvalueCutoff  = 0.05,
#   
#   readable      = TRUE)
# barplot(ego, showCategory =50)
# ggsave(paste0("barplot_GO_CC_all_",resname,".pdf"),width=20,height=20)
# 
# dotplot(ego, showCategory =50)
# ggsave(paste0("dotplot_GO_CC_all_",resname,".pdf"),width=20,height=20)
# write.csv(ego,file=paste0("GO_CC_all_",resname,".csv"))
# 
# 
# 
# ekegg <- enrichKEGG(
#   
#   gene          = gene,
#   
#   keyType     = "kegg",
#   
#   organism   = "mmu",
#   
#   pvalueCutoff      = 0.05,
#   
#   pAdjustMethod     = "BH",
#   
#   qvalueCutoff  = 0.05
#   
# )
# 
# 
# barplot(ekegg, showCategory =50)
# ggsave(paste0("barplot_KEGG_all_",resname,".pdf"),width=20,height=20)
# 
# dotplot(ekegg, showCategory =50)
# ggsave(paste0("dotplot_KEGG_all_",resname,".pdf"),width=20,height=20)
# write.csv(ekegg,file=paste0("KEGG_all_",resname,".csv"))
# 
# 
# ego <- enrichGO(gene = gene,
#                 OrgDb = org.Mm.eg.db, 
#                 pvalueCutoff =0.05, 
#                 qvalueCutoff = 0.05,
#                 ont="all",
#                 readable =T)
# 
# 
# barplot(ego, showCategory =50)
# ggsave(paste0("barplot_GO_MFBPCC_all_",resname,".pdf"),width=20,height=20)
# barplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
# ggsave(paste0("barplot_GO_MFBPCC_all_facet_grid_",resname,".pdf"),width=20,height=20)
# dotplot(ego, showCategory =50)
# ggsave(paste0("dotplot_GO_MFBPCC_all_",resname,".pdf"),width=20,height=20)
# dotplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
# ggsave(paste0("dotplot_GO_MFBPCC_all_facet_grid_",resname,".pdf"),width=20,height=20)
# write.csv(ego,file=paste0("GO_MFBPCC_all_",resname,".csv"))
# # 


#下调基因
library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)
b<-as.data.frame(rownames(diff_gene_down))
eg = bitr(b$`rownames(diff_gene_down)`, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "MF",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)

# barplot(ego, showCategory =50)
# ggsave(paste0("barplot_GO_MF_down_",resname,".pdf"),width=20,height=20)

dotplot(ego, showCategory =50)
ggsave(paste0("dotplot_GO_MF_down_",resname,".pdf"),width=20,height=20)
write.csv(ego,file=paste0("GO_MF_down_",resname,".csv"))



#生物过程(biologicalprocess)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "BP",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)

# barplot(ego, showCategory =50)
# ggsave(paste0("barplot_GO_BP_down_",resname,".pdf"),width=20,height=20)

dotplot(ego, showCategory =50)
ggsave(paste0("dotplot_GO_BP_down_",resname,".pdf"),width=20,height=20)
write.csv(ego,file=paste0("GO_BP_down_",resname,".csv"))


#细胞组成(cellularcomponent)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "CC",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)
# barplot(ego, showCategory =50)
# ggsave(paste0("barplot_GO_CC_down_",resname,".pdf"),width=20,height=20)

dotplot(ego, showCategory =50)
ggsave(paste0("dotplot_GO_CC_down_",resname,".pdf"),width=20,height=20)
write.csv(ego,file=paste0("GO_CC_down_",resname,".csv"))



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)


# barplot(ekegg, showCategory =50)
# ggsave(paste0("barplot_KEGG_down_",resname,".pdf"),width=20,height=20)
ekegg1<-setReadable(ekegg, OrgDb = "org.Mm.eg.db", keyType="ENTREZID")

dotplot(ekegg1, showCategory =50)
ggsave(paste0("dotplot_KEGG_down_",resname,".pdf"),width=20,height=20)
write.csv(ekegg1,file=paste0("KEGG_down_",resname,".csv"))


ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)


# barplot(ego, showCategory =50)
# ggsave(paste0("barplot_GO_MFBPCC_down_",resname,".pdf"),width=20,height=20)
# barplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
# ggsave(paste0("barplot_GO_MFBPCC_down_facet_grid_",resname,".pdf"),width=20,height=20)
dotplot(ego, showCategory =50)
ggsave(paste0("dotplot_GO_MFBPCC_down_",resname,".pdf"),width=20,height=20)
dotplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave(paste0("dotplot_GO_MFBPCC_down_facet_grid_",resname,".pdf"),width=20,height=20)
write.csv(ego,file=paste0("GO_MFBPCC_down_",resname,".csv"))
# 

#上调基因
library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)
b<-as.data.frame(rownames(diff_gene_up))
eg = bitr(b$`rownames(diff_gene_up)`, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(gene)
#分子功能(MolecularFunction)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "MF",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)

# barplot(ego, showCategory =50)
# ggsave(paste0("barplot_GO_MF_up_",resname,".pdf"),width=20,height=20)

dotplot(ego, showCategory =50)
ggsave(paste0("dotplot_GO_MF_up_",resname,".pdf"),width=20,height=20)
write.csv(ego,file=paste0("GO_MF_up_",resname,".csv"))



#生物过程(biologicalprocess)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "BP",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)

# barplot(ego, showCategory =50)
# ggsave(paste0("barplot_GO_BP_up_",resname,".pdf"),width=20,height=20)

dotplot(ego, showCategory =50)
ggsave(paste0("dotplot_GO_BP_up_",resname,".pdf"),width=20,height=20)
write.csv(ego,file=paste0("GO_BP_up_",resname,".csv"))


#细胞组成(cellularcomponent)

ego <- enrichGO(
  
  gene          = gene,
  
  keyType = "ENTREZID",
  
  OrgDb         = org.Mm.eg.db,
  
  ont           = "CC",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff  = 0.05,
  
  qvalueCutoff  = 0.05,
  
  readable      = TRUE)
# barplot(ego, showCategory =50)
# ggsave(paste0("barplot_GO_CC_up_",resname,".pdf"),width=20,height=20)

dotplot(ego, showCategory =50)
ggsave(paste0("dotplot_GO_CC_up_",resname,".pdf"),width=20,height=20)
write.csv(ego,file=paste0("GO_CC_up_",resname,".csv"))



ekegg <- enrichKEGG(
  
  gene          = gene,
  
  keyType     = "kegg",
  
  organism   = "mmu",
  
  pvalueCutoff      = 0.05,
  
  pAdjustMethod     = "BH",
  
  qvalueCutoff  = 0.05
  
)


# barplot(ekegg, showCategory =50)
# ggsave(paste0("barplot_KEGG_up_",resname,".pdf"),width=20,height=20)
ekegg1<-setReadable(ekegg, OrgDb = "org.Mm.eg.db", keyType="ENTREZID")

dotplot(ekegg1, showCategory =50)
ggsave(paste0("dotplot_KEGG_up_",resname,".pdf"),width=20,height=20)
write.csv(ekegg1,file=paste0("KEGG_up_",resname,".csv"))


ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)

# 
# barplot(ego, showCategory =50)
# ggsave(paste0("barplot_GO_MFBPCC_up_",resname,".pdf"),width=20,height=20)
# barplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
# ggsave(paste0("barplot_GO_MFBPCC_up_facet_grid_",resname,".pdf"),width=20,height=20)
dotplot(ego, showCategory =50)
ggsave(paste0("dotplot_GO_MFBPCC_up_",resname,".pdf"),width=20,height=20)
dotplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
ggsave(paste0("dotplot_GO_MFBPCC_up_facet_grid_",resname,".pdf"),width=20,height=20)
write.csv(ego,file=paste0("GO_MFBPCC_up_",resname,".csv"))
