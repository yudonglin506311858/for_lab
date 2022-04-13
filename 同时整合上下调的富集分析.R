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
setwd("D:/anxiiuli_data/count_炎红/CD44")


data_up<-read.csv("GO_up.csv")
data_down<-read.csv("GO_down.csv")

data_up$change<-"up"
data_down$change<-"down"
data<-rbind(data_up,data_down)
data<-data[,c("Description","pvalue","change")]
library(ggplot2)
library(dplyr)
set.seed(123)

head(data)

data <- data %>% 
  mutate(p2 = ifelse(change == "up", -log10(pvalue), log10(pvalue))) %>% 
  arrange(change,p2) 
data$Description = factor(data$Description,levels = unique(data$Description),ordered = T)
head(data)

tmp = with(data, labeling::extended(range(p2)[1], range(p2)[2], m = 5));tmp
## [1] -3 -2 -1  0  1  2
lm = tmp[c(1,length(tmp))];lm

ggplot(data, aes(x=Description, y=p2)) +
  geom_segment( aes(x=Description, xend=Description, y=0, yend=p2, color=change), size=5, alpha=0.9) +
  theme_light() +
  theme(
    panel.border = element_blank()
  ) +
  xlab("") +
  ylab("-log10(pvalue)")+
  ylim(lm)+
  scale_y_continuous(breaks = tmp,
                     labels = abs(tmp))+
  coord_flip()


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
typeColor <- structure(
  c(pal_nejm()(2)),
  names = c("up", "down")
)
ggplot(data, aes(x=Description, y=p2)) +
  geom_segment( aes(x=Description, xend=Description, y=0, yend=p2, color=change), size=5, alpha=0.9) +
  theme_light() +
  theme(
    panel.border = element_blank()
  ) +
  xlab("") +
  ylab("-log10(pvalue)")+
  ylim(lm)+
  scale_y_continuous(breaks = tmp,
                     labels = abs(tmp))+
  coord_flip()+
  scale_color_manual(values = typeColor)




#参考资料：https://www.jianshu.com/p/b763cbc92fd7
