library(stringr)
library(clusterProfiler)
#小鼠的是org.Mm.eg.db
library("org.Mm.eg.db")
library("enrichplot")
library("ggplot2")
library(dplyr)
library(msigdbr)
#各个细胞类型差异分析##############
setwd('D:\\workspace\\KOmouse\\Zbtb18Eybro\\HIP\\enrichment')
DEGs <- read.table('D:\\workspace\\KOmouse\\Zbtb18Eybro\\HIP\\DEGs\\DEGfilt_hipfinal.txt',
                   header=T,check.names=F,sep='\t')
lfc <- 0.1
#上调
for (i in unique(DEGs$Celltype)){
  tryCatch({
    gene <- DEGs[DEGs$Celltype==i & DEGs$avg_log2FC>lfc,'gene'] %>% unique(.)
    if (length(gene)>0){
      geneID <- mapIds(org.Mm.eg.db,keys = gene,column = 'ENTREZID',
                       keytype = 'SYMBOL',multiVals='filter')
      kkbp <- enrichGO(gene = gene,
                       OrgDb = org.Mm.eg.db, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff = 0.2,
                       ont="BP",
                       readable =F,
                       keyType = 'SYMBOL'
      ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min)
      
      kkcc <- enrichGO(gene = gene,
                       OrgDb = org.Mm.eg.db, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff = 0.2,
                       ont="CC",
                       readable =F,
                       keyType = 'SYMBOL'
      ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min)
      
      kkmf <- enrichGO(gene = gene,
                       OrgDb = org.Mm.eg.db, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff = 0.2,
                       ont="MF",
                       readable =F,
                       keyType = 'SYMBOL'
      ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min) 
      
      if (nrow(kkbp)>0 & nrow(kkmf)>0 & nrow(kkcc)>0){
        kkbp <- data.frame(kkbp,check.names = F)
        kkbp$Category <- 'BP'
        kkcc <-  data.frame(kkcc,check.names = F)
        kkcc$Category <- 'CC'
        kkmf <-  data.frame(kkmf,check.names = F)
        kkmf$Category <- 'MF'
        kk4 <- rbind(kkbp,kkcc,kkmf)
        write.table(kk4,file=paste0(i,"_","upDEGenrichGO.txt"),
                    quote=F,row.names = F,sep = '\t')
      } else {
        kk <- enrichGO(gene = gene,
                       OrgDb = org.Mm.eg.db, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff = 0.2,
                       ont="ALL",
                       readable =F,
                       keyType = 'SYMBOL'
        )
        write.table(kk,file=paste0(i,"_","upDEGenrichGO.txt"),
                    quote=F,row.names = F,sep = '\t')
      }
      
      #KEGG
      k <- enrichKEGG(gene = geneID, organism = "mmu",
                      pvalueCutoff =0.05, qvalueCutoff =0.2,
                      keyType = 'kegg',minGSSize=5) %>%
        setReadable(., OrgDb = org.Mm.eg.db, keyType="ENTREZID") %>%
        data.frame(.,check.names = F)
      if('Apoptosis - Mus musculus (house mouse)' %in% k$Description){
        print(paste0('up',i))
      }
      write.table(k,file=paste0(i,"_upDEGenrichKEGG.txt"),
                  quote=F,row.names = F,sep = '\t')
    }
  },error=function(e){})
}

#下调
for (i in unique(DEGs$Celltype)){
  tryCatch({
    gene <- DEGs[DEGs$Celltype==i & DEGs$avg_log2FC<(-lfc),'gene'] %>% unique(.)
    if (length(gene)>0){
      geneID <- mapIds(org.Mm.eg.db,keys = gene,column = 'ENTREZID',
                       keytype = 'SYMBOL',multiVals='filter')
      kkbp <- enrichGO(gene = gene,
                       OrgDb = org.Mm.eg.db, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff = 0.2,
                       ont="BP",
                       readable =F,
                       keyType = 'SYMBOL'
      ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min)
      
      kkcc <- enrichGO(gene = gene,
                       OrgDb = org.Mm.eg.db, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff = 0.2,
                       ont="CC",
                       readable =F,
                       keyType = 'SYMBOL'
      ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min)
      
      kkmf <- enrichGO(gene = gene,
                       OrgDb = org.Mm.eg.db, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff = 0.2,
                       ont="MF",
                       readable =F,
                       keyType = 'SYMBOL'
      ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min) 
      
      if (nrow(kkbp)>0 & nrow(kkmf)>0 & nrow(kkcc)>0){
        kkbp <- data.frame(kkbp,check.names = F)
        kkbp$Category <- 'BP'
        kkcc <-  data.frame(kkcc,check.names = F)
        kkcc$Category <- 'CC'
        kkmf <-  data.frame(kkmf,check.names = F)
        kkmf$Category <- 'MF'
        kk4 <- rbind(kkbp,kkcc,kkmf)
        write.table(kk4,file=paste0(i,"_","downDEGenrichGO.txt"),
                    quote=F,row.names = F,sep = '\t')
      } else {
        kk <- enrichGO(gene = gene,
                       OrgDb = org.Mm.eg.db, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff = 0.2,
                       ont="ALL",
                       readable =F,
                       keyType = 'SYMBOL'
        )
        write.table(kk,file=paste0(i,"_","downDEGenrichGO.txt"),
                    quote=F,row.names = F,sep = '\t')
      }
      #KEGG
      k <- enrichKEGG(gene = geneID, organism = "mmu",
                      pvalueCutoff =0.05, qvalueCutoff =0.2,
                      keyType = 'kegg',minGSSize=5) %>%
        setReadable(., OrgDb = org.Mm.eg.db, keyType="ENTREZID") %>%
        data.frame(.,check.names = F)
      if('Apoptosis - Mus musculus (house mouse)' %in% k$Description){
        print(paste0('up',i))
      }
      write.table(k,file=paste0(i,"_downDEGenrichKEGG.txt"),
                  quote=F,row.names = F,sep = '\t')
    }
  },error=function(e){})
}


#可视化####
source('D:/workspace/KOmouse/Zbtb18brain/enrichmentSourceV3.R')
#DG
boxdata <- read.delim2('DG Pyramidal neurons_downDEGenrichGO-selected.txt',header = T,sep='\t',check.names = F)
p1 <- GOvisual(boxdata,10,0,50,'#836FFF')
boxdata <- read.delim2('DG Pyramidal neurons_upDEGenrichGO-selected.txt',header = T,sep='\t',check.names = F)
p2 <- GOvisual(boxdata,10,0,50,'#FF6347')
#GO合并上下调
p3 <- p1+p2+plot_layout(ncol=1,heights = c(1,1))
ggsave(p3,filename = "DG Pyramidal neurons_enrichGO.pdf",width = 14,height = 22)

#条形图，KEGG
library(scales)
allenrich <- read.table('DG Pyramidal neurons_downDEGenrichKEGG-selected.txt',header=T,check.names=F,sep='\t')
p1 <- keggVisual(allenrich,20,'#836FFF',wrapwith = 30)+
  theme(axis.text.y = element_text(size=34))
allenrich <- read.table('DG Pyramidal neurons_upDEGenrichKEGG.txt',header=T,check.names=F,sep='\t')
p2 <- keggVisual(allenrich,20,'#FF6347',wrapwith = 30)+
  theme(axis.text.y = element_text(size=34))
p3 <- p1+p2+plot_layout(ncol=1,heights = c(1.7,1))#哪个多，哪个大
ggsave(p3,filename = "DG Pyramidal neurons_enrichKEGG.pdf",width = 13,height = 28)

#Intermediate progenitors
boxdata <- read.delim2('Intermediate progenitors_downDEGenrichGO.txt',header = T,sep='\t',check.names = F)
p1 <- GOvisual(boxdata,10,0,0.13,'#836FFF')
boxdata <- read.delim2('Intermediate progenitors_upDEGenrichGO.txt',header = T,sep='\t',check.names = F)
p2 <- GOvisual(boxdata,10,0,0.13,'#FF6347')
#GO合并上下调
p3 <- p1+p2+plot_layout(ncol=1,heights = c(1,1))
ggsave(p3,filename = "Intermediate progenitors_enrichGO.pdf",width = 14,height = 20)

#条形图，KEGG
library(scales)
allenrich <- read.table('Intermediate progenitors_downDEGenrichKEGG.txt',header=T,check.names=F,sep='\t')
p1 <- keggVisual(allenrich,20,'#836FFF')
ggsave(p1,filename = "Intermediate progenitors_enrichKEGG.pdf",width = 17,height = 6)


#Cycling glial cells
boxdata <- read.delim2('Cycling glial cells_downDEGenrichGO.txt',header = T,sep='\t',check.names = F)
p1 <- GOvisual(boxdata,10,0,0.13,'#836FFF')
boxdata <- read.delim2('Cycling glial cells_upDEGenrichGO.txt',header = T,sep='\t',check.names = F)
p2 <- GOvisual(boxdata,10,0,0.2,'#FF6347')
#GO合并上下调
p3 <- p1+p2+plot_layout(ncol=1,heights = c(1,1))
ggsave(p3,filename = "Cycling glial cells_enrichGO.pdf",width = 14,height = 22)

#条形图，KEGG
library(scales)
allenrich <- read.table('Cycling glial cells_downDEGenrichKEGG.txt',header=T,check.names=F,sep='\t')
p1 <- keggVisual(allenrich,20,'#836FFF')
ggsave(p1,filename = "Cycling glial cells_enrichKEGG.pdf",width = 17,height = 3)


#CA1 Pyramidal neurons
boxdata <- read.delim2('CA1 Pyramidal neurons_downDEGenrichGO.txt',header = T,sep='\t',check.names = F)
p1 <- GOvisual(boxdata,10,0,0.12,'#836FFF')
boxdata <- read.delim2('CA1 Pyramidal neurons_upDEGenrichGO.txt',header = T,sep='\t',check.names = F)
p2 <- GOvisual(boxdata,10,0,0.1,'#FF6347')
#GO合并上下调
p3 <- p1+p2+plot_layout(ncol=1,heights = c(1,1))
ggsave(p3,filename = "CA1 Pyramidal neurons_enrichGO.pdf",width = 14,height = 22)

#条形图，KEGG
library(scales)
allenrich <- read.table('CA1 Pyramidal neurons_downDEGenrichKEGG.txt',header=T,check.names=F,sep='\t')
p1 <- keggVisual(allenrich,20,'#836FFF')
allenrich <- read.table('CA1 Pyramidal neurons_upDEGenrichKEGG.txt',header=T,check.names=F,sep='\t')
p2 <- keggVisual(allenrich,20,'#FF6347')
p3 <- p1+p2+plot_layout(ncol=1,heights = c(2,1))#哪个多，哪个大
ggsave(p3,filename = "CA1 Pyramidal neurons_enrichKEGG.pdf",width = 17,height = 20)


#GABAergic neurons
boxdata <- read.delim2('GABAergic neurons_downDEGenrichGO.txt',header = T,sep='\t',check.names = F)
p1 <- GOvisual(boxdata,10,0,0.12,'#836FFF')
boxdata <- read.delim2('GABAergic neurons_upDEGenrichGO.txt',header = T,sep='\t',check.names = F)
p2 <- GOvisual(boxdata,10,0,0.2,'#FF6347')
#GO合并上下调
p3 <- p1+p2+plot_layout(ncol=1,heights = c(1,1))
ggsave(p3,filename = "GABAergic neurons_enrichGO.pdf",width = 14,height = 22)

#条形图，KEGG
library(scales)
allenrich <- read.table('GABAergic neurons_downDEGenrichKEGG.txt',header=T,check.names=F,sep='\t')
p1 <- keggVisual(allenrich,20,'#836FFF')
allenrich <- read.table('GABAergic neurons_upDEGenrichKEGG.txt',header=T,check.names=F,sep='\t')
p2 <- keggVisual(allenrich,20,'#FF6347')
p3 <- p1+p2+plot_layout(ncol=1,heights = c(2,1))#哪个多，哪个大
ggsave(p3,filename = "GABAergic neurons_enrichKEGG.pdf",width = 17,height = 20)


#CA2,3 Pyramidal neurons
boxdata <- read.delim2('CA2,3 Pyramidal neurons_downDEGenrichGO.txt',header = T,sep='\t',check.names = F)
p1 <- GOvisual(boxdata,10,0,0.1,'#836FFF')
boxdata <- read.delim2('CA2,3 Pyramidal neurons_upDEGenrichGO.txt',header = T,sep='\t',check.names = F)
p2 <- GOvisual(boxdata,10,0,0.1,'#FF6347')
#GO合并上下调
p3 <- p1+p2+plot_layout(ncol=1,heights = c(1,1))
ggsave(p3,filename = "CA2,3 Pyramidal neurons_enrichGO.pdf",width = 14,height = 22)

#条形图，KEGG
library(scales)
allenrich <- read.table('CA2,3 Pyramidal neurons_downDEGenrichKEGG.txt',header=T,check.names=F,sep='\t')
p1 <- keggVisual(allenrich,20,'#836FFF')
allenrich <- read.table('CA2,3 Pyramidal neurons_upDEGenrichKEGG.txt',header=T,check.names=F,sep='\t')
p2 <- keggVisual(allenrich,20,'#FF6347')
p3 <- p1+p2+plot_layout(ncol=1,heights = c(1.5,1))#哪个多，哪个大
ggsave(p3,filename = "CA2,3 Pyramidal neurons_enrichKEGG.pdf",width = 17,height = 20)


#Interneurons
boxdata <- read.delim2('Interneurons_downDEGenrichGO.txt',header = T,sep='\t',check.names = F)
p1 <- GOvisual(boxdata,10,0,0.25,'#836FFF')
boxdata <- read.delim2('Interneurons_upDEGenrichGO.txt',header = T,sep='\t',check.names = F)
p2 <- GOvisual(boxdata,10,0,0.25,'#FF6347')
#GO合并上下调
p3 <- p1+p2+plot_layout(ncol=1,heights = c(1,1))
ggsave(p3,filename = "Interneurons_enrichGO.pdf",width = 14,height = 22)

#条形图，KEGG
library(scales)
allenrich <- read.table('Interneurons_downDEGenrichKEGG.txt',header=T,check.names=F,sep='\t')
p1 <- keggVisual(allenrich,20,'#836FFF')
allenrich <- read.table('Interneurons_upDEGenrichKEGG.txt',header=T,check.names=F,sep='\t')
p2 <- keggVisual(allenrich,20,'#FF6347')
p3 <- p1+p2+plot_layout(ncol=1,heights = c(2,1))#哪个多，哪个大
ggsave(p3,filename = "Interneurons_enrichKEGG.pdf",width = 17,height = 20)


#CA2,3 Pyramidal neurons
boxdata <- read.delim2('Astrocytes_downDEGenrichGO.txt',header = T,sep='\t',check.names = F)
p1 <- GOvisual(boxdata,10,0,0.13,'#836FFF')
boxdata <- read.delim2('Astrocytes_upDEGenrichGO.txt',header = T,sep='\t',check.names = F)
p2 <- GOvisual(boxdata,10,0,0.1,'#FF6347')
#GO合并上下调
p3 <- p1+p2+plot_layout(ncol=1,heights = c(1,1))
ggsave(p3,filename = "Astrocytes_enrichGO.pdf",width = 14,height = 22)

#条形图，KEGG
library(scales)
allenrich <- read.table('Astrocytes_downDEGenrichKEGG.txt',header=T,check.names=F,sep='\t')
p1 <- keggVisual(allenrich,20,'#836FFF')
ggsave(p1,filename = "Astrocytes_enrichKEGG.pdf",width = 17,height = 12)

#合并细胞类型的富集分析################
setwd('D:\\workspace\\KOmouse\\Zbtb18Eybro\\HIP\\enrichment\\Mergetp/')
load('D:\\workspace\\KOmouse\\Zbtb18Eybro\\HIP\\DEGs\\MergetpDEGfilt.rda')
lfc <- 0.1
#上调
for (i in unique(MergetpDEGfilt$NewType)){
  tryCatch({
    gene <- MergetpDEGfilt[MergetpDEGfilt$NewType==i & MergetpDEGfilt$avg_log2FC>lfc,'gene'] %>% unique(.)
    if (length(gene)>0){
      geneID <- mapIds(org.Mm.eg.db,keys = gene,column = 'ENTREZID',
                       keytype = 'SYMBOL',multiVals='filter')
      kkbp <- enrichGO(gene = gene,
                       OrgDb = org.Mm.eg.db, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff = 0.2,
                       ont="BP",
                       readable =F,
                       keyType = 'SYMBOL'
      ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min)
      
      kkcc <- enrichGO(gene = gene,
                       OrgDb = org.Mm.eg.db, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff = 0.2,
                       ont="CC",
                       readable =F,
                       keyType = 'SYMBOL'
      ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min)
      
      kkmf <- enrichGO(gene = gene,
                       OrgDb = org.Mm.eg.db, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff = 0.2,
                       ont="MF",
                       readable =F,
                       keyType = 'SYMBOL'
      ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min) 
      
      if (nrow(kkbp)>0 & nrow(kkmf)>0 & nrow(kkcc)>0){
        kkbp <- data.frame(kkbp,check.names = F)
        kkbp$Category <- 'BP'
        kkcc <-  data.frame(kkcc,check.names = F)
        kkcc$Category <- 'CC'
        kkmf <-  data.frame(kkmf,check.names = F)
        kkmf$Category <- 'MF'
        kk4 <- rbind(kkbp,kkcc,kkmf)
        write.table(kk4,file=paste0(i,"_","upDEGenrichGO.txt"),
                    quote=F,row.names = F,sep = '\t')
      } else {
        kk <- enrichGO(gene = gene,
                       OrgDb = org.Mm.eg.db, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff = 0.2,
                       ont="ALL",
                       readable =F,
                       keyType = 'SYMBOL'
        )
        write.table(kk,file=paste0(i,"_","upDEGenrichGO.txt"),
                    quote=F,row.names = F,sep = '\t')
      }
      
      #KEGG
      k <- enrichKEGG(gene = geneID, organism = "mmu",
                      pvalueCutoff =0.05, qvalueCutoff =0.2,
                      keyType = 'kegg',minGSSize=5) %>%
        setReadable(., OrgDb = org.Mm.eg.db, keyType="ENTREZID") %>%
        data.frame(.,check.names = F)
      if('Apoptosis - Mus musculus (house mouse)' %in% k$Description){
        print(paste0('up',i))
      }
      write.table(k,file=paste0(i,"_upDEGenrichKEGG.txt"),
                  quote=F,row.names = F,sep = '\t')
    }
  },error=function(e){})
}

#下调
for (i in unique(MergetpDEGfilt$NewType)){
  tryCatch({
    gene <- MergetpDEGfilt[MergetpDEGfilt$NewType==i & MergetpDEGfilt$avg_log2FC<(-lfc),'gene'] %>% unique(.)
    if (length(gene)>0){
      geneID <- mapIds(org.Mm.eg.db,keys = gene,column = 'ENTREZID',
                       keytype = 'SYMBOL',multiVals='filter')
      kkbp <- enrichGO(gene = gene,
                       OrgDb = org.Mm.eg.db, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff = 0.2,
                       ont="BP",
                       readable =F,
                       keyType = 'SYMBOL'
      ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min)
      
      kkcc <- enrichGO(gene = gene,
                       OrgDb = org.Mm.eg.db, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff = 0.2,
                       ont="CC",
                       readable =F,
                       keyType = 'SYMBOL'
      ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min)
      
      kkmf <- enrichGO(gene = gene,
                       OrgDb = org.Mm.eg.db, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff = 0.2,
                       ont="MF",
                       readable =F,
                       keyType = 'SYMBOL'
      ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min) 
      
      if (nrow(kkbp)>0 & nrow(kkmf)>0 & nrow(kkcc)>0){
        kkbp <- data.frame(kkbp,check.names = F)
        kkbp$Category <- 'BP'
        kkcc <-  data.frame(kkcc,check.names = F)
        kkcc$Category <- 'CC'
        kkmf <-  data.frame(kkmf,check.names = F)
        kkmf$Category <- 'MF'
        kk4 <- rbind(kkbp,kkcc,kkmf)
        write.table(kk4,file=paste0(i,"_","downDEGenrichGO.txt"),
                    quote=F,row.names = F,sep = '\t')
      } else {
        kk <- enrichGO(gene = gene,
                       OrgDb = org.Mm.eg.db, 
                       pvalueCutoff =0.05, 
                       qvalueCutoff = 0.2,
                       ont="ALL",
                       readable =F,
                       keyType = 'SYMBOL'
        )
        write.table(kk,file=paste0(i,"_","downDEGenrichGO.txt"),
                    quote=F,row.names = F,sep = '\t')
      }
      #KEGG
      k <- enrichKEGG(gene = geneID, organism = "mmu",
                      pvalueCutoff =0.05, qvalueCutoff =0.2,
                      keyType = 'kegg',minGSSize=5) %>%
        setReadable(., OrgDb = org.Mm.eg.db, keyType="ENTREZID") %>%
        data.frame(.,check.names = F)
      if('Apoptosis - Mus musculus (house mouse)' %in% k$Description){
        print(paste0('up',i))
      }
      write.table(k,file=paste0(i,"_downDEGenrichKEGG.txt"),
                  quote=F,row.names = F,sep = '\t')
    }
  },error=function(e){})
}

#画条形图
#画CA的
source('D:/workspace/KOmouse/Zbtb18brain/enrichmentSourceV3.R')
setwd('D:/workspace/KOmouse/Zbtb18Eybro/HIP/enrichment/Mergetp')
boxdata <- read.delim2('CA Pyramidal neurons_downDEGenrichGO-selected.txt',header = T,sep='\t',check.names = F)
p1 <- GOvisual(boxdata,20,0,130,'#836FFF')
boxdata <- read.delim2('CA Pyramidal neurons_upDEGenrichGO-selected.txt',header = T,sep='\t',check.names = F)
p2 <- GOvisual(boxdata,20,0,40,'#FF6347')
#GO合并上下调
p3 <- p1+p2+plot_layout(ncol=1,heights = c(1,1))
ggsave(p3,filename = "CA Pyramidal neurons_enrichGO.pdf",width = 11,height = 18)

#条形图，KEGG
library(scales)
allenrich <- read.table('CA Pyramidal neurons_downDEGenrichKEGG-selected.txt',header=T,check.names=F,sep='\t')
p1 <- keggVisual(allenrich,20,'#836FFF',wrapwith = 30)+
  theme(axis.text.y = element_text(size=34))
allenrich <- read.table('CA Pyramidal neurons_upDEGenrichKEGG.txt',header=T,check.names=F,sep='\t')
p2 <- keggVisual(allenrich,20,'#FF6347',wrapwith = 30)+
  theme(axis.text.y = element_text(size=34))
p3 <- p1+p2+plot_layout(ncol=1,heights = c(1.1,1))#哪个多，哪个大
ggsave(p3,filename = "CA Pyramidal neurons_enrichKEGG.pdf",width = 11,height = 25)


# #CA DG一致下调富集###########
# DEGs <- read.table('D:\\workspace\\KOmouse\\Zbtb18Eybro\\HIP\\DEGs\\CADGdown.txt',
#                    header=T,check.names=F,sep='\t')
# gene <- DEGs$gene
# geneID <- mapIds(org.Mm.eg.db,keys = gene,column = 'ENTREZID',
#                  keytype = 'SYMBOL',multiVals='filter')
# kkbp <- enrichGO(gene = gene,
#                  OrgDb = org.Mm.eg.db, 
#                  pvalueCutoff =0.05, 
#                  qvalueCutoff = 0.2,
#                  ont="BP",
#                  readable =F,
#                  keyType = 'SYMBOL'
# ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min)
# 
# kkcc <- enrichGO(gene = gene,
#                  OrgDb = org.Mm.eg.db, 
#                  pvalueCutoff =0.05, 
#                  qvalueCutoff = 0.2,
#                  ont="CC",
#                  readable =F,
#                  keyType = 'SYMBOL'
# ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min)
# 
# kkmf <- enrichGO(gene = gene,
#                  OrgDb = org.Mm.eg.db, 
#                  pvalueCutoff =0.05, 
#                  qvalueCutoff = 0.2,
#                  ont="MF",
#                  readable =F,
#                  keyType = 'SYMBOL'
# ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min) 
# 
# if (nrow(kkbp)>0 & nrow(kkmf)>0 & nrow(kkcc)>0){
#   kkbp <- data.frame(kkbp,check.names = F)
#   kkbp$Category <- 'BP'
#   kkcc <-  data.frame(kkcc,check.names = F)
#   kkcc$Category <- 'CC'
#   kkmf <-  data.frame(kkmf,check.names = F)
#   kkmf$Category <- 'MF'
#   kk4 <- rbind(kkbp,kkcc,kkmf)
#   write.table(kk4,file="CADG-down_enrichGO.txt",
#               quote=F,row.names = F,sep = '\t')
# } else {
#   kk <- enrichGO(gene = gene,
#                  OrgDb = org.Mm.eg.db, 
#                  pvalueCutoff =0.05, 
#                  qvalueCutoff = 0.2,
#                  ont="ALL",
#                  readable =F,
#                  keyType = 'SYMBOL'
#   )
#   write.table(kk,file="CADG-down_enrichGO.txt",
#               quote=F,row.names = F,sep = '\t')
# }
# 
# #KEGG
# k <- enrichKEGG(gene = geneID, organism = "mmu",
#                 pvalueCutoff =0.05, qvalueCutoff =0.2,
#                 keyType = 'kegg',minGSSize=5) %>%
#   setReadable(., OrgDb = org.Mm.eg.db, keyType="ENTREZID") %>%
#   data.frame(.,check.names = F)
# write.table(k,file="CADG-down_enrichKEGG.txt",
#             quote=F,row.names = F,sep = '\t')
# 
# 
# #CA DG一致上富集###########
# DEGs <- read.table('D:\\workspace\\KOmouse\\Zbtb18Eybro\\HIP\\DEGs\\CADGup.txt',
#                    header=T,check.names=F,sep='\t')
# gene <- DEGs$gene
# geneID <- mapIds(org.Mm.eg.db,keys = gene,column = 'ENTREZID',
#                  keytype = 'SYMBOL',multiVals='filter')
# kkbp <- enrichGO(gene = gene,
#                  OrgDb = org.Mm.eg.db, 
#                  pvalueCutoff =0.05, 
#                  qvalueCutoff = 0.2,
#                  ont="BP",
#                  readable =F,
#                  keyType = 'SYMBOL'
# ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min)
# 
# kkcc <- enrichGO(gene = gene,
#                  OrgDb = org.Mm.eg.db, 
#                  pvalueCutoff =0.05, 
#                  qvalueCutoff = 0.2,
#                  ont="CC",
#                  readable =F,
#                  keyType = 'SYMBOL'
# ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min)
# 
# kkmf <- enrichGO(gene = gene,
#                  OrgDb = org.Mm.eg.db, 
#                  pvalueCutoff =0.05, 
#                  qvalueCutoff = 0.2,
#                  ont="MF",
#                  readable =F,
#                  keyType = 'SYMBOL'
# ) %>% simplify(.,cutoff=0.7,by='p.adjust',select_fun=min) 
# 
# if (nrow(kkbp)>0 & nrow(kkmf)>0 & nrow(kkcc)>0){
#   kkbp <- data.frame(kkbp,check.names = F)
#   kkbp$Category <- 'BP'
#   kkcc <-  data.frame(kkcc,check.names = F)
#   kkcc$Category <- 'CC'
#   kkmf <-  data.frame(kkmf,check.names = F)
#   kkmf$Category <- 'MF'
#   kk4 <- rbind(kkbp,kkcc,kkmf)
#   write.table(kk4,file="CADG-up_enrichGO.txt",
#               quote=F,row.names = F,sep = '\t')
# } else {
#   kk <- enrichGO(gene = gene,
#                  OrgDb = org.Mm.eg.db, 
#                  pvalueCutoff =0.05, 
#                  qvalueCutoff = 0.2,
#                  ont="ALL",
#                  readable =F,
#                  keyType = 'SYMBOL'
#   )
#   write.table(kk,file="CADG-up_enrichGO.txt",
#               quote=F,row.names = F,sep = '\t')
# }
# 
# #KEGG
# k <- enrichKEGG(gene = geneID, organism = "mmu",
#                 pvalueCutoff =0.05, qvalueCutoff =0.2,
#                 keyType = 'kegg',minGSSize=5) %>%
#   setReadable(., OrgDb = org.Mm.eg.db, keyType="ENTREZID") %>%
#   data.frame(.,check.names = F)
# write.table(k,file="CADG-up_enrichKEGG.txt",
#             quote=F,row.names = F,sep = '\t')
