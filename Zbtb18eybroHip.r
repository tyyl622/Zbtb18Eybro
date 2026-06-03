library(patchwork)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library(tidyr)
library(stringr)
library(stringi)
library(vctrs)
library(cluster)
library(clustree)
library(ggGenshin)
library(Seurat)

samples <- c('HipWt','HipMut')
second <- '/02.count/filter_matrix'
for (sam in samples){
    print(sam)
    data.dir <- paste0(path,sam,second)
    exp <- Read10X(
    data.dir,
    gene.column = 1,
    cell.column = 1,
    unique.features = TRUE,
    strip.suffix = FALSE
    )
    object <- CreateSeuratObject(counts = exp, project = sam,min.cells = 10, min.features = 200)
    Idents(object) <- sam
    object@meta.data$Sample <- sam
    save(object,file=paste0(sam,'_ob.rda')) 
}

#过滤及细胞重命名
#HIP
load('HipWt_ob.rda')
HipWt <- object
HipWt$Type <- 'wt'

load('HipMut_ob.rda')
HipMut <- object
HipMut$Type <- 'mut'

#过滤
filtob <- function(object,nfeamin,nfeamax,nCountmin,nCountmax,permt,Hb){
    obnew <- PercentageFeatureSet(object ,"^mt-", col.name = "percent.mt")
    obnew <- PercentageFeatureSet(obnew ,"^Hb", col.name = "percent.Hb")
    obnew <- PercentageFeatureSet(obnew ,"^Rp[sl]", col.name = "percent.Rpsl")
    obnew <- obnew[,obnew$nFeature_RNA > nfeamin & 
                    obnew$nFeature_RNA < nfeamax &
                    obnew$nCount_RNA > nCountmin &
                    obnew$nCount_RNA < nCountmax &
                    obnew$percent.mt < permt & 
                    obnew$percent.Hb < Hb]
    counts <- GetAssayData(obnew, assay = "RNA")
    de1 <- grep('^mt-',rownames(counts))
    de2 <- grep('^Hb',rownames(counts))
    de3 <- grep('^Rp[sl]',rownames(counts))
    de <- c(de1,de2,de3)
    counts <- counts[-de,]
    obnew <- subset(obnew, features = rownames(counts))
    #修改细胞名字，以便和loupe brower对应
    sampName <- obnew$Sample 
    newname <- paste0(sampName,'_',colnames(obnew))
    obnew <- RenameCells(obnew,new.names=newname)   
    return(obnew)
}
#filtob(H3HIP1,100,7500,200,30000,20,5)
#HIP
HipWt <- filtob(HipWt,100,7500,200,20000,10,5)
HipMut <- filtob(HipMut,100,7500,200,20000,10,5)
save(HipWt,HipMut,file='Object/filtob.rda')

#按照批次，也就是样本，分别进行标准流程,包括质控，要先画一下质控图，再决定nCount_RNA阈值
hipmerge <- merge(HipWt,HipMut)
#HIP
gene_counts_per_cell <- hipmerge$nFeature_RNA
gene_counts_summary <- summary(gene_counts_per_cell)
gene_counts_summary.show <- data.frame(Quantile=names(gene_counts_summary),
                                       summary=as.numeric(gene_counts_summary))
write.table(gene_counts_summary.show,file='HIPgene_counts_summary.show',sep = '\t',quote = F,row.names = F)

#plot
data <- data.frame(gene_counts_per_cell)
colnames(data) <- 'gene_counts_per_cell'
data$Cell <- rownames(data)
png('HIPCellCount.png',width=300,height=200,unit='mm',res=300)
ggplot(data,aes(gene_counts_per_cell))+
       geom_histogram(binwidth = 50,fill='blue',color='blue')+
       ggtitle("Distribution of Gene counts")+
       xlab("Gene counts")+ylab("Cell number")+theme_bw()+
       theme(axis.text=element_text(size=30),
             title=element_text(size=30),
             axis.title=element_text(size=30))+
       theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
dev.off()

#pcor
Idents(hipmerge) <- hipmerge$Sample
# cols4 <- c("DarkTurquoise",'Gold',"red", "DodgerBlue3")
p <- VlnPlot(hipmerge, features = c("nFeature_RNA", "nCount_RNA",
                                "percent.Hb",'percent.mt','percent.Rpsl'), 
                                ncol = 3,pt.size = 0)
p1 <- lapply(list(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]]),function(x){
      return(x+theme(axis.title.x=element_blank(),
                    axis.text=element_text(size=20),
                    title=element_text(size=20)))
})        
p2 <- plot_grid(plotlist=p1,ncol=3)
ggsave(p2,filename='merge_VlnPlotForHbMt_hip.pdf',width=15,height=12)

#只画pcor
Idents(hipmerge) <- hipmerge$Sample
pcor <- FeatureScatter(hipmerge,feature1='nFeature_RNA',
                    feature2 = 'nCount_RNA')+
                    theme(axis.title.x=element_blank(),
                    axis.text=element_text(size=20),
                    axis.text.x=element_text(angle=45,hjust=1,vjust=1),
                    title=element_text(size=20),
                    legend.title=element_text(size=20),
                    legend.text=element_text(size=20))
ggsave(pcor,filename='pcorhip.pdf',width=6.5,height=6)

##标准化
DefaultAssay(hipmerge) <- 'RNA'
hipmerge <- NormalizeData(hipmerge)
hipmerge <- FindVariableFeatures(hipmerge,nfeatures=3000)
#处理高变基因
VariableFeatures <- VariableFeatures(hipmerge)[1:3000]
#个性化scale
hipmerge <- ScaleData(hipmerge, features = VariableFeatures,vars.to.regress = c("percent.mt", "percent.Rpsl", "percent.Hb"))
hipmerge <- RunPCA(hipmerge)
save(hipmerge,file='hipmerge.rda')

#尝试用各种方法做整合，注意，整合数据是为了聚类，不是做基因的批次校正，因为只挑选HVG基因做整合
options(future.globals.maxSize = 8000 * 1024^2)#设大一点，否则会报错
hipInte <- IntegrateLayers(
  object = hipmerge, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
#全部整合模型运行完之后，再合并layers,才能做findallmarker
print('Joint')
hipInte[["RNA"]] <- JoinLayers(hipInte[["RNA"]])
save(hipInte,file='hipInteJoint.rda')

#library size 的图Joint合并后才能画
data <- data.frame(Matrix::colSums(hipInte))
colnames(data) <- 'library_size'
data$Cell <- rownames(data)
png('distri_lib_size_hip.png',width=300,height=200,unit='mm',res=300)
ggplot(data,aes(library_size))+
       geom_histogram(binwidth = 50,fill='blue',color='blue')+
       ggtitle("Distribution of library size")+
       xlab("library size")+ylab("Cell number")+theme_bw()+
       theme(axis.text=element_text(size=30),
             title=element_text(size=30),
             axis.title=element_text(size=30))+
       theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
dev.off()

#循环resolution，可以大一点，然后每个小类再聚亚类
hipfinal <- FindNeighbors(hipInte, reduction = 'harmony', dims = 1:30)
hipfinal <- FindClusters(hipfinal, resolution = seq(0.04,0.1,0.01))
hipfinal <- RunUMAP(hipfinal, reduction = 'harmony', dims = 1:30, reduction.name = 'harmony_umap')
hipfinal <- RunTSNE(hipfinal, reduction = 'harmony', dims = 1:30, reduction.name = 'harmony_tsne')
save(hipfinal,file='hipfinal_harmony.rda')
hipmeta <- hipfinal@meta.data
hipmeta$Cell <- rownames(hipmeta)
hipexp <- hipfinal@assays$RNA$data

#降维可视化
myreso <- seq(0.04,0.1,0.01)
for (reso in myreso){
    cluster.name <- paste0('RNA_snn_res.',reso)
    p1 <- DimPlot(
    hipfinal,label=TRUE,repel=FALSE,
    reduction = 'harmony_tsne',
    group.by = c("Sample", cluster.name),
    combine = TRUE, label.size = 2
    )
    ggsave(p1,filename=paste0('reso',reso,'_inteTSNE.png'),width=15,height=6)

    p2 <- DimPlot(
    hipfinal,label=TRUE,repel=FALSE,
    reduction = 'harmony_umap',
    group.by = c("Sample", cluster.name),
    combine = TRUE, label.size = 2
    )
    ggsave(p2,filename=paste0('reso',reso,'_inteUMAP.png'),width=15,height=6)
}

#allmarer##############
load('hipfinal_harmony.rda')
reso <- 0.25
Idents(hipfinal) <- hipfinal@meta.data[,'RNA_snn_res.0.25']
allmarkers <- FindAllMarkers(hipfinal, only.pos = FALSE, min.pct = 0.2, 
                        logfc.threshold = 0.1,assay='RNA',slot='data')
save(allmarkers,file=paste0('allmarkerharmony.reso',reso,'.rda'))
allmarkerfilt <- allmarkers[allmarkers$p_val_adj<0.05,]
save(allmarkerfilt,file=paste0('allmarkerfiltharmony.reso',reso,'.rda'))
write.csv(allmarkerfilt,file=paste0('allmarkerfiltharmony.reso',reso,'.csv'),row.names=F,quote=T)
table(hipfinal$RNA_snn_res.0.25)

MAKERMAP <- function(myclus,refsub,refs,mymarkerfilt,logfc){
  MAmapsub <- data.frame()
  for(mycl in myclus){
  for (subclass in refsub){
      #subclass marker
      submarker <- refs[refs$cluster==subclass,]
      submarker <- submarker[submarker$avg_log2FC>logfc,'gene']
      #class marker
      mymarker <- mymarkerfilt[mymarkerfilt$cluster==mycl,]

      mymarker <- mymarker[mymarker$avg_log2FC>logfc,'gene']
      #取交集
      intergenes <- intersect(submarker,mymarker)
      numbers <- length(intergenes)
      myper <- numbers/length(mymarker) %>% round(.,2)
      subper <- numbers/length(submarker) %>% round(.,2)
      #取平均logfc
      subfc <- mean(refs[refs$cluster==subclass & refs$gene %in% intergenes,'avg_log2FC'])
      myfc <- mean(mymarkerfilt[mymarkerfilt$cluster==mycl & mymarkerfilt$gene %in% intergenes,'avg_log2FC'])
      tmp <- data.frame(mycluster=mycl,sub=subclass,inter=numbers,
                  myper=myper,subper=subper,
                  mymarker=paste0(mymarker,collapse=','),
                  intergene=paste0(intergenes,collapse=','),
                  outgene=paste0(setdiff(mymarker,intergenes),collapse=','),
                  myfc=myfc,subfc=subfc)
      MAmapsub <- rbind(MAmapsub,tmp)
      }
    }
  return(MAmapsub)
}

#Nature2021 marker map
#先参考Nature2021的结果注释
load('Nature2021cortex_allmarkerfilt.rda')
cortexref <- allmarkerfilt
cortexref <- cortexref[cortexref$avg_log2FC>0,]
cortexrefsub <- as.vector(unique(cortexref$cluster))
cortexrefsub <- cortexrefsub[-c(1,4,8,10,24,25,26,9,19,22)]

load('allmarkerfiltharmony.reso0.25.rda')
hipmyclus <- as.vector(unique(allmarkerfilt$cluster))
hipMarkermap <- MAKERMAPSUB(hipmyclus,refsub=cortexrefsub,refs=cortexref,allmarkerfilt,1)
write.table(hipMarkermap,'hip_reso0.25_MakerMap_t1.txt',sep='\t',row.names=F,quote=F)

#参考成鼠脑注释
load('hipclassallmarkerfilt.rda')
hipclass <- allmarkerfilt
hipclass <- hipclass[hipclass$avg_log2FC>0,]
allenhipcla <- as.vector(unique(hipclass$cluster))

load('hipallmarkerfilt.rda')
hipsubclass <- allmarkerfilt
hipsubclass <- hipsubclass[hipsubclass$avg_log2FC>0,]
allenhipsub <- as.vector(unique(hipsubclass$cluster))
allenhipsub <- allenhipsub[-c(17,18,22,24)]

#allen subclass
load('allmarkerfiltharmony.reso0.25.rda')
hipmyclus <- as.vector(unique(allmarkerfilt$cluster))
hipMarkermap <- MAKERMAP(hipmyclus,refsub=allenhipsub,refs=hipsubclass,allmarkerfilt,1)
write.table(hipMarkermap,'hip_reso0.25_MakerMap_allensubclass_t1.txt',sep='\t',row.names=F,quote=F)

#allen class
load('allmarkerfiltharmony.reso0.25.rda')
hipmyclus <- as.vector(unique(allmarkerfilt$cluster))
hipMarkermap <- MAKERMAP(hipmyclus,refsub=allenhipcla,refs=hipclass,allmarkerfilt,1)
write.table(hipMarkermap,'hip_reso0.25_MakerMap_allenclass_t1.txt',sep='\t',row.names=F,quote=F)

#合并两种注释结果
hipmeta <- hipfinal@meta.data
hipmeta$Celltype <- hipmeta$RNA_snn_res.0.25
hipmeta$Celltype <- gsub('^0$','Intermediate progenitors',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^1$','Astrocytes',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^2$','CA2,3 Pyramidal neurons',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^3$','CA2,3 Pyramidal neurons',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^4$','Inhibitory neurons',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^5$','CA1 Pyramidal neurons',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^6$','DG Pyramidal neurons',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^7$','CA1 Pyramidal neurons',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^8$','Cycling glial cells',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^9$','Inhibitory neurons',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^10$','Inhibitory neurons',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^11$','DG Pyramidal neurons',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^12$','VLMC',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^13$','Inhibitory neurons',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^14$','Endothelial cells',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^15$','Astrocytes',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^16$','Inhibitory neurons',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^17$','Astrocytes',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^18$','Pericytes',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^19$','Microglia',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^20$','Oligodendrocytes',hipmeta$Celltype)
hipmeta$Celltype <- gsub('^21$','Astrocytes',hipmeta$Celltype)
hipmeta$Cell <- rownames(hipmeta)
hipfinal$Celltype <- hipmeta[colnames(hipfinal),'Celltype']
save(hipmeta,file='hipmeta.rda')
save(hipfinal,file='hipfinal.rda')

#cluster0聚亚类,从中获得另外的DG和CA1
reducMark <- function(object,resolutions,reducname,allmarkername){
  for ( reso in resolutions){
    prefix <- paste0('RNA_snn_res.',reso)
    #指定颜色
    colourCount <- length(unique(object@meta.data[,prefix]))
    getPalette = colorRampPalette(brewer.pal(colourCount, "Set2"))
    mycolsub <- getPalette(colourCount)
    mycolsub <- c(mycolsub,'#4DBBD5','Gold')
    names(mycolsub) <- c(as.vector(unique(object@meta.data[,prefix])),'wt','mut')
    #UMAP
    p1 <- DimPlot(object,cols=mycolsub,reduction = 'harmony_umap',group.by = c(prefix,'Type'),
                    combine = FALSE, label.size = 5)
    #TSNE
    p2 <- DimPlot(object,cols=mycolsub,reduction = 'harmony_tsne',group.by = c(prefix,'Type'),
                    combine = FALSE, label.size = 5)
    p <- c(p1,p2)
    p <- aplot::plot_list(gglist=p,ncol=4,nrow=1)
    ggsave(p,filename=paste0(reducname,reso,'_tsne+umap.pdf'),width=20,height=4.5)

    #allmarer
    Idents(object) <- object@meta.data[,prefix]
    allmarkers <- FindAllMarkers(object, only.pos = FALSE, min.pct = 0.2, 
                            logfc.threshold = 0.1,assay='RNA',slot='data')
    allmarkerfilt <- allmarkers[allmarkers$p_val_adj<0.05,]
    save(allmarkerfilt,file=paste0(allmarkername,reso,'.rda'))
    write.csv(allmarkerfilt,file=paste0(allmarkername,reso,'.csv'),row.names=F,quote=T)
  }
}
inps <- subset(hipfinal,Celltype=='Intermediate progenitors')
inps <- FindNeighbors(inps, reduction = 'harmony', dims = 1:30)
inps <- FindClusters(inps, resolution = 0.21)
inps <- RunUMAP(inps, reduction = 'harmony', dims = 1:30,reduction.name='harmony_umap')
inps <- RunTSNE(inps, reduction = 'harmony', dims = 1:30,reduction.name='harmony_tsne')
#tsne umap allmarker
reducname <- 'hipreinps_reso'
allmarkername <- 'hipallmarkerfiltharmony-inps.reso'
reducMark(inps,0.21,reducname,allmarkername)

#marker map
load('Nature2021cortex_allmarkerfilt.rda')
cortexref <- allmarkerfilt
cortexref <- cortexref[cortexref$avg_log2FC>0,]
cortexrefsub <- as.vector(unique(cortexref$cluster))
cortexrefsub <- cortexrefsub[-c(1,4,8,10,24,25,26,9,19,22)]

load('Inpsub/hipallmarkerfiltharmony-inps.reso0.21.rda')
hipmyclus <- as.vector(unique(allmarkerfilt$cluster))
InpsMarkermap <- MAKERMAPSUB(hipmyclus,refsub=cortexrefsub,refs=cortexref,allmarkerfilt,1)
write.table(InpsMarkermap,'Inps_reso0.21_MakerMap_t1.txt',sep='\t',row.names=F,quote=F)

load('hipallmarkerfilt.rda')
hipsubclass <- allmarkerfilt
hipsubclass <- hipsubclass[hipsubclass$avg_log2FC>0,]
allenhipsub <- as.vector(unique(hipsubclass$cluster))
allenhipsub <- allenhipsub[-c(17,18,22,24)]
#allen subclass
load('hipallmarkerfiltharmony-inps.reso0.21.rda')
hipmyclus <- as.vector(unique(allmarkerfilt$cluster))
InpsMarkermap <- MAKERMAPSUB(hipmyclus,refsub=allenhipsub,refs=hipsubclass,allmarkerfilt,1)
write.table(InpsMarkermap,'Inps_reso0.21_MakerMap_allensubclass_t1.txt',sep='\t',row.names=F,quote=F)

#替换亚型
table(inps$RNA_snn_res.0.21,inps$Type)
CelltypeInps <- data.frame(Cell=rownames(inps@meta.data),
                      Celltype=inps$RNA_snn_res.0.21)
CelltypeInps$Celltype <- gsub('^0$','CA1 Pyramidal neurons',CelltypeInps$Celltype)
CelltypeInps$Celltype <- gsub('^1$','DG Pyramidal neurons',CelltypeInps$Celltype)
CelltypeInps$Celltype <- gsub('^2$','Intermediate progenitors',CelltypeInps$Celltype)
rownames(CelltypeInps) <- CelltypeInps$Cell
inps$Celltype <- CelltypeInps[colnames(inps),'Celltype']
save(inps,file='hipreclusterinps.rda')

# 新建Subtype替换inps亚型
load('Inpsub/hipreclusterinps.rda')
hipmeta[rownames(CelltypeInps),'Celltype'] <- CelltypeInps$Celltype
hipmeta$Cell <- rownames(hipmeta)
table(hipmeta$Celltype)
#结合细胞类型和样本类型
hipmeta$CellSamtype <- paste0(hipmeta$Celltype,' ',hipmeta$Type)
hipfinal$Celltype <- hipmeta[colnames(hipfinal),'Celltype']
hipfinal$CellSamtype <- hipmeta[colnames(hipfinal),'CellSamtype']
table(hipfinal$Celltype)
table(hipfinal$CellSamtype)
save(hipmeta,file='hipmeta.rda')
save(hipfinal,file='hipfinal.rda')

#最终的markermap
load('hipallmarkerfilt.rda')
hipsubclass <- allmarkerfilt
hipsubclass <- hipsubclass[hipsubclass$avg_log2FC>2,]
allenhipsub <- as.vector(unique(hipsubclass$cluster))
allenhipsub <- allenhipsub[-c(17,18,22,24)]
load('hipfinal_allmarkerfilt.rda')
hipmyclus <- as.vector(unique(allmarkerfilt$cluster))
Markermap <- MAKERMAP(hipmyclus,refsub=allenhipsub,refs=hipsubclass,allmarkerfilt,1)
write.table(Markermap,'Celltypemarkermap_allenclass_t1.txt',sep='\t',row.names=F,quote=F)

#合并细胞类型
multiType <- list(c('CA1 Pyramidal neurons','CA2,3 Pyramidal neurons'),
                c('GABAergic neurons','Interneurons'),
                c('Cycling glial cells','Intermediate progenitors'))
NewType <- data.frame(Cell=colnames(hipfinal),NewType=hipfinal$Celltype,Type=hipfinal$Type,Sample=hipfinal$Sample)
NewType$NewType <- ifelse(NewType$NewType %in% multiType[[1]],'CA Pyramidal neurons',NewType$NewType)
NewType$NewType <- ifelse(NewType$NewType %in% multiType[[2]],'Inhibitory neurons',NewType$NewType)
NewType$NewType <- ifelse(NewType$NewType %in% multiType[[3]],'Cycling cells',NewType$NewType)
rownames(NewType) <- NewType$Cell
unique(NewType$NewType)
hipmeta$NewType <- NewType[rownames(hipmeta),'NewType']
hipfinal$NewType <- NewType[colnames(hipfinal),'NewType']
save(hipmeta,file='hipmeta.rda')
save(hipfinal,file='hipfinal.rda')
hipfinalexp <- hipfinal@assays$RNA$data
save(hipfinalexp,file='hipfinalexp.rda')

#指定颜色
colourCount <- length(unique(hipmeta$Celltype))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
mycol <- getPalette(colourCount)
mycol <- c(mycol,'#4DBBD5','Gold')
names(mycol) <- c(unique(hipmeta$Celltype),'wt','mut')
mycol['Intermediate progenitors'] <- 'MediumPurple1'
mycol['DG Pyramidal neurons'] <- 'Cyan'
mycol['Astrocytes'] <- 'Blue'
mycol['CA1 Pyramidal neurons'] <- 'GreenYellow'
mycol['Pericytes'] <- '#FFF8DC'
mycol['Inhibitory neurons'] <- '#ffbd15'
mycol['Microglia'] <- '#FACBE5'
save(mycol,file='Hipmycol.rda')

#重新画降维图
#UMAP
p2 <- DimPlot(hipfinal,cols=mycol,reduction = 'harmony_umap',group.by = c('Type','Celltype'),combine = TRUE, label.size = 0.1)
ggsave(p2,filename='Hipfinalumap.pdf',width=13,height=5)

#TSNE
p2 <- DimPlot(hipfinal,cols=mycol,reduction = 'harmony_tsne',group.by = c('Type','Celltype'),combine = TRUE, label.size = 0.1)
ggsave(p2,filename='Hipfinaltsne.pdf',width=13,height=5)

#分野生型和突变体画TSNE
wtcell <- subset(hipmeta,Type=='wt')[,'Cell']
hipwt <- hipfinal[,wtcell]
mutcell <- subset(hipmeta,Type=='mut')[,'Cell']
hipmut <- hipfinal[,mutcell]
p1 <- DimPlot(hipwt,cols=mycol,reduction = 'harmony_tsne',group.by = c('Celltype'),combine = TRUE, label.size = 0.1)+ggtitle('Zbtb18+/+\n(n=9875)')
p2 <- DimPlot(hipmut,cols=mycol,reduction = 'harmony_tsne',group.by = c('Celltype'),combine = TRUE, label.size = 0.1)+ggtitle('Zbtb18-/-\n(n=12081)')
p <- p1+p2+patchwork::plot_layout(ncol = 2, nrow = 1,guides='collect')
ggsave(p,filename='Hipfinaltsne2.pdf',width=11,height=5)


#重新findmarker
Idents(hipfinal) <- hipfinal$Celltype
allmarkers <- FindAllMarkers(hipfinal, only.pos = FALSE, min.pct = 0.2, 
                        logfc.threshold = 0.1,assay='RNA',slot='data')
save(allmarkers,file='hipfinal_allmarker.rda')
allmarkerfilt <- allmarkers[allmarkers$p_val_adj<0.05,]
save(allmarkerfilt,file='hipfinal_allmarkerfilt.rda')
write.csv(allmarkerfilt,file='hipfinal_allmarkerfilt.csv',row.names=F,quote=T)

#complex heatmap
#多scale一些marker，要不画图不好看
#处理高变基因
library(ComplexHeatmap)
hipfinal <- ScaleData(hipfinal, features = rownames(hipfinal@assays$RNA$data),vars.to.regress = c("percent.mt", "percent.Rpsl", "percent.Hb"))
hipfinalnewscale <- hipfinal@assays$RNA$scale.data#这个hipfinal就不用保存了
save(hipfinalnewscale,file='hipfinalExpForplot.rda')

load('hipfinal_allmarkerfilt.rda')
Idents(hipfinal) <- hipfinal$Celltype
celltypes <- unique(hipfinal$Celltype)[order(unique(hipfinal$Celltype))]
#按照细胞数量多少排序
celltypes1 <- names(table(hipfinal$Celltype)[order(table(hipfinal$Celltype),decreasing=T)])
markplot3 <- lapply(celltypes1,function(x){
  tmp <- allmarkerfilt[allmarkerfilt$cluster==x,] 
  tmp <- tmp[order(tmp$avg_log2FC,decreasing=T),]
  tmp <- tmp$gene %>% .[1:15]
  return(tmp)
}) %>% unlist(.)

annogene <- lapply(celltypes1,function(x){
  tmp <- allmarkerfilt[allmarkerfilt$cluster==x,] 
  tmp <- tmp[order(tmp$avg_log2FC,decreasing=T),]
  tmp <- tmp$gene %>% .[1:4]
  return(tmp)
}) %>% unlist(.)
group.use <- c()
for (re in celltypes1){
  tmp <- hipmeta[hipmeta$Celltype==re,'Cell']
  group.use <- c(group.use,tmp)
}

heatbox <- hipfinalnewscale[markplot3,group.use] %>% data.frame(.,check.names=F)
ann <- data.frame(hipmeta[colnames(heatbox),'Celltype'])
colnames(ann) <- 'Celltype'
#颜色
load('Hipmycol.rda')
colours <- list('Celltype' = mycol)
colAnn  <- ComplexHeatmap::HeatmapAnnotation(df = ann,
  which = 'column',
  col = colours)
theatbox <-  t(heatbox) %>% pheatmap:::scale_rows(.)
theatbox <- t(theatbox)
colors2 <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
colors <- colorRampPalette(colors =colors2 )(100)
n_col = length(colors)
lim = max(abs(theatbox), na.rm = TRUE)
library(colorRamp2)
colors2 <- colorRamp2(seq(-lim, lim, length = n_col), colors)
colsplit <- factor(ann$Celltype,levels=celltypes1,labels=c(1:length(unique(hipmeta$Celltype))))
p <- ComplexHeatmap::Heatmap(theatbox,show_column_names=FALSE,cluster_columns = FALSE,cluster_rows = FALSE,
                        show_row_dend = FALSE,show_column_dend = FALSE,#不展示聚类树
                        show_row_names = FALSE,#不展示横坐标
                        col = colors2, 
                        border=TRUE,border_gp = gpar(col = "white", lty = 1,lwd=0.5),
                        top_annotation=colAnn,
                        column_split=colsplit,
                        column_title = NULL,
                        )
pos <-  which(rownames(heatbox) %in% annogene)
p <- p + rowAnnotation(link = anno_mark(at = pos, 
    labels = rownames(heatbox)[pos], labels_gp = gpar(fontsize = 25)))#labels_rot可以调整角度
pdf('allmarkerComplexHeatTop15.pdf',width=12,height=20)
p
dev.off()

#画marker基因的气泡图
genes_to_plot <- c('Mlc1','Slc1a3','Satb2','Dlg2','Slit3','Grin2a','Pif1','Depdc1a','Prox1os','Prox1','Zfp366','Rassf9','Gad1','Gad2',
                    'Eomes','Pax6','Mrc1','Vav1','Sox10','Gpr17','Slc38a11','Atp13a5','Slc6a13','Col1a2')
hipfinal$Celltype <- factor(hipfinal$Celltype)
Idents(hipfinal) <- hipfinal$Celltype
p <- DotPlot(hipfinal, features = genes_to_plot) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=15))
ggsave(p,filename='allmarkerdotplot.pdf',width=11,height=4)

#所有细胞类型的比例########
load('hipmeta.rda')
mutall <- nrow(hipmeta[hipmeta$Type=='mut',])
wtalll <- nrow(hipmeta[hipmeta$Type=='wt',])
number <- table(hipmeta$Type,hipmeta$Celltype)
number <- data.frame(number,check.names=F)
number$Percent <- ifelse(number$Var1=='mut',number$Freq/mutall*100,number$Freq/wtalll*100)
number$Percent <- round(number$Percent,2)
colnames(number)[1:3] <- c('Type','Celltype','Number')
write.table(number,file='Cellnumber-Zbtb18Hip.txt',sep='\t',quote=F,row.names=F)

cellnumber <- table(hipmeta$Celltype)
cellnumber <- data.frame(Celltype=names(cellnumber),cellnumber=as.vector(cellnumber))
write.table(cellnumber,file='Cellnumber2-Zbtb18Hip.txt',sep='\t',quote=F,row.names=F)

#堆叠图，总和是wt或者mut
load('Hipmycol.rda')
number$Type <- factor(number$Type,levels=c('wt','mut'))
p <- ggplot( number, aes( x = Type, weight = Percent, fill = Celltype))+
 geom_bar( position = "stack") + 
 scale_fill_manual(values = mycol)+ #设置填充的颜色
 theme_bw()+ #背景变为白色
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size = 15), 
        panel.border = element_blank(),axis.line = element_line(colour = "black",linewidth=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.position = "right",
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),#不显示网格线
        panel.grid.minor = element_blank())+#显示显著性
  ylab("Percent of celltype")+xlab("")
ggsave(p,filename='celltypePercentAll.pdf',width=6,height=10)

#堆叠图，总和是每个细胞类型
data <- lapply(unique(hipmeta$Subtype),function(i){
  metatmp <- subset(hipmeta,Subtype==i)
  per <- table(metatmp$Subtype,metatmp$Type) %>% data.frame(.,check.names=F)
  per$percent <- per$Freq/sum(per$Freq)
  return(per)
})
data <- do.call(rbind,data)
colnames(data) <- c('Celltype','Type','count','percent')
data$Celltype <- factor(data$Celltype,levels=unique(data$Celltype)[c(6,5,11,7,1,8,13,2:4,9,10,12,14,15)])
p <- ggplot( data, aes( x = Celltype, weight = percent, fill = Type))+
 geom_bar( position = "stack") + 
 scale_fill_manual(values = c('Gold','#4DBBD5'))+ #设置填充的颜色
 theme_bw()+ #背景变为白色
  theme(axis.text=element_text(size=15),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1),
        axis.title=element_text(size = 15), 
        panel.border = element_blank(),axis.line = element_line(colour = "black",linewidth=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.background = element_blank(),#不显示网格线
        panel.grid.minor = element_blank())+#显示显著性
  ylab("Percent of celltype")+xlab("")
ggsave(p,filename='celltypePercent.pdf',width=10,height=5)

#差异表达
#不计算subtype了，有的subtype，wu和mut数量差异太大，无法计算差异基因
minpct <- 0.1
logfc <- 0.1
Idents(hipfinal) <- hipfinal$Celltype
DEGs <- data.frame()
for (type in unique(hipfinal$Celltype)){
  print(type)
  DEGtmp <- FindMarkers(hipfinal, only.pos = FALSE, min.pct = minpct,
                        test.use = "wilcox",logfc.threshold = logfc,assay='RNA',slot='data',
                        ident.1 = "mut", group.by = 'Type', subset.ident = type)#ident.1=case
  DEGtmp$Celltype <- type
  DEGtmp$gene <- rownames(DEGtmp)
  DEGs <- rbind(DEGs,DEGtmp)
}
save(DEGs,file='DEG_hipfinal.rda')
DEGfilt <- DEGs[DEGs$p_val_adj<0.05,]
save(DEGfilt,file='DEGfilt_hipfinal.rda')
write.table(DEGfilt,file='DEGfilt_hipfinal.txt',row.names=F,quote=F,sep='\t') 

DEGsta <- table(DEGfilt$Celltype)
DEGsta <- data.frame(Celltype=names(DEGsta),number=as.vector(DEGsta))
write.table(DEGsta,file='DEGsta_hipfinal.txt',row.names=F,quote=F,sep='\t') 

#目标基因
load('DEG_hipfinal.rda')
load('hipfinal.rda')
test <- subset(DEGs,Celltype %in% c('Cycling glial cells','Intermediate progenitors') & gene %in% c('Ctnnb1','Wnt1','Wnt2','Wnt2b','Wnt3','Wnt3a','Wnt4','Wnt5a','Wnt5b','Wnt6','Wnt7a','Wnt7b','Wnt8a','Wnt8b','Wnt9a','Wnt9b','Wnt10b','Wnt10a','Wnt11','Wnt16'))
write.table(test,file='Ctnnb1+Wnt.txt',sep='\t',row.names=F,quote=F)

#1216目标基因
load('DEG_hipfinal.rda')
genes <- read.table('gene1216.txt',sep='\t',check.names=F,header=T) %>% .$gene
test <- subset(DEGs,gene %in% genes)
write.table(test,file='HipgenesDEinfo_1216.txt',sep='\t',row.names=F,quote=F)


#DG椎体神经元差异基因火山图
library(ggrepel)
load('DEG_hipfinal.rda')
dat <- DEGs[DEGs$Celltype=='DG Pyramidal neurons',]
log2fc <- 0.1
dat$Significant[dat$avg_log2FC > log2fc & dat$p_val_adj < 0.05] <- 'Up'
dat$Significant[dat$avg_log2FC < (-log2fc) & dat$p_val_adj < 0.05] <- 'Down'
dat$Significant[!(abs(dat$avg_log2FC) > log2fc  & dat$p_val_adj < 0.05)] <- 'No'
mycol <-c('Blue','#3D3D3D','Coral')
names(mycol) <- c('Down','No','Up')
dat$label <- ifelse(dat$gene %in% c('Pax6','Grik1','Grin2a','Slit3'),dat$gene,'')

# dat$label <- unlist(sapply(dat$gene,function(x){
#                     ifelse(x %in% astroDEmarker,x,"")}))
p <-ggplot(dat,aes(x=avg_log2FC,y=-log10(p_val_adj))) + 
        geom_point(size=4,aes(color=Significant,fill=Significant),shape = 21) + 
        scale_fill_manual(values=mycol,name = 'Significant') +
        scale_color_manual(values=mycol,name = 'Significant')+
        geom_hline(yintercept=-log10(0.05),linetype=3,lwd = 1) +geom_vline(xintercept=c(-log2fc,log2fc),linetype=3,lwd = 1) +
        theme_classic()+
        theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
              legend.text = element_text(size = 15),text = element_text(size = 15),
              legend.position='left')+
        geom_text_repel(data = dat, size=6,aes(x = avg_log2FC, #如果不想要框框，就用geom_text_repel,否则就用label_repel
                              y = -log10(p_val_adj), 
                              label = label), color="black", 
                          box.padding=unit(0.35, "lines"), point.padding=unit(1.6, "lines"),
                          segment.colour = "black",max.overlaps=20000)#图像尺寸大一些，否则就会画不出来框
ggsave(p,filename='DG volcano.pdf',width=6,height=4)

#其他细胞类型，不标注基因
for(ct in unique(DEGs$Celltype)[-3]){
  dat <- DEGs[DEGs$Celltype==ct,]
  log2fc <- 0.1
  dat$Significant[dat$avg_log2FC > log2fc & dat$p_val_adj < 0.05] <- 'Up'
  dat$Significant[dat$avg_log2FC < (-log2fc) & dat$p_val_adj < 0.05] <- 'Down'
  dat$Significant[!(abs(dat$avg_log2FC) > log2fc  & dat$p_val_adj < 0.05)] <- 'No'
  mycol <-c('Blue','#3D3D3D','Coral')
  names(mycol) <- c('Down','No','Up')
  # dat$label <- unlist(sapply(dat$gene,function(x){
  #                     ifelse(x %in% astroDEmarker,x,"")}))
  p <-ggplot(dat,aes(x=avg_log2FC,y=-log10(p_val_adj))) + 
        geom_point(size=4,aes(color=Significant,fill=Significant),shape = 21) + 
        scale_fill_manual(values=mycol,name = 'Significant') +
        scale_color_manual(values=mycol,name = 'Significant')+
          geom_hline(yintercept=-log10(0.05),linetype=3,lwd = 1) +geom_vline(xintercept=c(-log2fc,log2fc),linetype=3,lwd = 1) +
          theme_classic()+
          theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
                legend.text = element_text(size = 15),text = element_text(size = 15),
                legend.position='left')
  ggsave(p,filename=paste0(ct,' volcano.pdf'),width=6,height=4)
}

#细胞特异性的差异基因小提琴图
library(stringi)
library(stringr)
library(tidyr)
library(dplyr)
library(Seurat)
library(ggplot2)
load('hipfinal.rda')
load('hipmeta.rda')
load('DEGfilt_hipfinal.rda')
load('hipfinal_allmarkerfilt.rda')
mycelltype <- unique(DEGfilt$Celltype)[c(2:10)]
mycelltype <- mycelltype[order(mycelltype)]
mycell <- subset(hipmeta,Celltype %in% mycelltype)[,'Cell']
plotgene <- c()
for(ct in mycelltype){
  tmp <- subset(DEGfilt,Celltype==ct)
  tmp <- subset(tmp,avg_log2FC>1 | avg_log2FC <(-1))
  mygene <- tmp[,'gene']
  mygene <- intersect(mygene,subset(allmarkerfilt,cluster==ct & avg_log2FC > 2)[,'gene'])
  mygene <- mygene[1:3]
  plotgene <- c(plotgene,mygene)
}
plotgene <- unique(plotgene)
plotgene <- na.omit(plotgene)

hipfinal$Type <- factor(hipfinal$Type,levels=c('wt','mut'))
hipfinal$Celltype <- gsub('CA1 Pyramidal neurons','CA1\nPyramidal neurons',hipfinal$Celltype)
hipfinal$Celltype <- gsub('CA2,3 Pyramidal neurons','CA2,3\nPyramidal neurons',hipfinal$Celltype)
hipfinal$Celltype <- gsub('Intermediate progenitors','Intermediate\nprogenitors',hipfinal$Celltype)
hipfinal$Celltype <- gsub('Cycling glial cells','Cycling\nglial cells',hipfinal$Celltype)

p <- VlnPlot(object = hipfinal[,mycell],features = plotgene,  # 你的基因名
  group.by = "Celltype",  # 分组变量（细胞类型）pt.size = 0,  # 不要散点
  split.by='Type',
  stack = TRUE,  # 堆叠模式
  flip = TRUE,   # 让基因在Y轴，celltype在X轴
  cols = c('#4DBBD5','#E64B35')) +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title.x=element_blank(),
        strip.text.y = element_text(size = 15,vjust=1,hjust=0,face='italic'))
ggsave(p,filename='HipDEGmarkerdot.pdf',width=6,height=8)

#画图#######################################
#marker基因或者感兴趣基因 某个基因 的降维图UMAP TSNE###########
colsedit <- c(brewer.pal(n = 11, name = "Spectral"))
colsedit[11] <- '#e8d4e8'
fix.sc <- scale_colour_gradientn(colours = rev(colsedit))
# features <- c('Prox1','Eomes')
# features <- 'Satb2'
features <- c('Synpr','Nrg2')
Idents(hipfinal) <- hipfinal$Celltype
DefaultAssay(hipfinal) <- 'RNA'
plotdata <- hipfinal
for ( i in features){
    pp <- quantile(plotdata@assays$RNA$data[i,],1)
    plotdata@assays$RNA$data[i,] <- sapply(plotdata@assays$RNA$data[i,],
                                        function(x) {ifelse(x>pp,pp,x)})
    p1 <- FeaturePlot(object = plotdata, reduction = "harmony_tsne",#slot='scale.data',
                  features = i,combine = FALSE,pt.size=0.1)
    p2 <- lapply(p1,function(x)x+fix.sc+theme(legend.position="none",
                                              axis.text=element_text(size=20))+
                                              xlab('')+ylab('')) 
    p3 <- aplot::plot_list(gglist=p2,ncol=1)
    ggsave(p3,filename=paste0('plot/tsne_',i,'.pdf'),width=5,height=5)
}
#画一个colorbar
p1 <- FeaturePlot(object = plotdata, reduction = "harmony_umap",#slot='scale.data',
                features = i,combine = FALSE,pt.size=0.1)
p2 <- lapply(p1,function(x)x+fix.sc+xlab('')+ylab('')) 
p3 <- aplot::plot_list(gglist=p2,ncol=1)
ggsave(p3,filename='colorbar.pdf',width=5,height=5)

#cellmarker画图，每个区域多画一些，挑最好的
load('hipfinal.rda')
load('hipfinal_allmarkerfilt.rda')
markermap <- read.table('Celltypemarkermap_allenclass_t1.txt',sep='\t',header=T,check.names=F)
Idents(hipfinal) <- hipfinal$Celltype
celltype <- 'VLMC'
subs = 'VLMC'
library(viridis)
colsedit <- c(brewer.pal(n = 9, name = "YlOrRd"))[2:8]
fix.sc <- scale_colour_gradientn(colours = colsedit)

features <- subset(markermap,mycluster==celltype & subs == sub)[,'intergene'] %>% str_split_fixed(.,pattern=',',n=Inf) %>% as.vector(.)
features <- subset(allmarkerfilt,cluster == celltype & gene %in% features)
features <- features[order(features$avg_log2FC,decreasing=T),]
selectn <- as.integer(nrow(features) * 0.4)
features <- features[1:selectn,'gene']
DefaultAssay(hipfinal) <- 'RNA'
plotdata <- hipfinal
for ( i in features){
    pp <- quantile(plotdata@assays$RNA$data[i,],1)
    plotdata@assays$RNA$data[i,] <- sapply(plotdata@assays$RNA$data[i,],
                                        function(x) {ifelse(x>pp,pp,x)})
    p1 <- FeaturePlot(object = plotdata, reduction = "harmony_tsne",#slot='scale.data',
                  features = i,combine = FALSE,pt.size=0.1)
    p2 <- lapply(p1,function(x)x+fix.sc+theme(legend.position="none",
                                              axis.text=element_text(size=20))+
                                              xlab('')+ylab('')) 
    p3 <- aplot::plot_list(gglist=p2,ncol=1)
    ggsave(p3,filename=paste0('tsne_',i,'_',celltype,'.pdf'),width=5,height=5)
}

#多个基因全局表达/半全局表达的小提琴图合并散点图,手动加log2fc
#小提琴图+散点图
library(ggpubr)
library(ggrepel)
gene <- c('Zbtb18')
mycell <- hipmeta$Cell
for (i in gene){
  myexp <- hipexp[i,mycell] %>% data.frame(.,check.names=F)
  colnames(myexp)[1] <- 'expression'
  myexp$gene <- i
  myexp$Type <- hipmeta[rownames(myexp),'Type']
  max <- max(myexp$expression)+1          
  log2FC <- aggregate(expression~Type,data=myexp,mean)
  log2FC <- log2(log2FC[log2FC$Type=='mut','expression']/log2FC[log2FC$Type=='wt','expression'])
  p1 <-  ggplot(myexp, aes(x=Type, y=expression,fill=Type)) + 
      geom_violin(width=0.6,trim=TRUE,color="#EED8AE",linewidth=0.6,
                  scale = "width",position=position_dodge(width=0.6)) + #小提琴图
      geom_jitter(size=0.3,colour='#EED8AE',
                  position=position_jitterdodge(jitter.width = 0.35, 
                                                jitter.height = 0, 
                                                dodge.width = 0.6))+#散点
      stat_summary(fun.data = 'mean_se',
                  size = 0.6, color = "Tomato",
                  position=position_dodge(0.6))+
    scale_fill_manual(values = c('#E64B35','#4DBBD5'))+ #设置填充的颜色
    theme_bw()+ #背景变为白色
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size = 15), 
          panel.border = element_blank(),axis.line = element_line(colour = "black",linewidth=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.background = element_blank(),#不显示网格线
          panel.grid.minor = element_blank())+#显示显著性
    ylab("Normalized Expression")+xlab("")+ 
    stat_compare_means(size = 5,label.x=1.5,label = 'p.signif',label.y=3)+
    annotate('text',x=1.5,y=max-1,label=paste0('log2FC = ',round(log2FC,2)))
  ggsave(p1,filename = paste0(i,'_violin.pdf'),height = 3,width=3)
}

#单个基因在显著差异的细胞类型中的表达,手动加pvalue
#pvalue可以用pval或者p.adj
load('DEGfilt_hipfinal.rda')
load('hipfinalexp.rda')
load('hipmeta.rda')
load('hipfinal.rda')
gene <- c('Zbtb18','Grin2a')
gene <- c('Pax6','Satb2')#3
gene <- 'Rbfox3'#4
gene <- 'Klf8'
thiswidth=3
thisheight=3.5
for (i in gene){
  x.pos <- 0
  feainfo <- data.frame()
  dbox <- data.frame()
  hipfinalexp <- hipfinal@assays$RNA$data
  thistype <- unique(DEGfilt[DEGfilt$gene==i,'Celltype'])
  thistype <- thistype[order(thistype)]#多个细胞类型时候，一定按照字母顺序排好序
  thislength <- length(thistype)
  print(paste0(i,thislength))
  if(thislength>0){
    for (celltype in thistype){
      x.pos <- x.pos+1
      thisexp <- hipfinalexp[,hipmeta[hipmeta$Celltype == celltype,'Cell']]
      dbox1 <- data.frame(expression=thisexp[i,])
      dbox1$gene <- i
      dbox1$Celltype <- hipmeta[rownames(dbox1),'Celltype']
      dbox1$class_label <- hipmeta[rownames(dbox1),'class_label']
      dbox1$Type <- hipmeta[rownames(dbox1),'Type']
      # print(min(dbox1$expression))
      #prepare stat.test
      stat.test <- DEGfilt[DEGfilt$gene==i & DEGfilt$Celltype==celltype,]
      stat.test <- cbind(data.frame(`.y.`=rep('expression',1)),stat.test,data.frame(method='Wilcox'))
      stat.test$group1 <- c('mut')
      stat.test$group2 <- c('wt')
      stat.test$p.format <- as.character(format(stat.test$p_val_adj,scientific = TRUE))
      stat.test$p.signif <- ifelse(stat.test$p_val_adj<0.05,
                                    ifelse(stat.test$p_val_adj<0.01,
                                      ifelse(stat.test$p_val_adj<0.001,'***','**'),'*'),'ns')                            
      stat.test$xpos <- c(x.pos+0)
      feainfo <- rbind(feainfo,stat.test)
      dbox <- rbind(dbox,dbox1)
    }
    max <- max(dbox$expression)
    feainfo$`y.position` <- max+1
    dbox$Type <- factor(dbox$Type,levels=c('mut','wt'))
    dbox$Celltype <- gsub('Intermediate progenitors','Intermediate\nprogenitors',dbox$Celltype)
    dbox$Celltype <- gsub('Endothelial cells','Endothelia\ncells',dbox$Celltype)
    dbox$Celltype <- gsub('Pyramidal neurons','Pyramidal\nneurons',dbox$Celltype)
    dbox$Celltype <- gsub('GABAergic neurons','GABAergic\nneurons',dbox$Celltype)
    dbox$Celltype <- gsub('Cycling glial cells','Cycling glial\ncells',dbox$Celltype)
    p <- ggplot(dbox, aes(x=Celltype, y=expression,fill=Type)) + 
    geom_violin(width=0.6,trim=TRUE,color="#EED8AE",linewidth=0.6,
                scale = "width",position=position_dodge(width=0.6)) + #小提琴图
    geom_jitter(size=0.3,colour='#EED8AE',
                position=position_jitterdodge(jitter.width = 0.35, 
                                              jitter.height = 0, 
                                              dodge.width = 0.6))+#散点
    stat_summary(fun.data = 'mean_se',
                size = 0.6, color = "Tomato",
                position=position_dodge(0.6))+
    # scale_fill_manual(values = c("DarkOrange","SpringGreen","RoyalBlue1","Yellow1","OrangeRed1"))+ #设置填充的颜色
    scale_fill_manual(values=c('#E64B35','#4DBBD5'))+
    theme_bw()+ #背景变为白色
    theme(axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10,angle=90,vjust=1,hjust=1),
          axis.title=element_text(size = 10), 
          panel.border = element_blank(),axis.line = element_line(colour = "black",linewidth=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
          legend.position = "top",
          panel.grid.major = element_blank(),
          panel.background = element_blank(),#不显示网格线
          panel.grid.minor = element_blank())+#显示显著性
    ylab(paste0(i," Expression"))+xlab("")+
    annotate(geom='text',size=5,x=feainfo$xpos,
            y=feainfo$y.position,label=feainfo$p.signif,
            color='black') 
    ggsave(p,filename = paste0(i,'_inSigCelltypes.pdf'),width=thiswidth,height=thisheight)
  }
}

#pvalue,pax6
thiswidth=3
thisheight=3
gene <- 'Pax6'
load('DEG_hipfinal.rda')
for (i in gene){
  x.pos <- 0
  feainfo <- data.frame()
  dbox <- data.frame()
  hipfinalexp <- hipfinal@assays$RNA$data
  thistype <- unique(DEGs[DEGs$gene==i,'Celltype'])
  thistype <- thistype[order(thistype)]#多个细胞类型时候，一定按照字母顺序排好序
  thislength <- length(thistype)
  print(paste0(i,thislength))
  if(thislength>0){
    for (celltype in thistype){
      x.pos <- x.pos+1
      thisexp <- hipfinalexp[,hipmeta[hipmeta$Celltype == celltype,'Cell']]
      dbox1 <- data.frame(expression=thisexp[i,])
      dbox1$gene <- i
      dbox1$Celltype <- hipmeta[rownames(dbox1),'Celltype']
      dbox1$class_label <- hipmeta[rownames(dbox1),'class_label']
      dbox1$Type <- hipmeta[rownames(dbox1),'Type']
      # print(min(dbox1$expression))
      #prepare stat.test
      stat.test <- DEGs[DEGs$gene==i & DEGs$Celltype==celltype,]
      stat.test <- cbind(data.frame(`.y.`=rep('expression',1)),stat.test,data.frame(method='Wilcox'))
      stat.test$group1 <- c('mut')
      stat.test$group2 <- c('wt')
      stat.test$p.format <- as.character(format(stat.test$p_val,scientific = TRUE))
      stat.test$p.signif <- ifelse(stat.test$p_val<0.05,
                                    ifelse(stat.test$p_val<0.01,
                                      ifelse(stat.test$p_val<0.001,'***','**'),'*'),'ns')                            
      stat.test$xpos <- c(x.pos+0)
      feainfo <- rbind(feainfo,stat.test)
      dbox <- rbind(dbox,dbox1)
    }
    max <- max(dbox$expression)
    feainfo$`y.position` <- max+1
    dbox$Type <- factor(dbox$Type,levels=c('mut','wt'))
    dbox$Celltype <- gsub('Intermediate progenitors','Intermediate\nprogenitors',dbox$Celltype)
    dbox$Celltype <- gsub('Endothelial cells','Endothelia\ncells',dbox$Celltype)
    dbox$Celltype <- gsub('Pyramidal neurons','Pyramidal\nneurons',dbox$Celltype)
    dbox$Celltype <- gsub('GABAergic neurons','GABAergic\nneurons',dbox$Celltype)
    dbox$Celltype <- gsub('Cycling glial cells','Cycling glial\ncells',dbox$Celltype)
    p <- ggplot(dbox, aes(x=Celltype, y=expression,fill=Type)) + 
    geom_violin(width=0.6,trim=TRUE,color="#EED8AE",linewidth=0.6,
                scale = "width",position=position_dodge(width=0.6)) + #小提琴图
    geom_jitter(size=0.3,colour='#EED8AE',
                position=position_jitterdodge(jitter.width = 0.35, 
                                              jitter.height = 0, 
                                              dodge.width = 0.6))+#散点
    stat_summary(fun.data = 'mean_se',
                size = 0.6, color = "Tomato",
                position=position_dodge(0.6))+
    # scale_fill_manual(values = c("DarkOrange","SpringGreen","RoyalBlue1","Yellow1","OrangeRed1"))+ #设置填充的颜色
    scale_fill_manual(values=c('Gold','#4DBBD5'))+
    theme_bw()+ #背景变为白色
    theme(axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10,angle=90,vjust=1,hjust=1),
          axis.title=element_text(size = 10), 
          panel.border = element_blank(),axis.line = element_line(colour = "black",linewidth=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.background = element_blank(),#不显示网格线
          panel.grid.minor = element_blank())+#显示显著性
    ylab(paste0(i," Expression"))+xlab("")+
    annotate(geom='text',size=5,x=feainfo$xpos,
            y=feainfo$y.position,label=feainfo$p.signif,
            color='black') 
    ggsave(p,filename = paste0(i,'_inCelltypes-pvalue.pdf'),width=thiswidth,height=thisheight)
  }
}

#单个基因，不在DEGs里（findallmarker检验），在多种细胞的小提琴图,用pvalue（秩和检验）
#单个基因在class label中的表达,Het与WT比较
i <- 'Gli2'
x.pos <- 0
feainfo <- data.frame()
dbox <- data.frame()
hipexp <- hipfinal@assays$RNA$data
for (celltype in c('DG Pyramidal neurons','CA1 Pyramidal neurons','CA2,3 Pyramidal neurons')){
  x.pos <- x.pos+1
  dbox1 <- data.frame(expression=hipexp[i,])
  dbox1$gene <- i
  dbox1$Celltype <- hipmeta[rownames(dbox1),'Celltype']
  dbox1$Type <- hipmeta[rownames(dbox1),'Type']
  dbox1 <- dbox1[dbox1$Celltype==celltype,]
  print(min(dbox1$expression))
  #prepare stat.test
  test <- wilcox.test(expression~Type,dbox1)
  p_val <- test$p.value
  stat.test <- data.frame(`.y.`=rep('express',1),method='Wilcox')
  stat.test$group1 <- c('mut')
  stat.test$group2 <- c('wt')
  stat.test$p.format <- as.character(format(p_val,scientific = TRUE))
  stat.test$p.signif <- ifelse(p_val<0.05,
                                ifelse(p_val<0.01,
                                  ifelse(p_val<0.001,'***','**'),'*'),'ns')                       
  max <- max(dbox1$expression)
  stat.test <- stat.test %>%
                mutate(y.position = max+1)
  stat.test$xpos <- c(x.pos+0)
  feainfo <- rbind(feainfo,stat.test)
  dbox <- rbind(dbox,dbox1)
}
dbox$Celltype <- gsub('Pyramidal neurons','Pyramidal\nneurons',dbox$Celltype)
p <- ggplot(dbox, aes(x=Celltype, y=expression,fill=Type)) + 
      geom_violin(width=0.6,trim=TRUE,color="#EED8AE",linewidth=0.6,
                  scale = "width",position=position_dodge(width=0.6)) + #小提琴图
      geom_jitter(size=0.2,colour='#EED8AE',
                  position=position_jitterdodge(jitter.width = 0.35, 
                                                jitter.height = 0, 
                                                dodge.width = 0.6))+#散点
      stat_summary(fun.data = 'mean_se',
                  size = 0.6, color = "Tomato",
                  position=position_dodge(0.6))+
      scale_fill_manual(values=c('Gold','#4DBBD5'))+
      theme_bw()+ #背景变为白色
      theme(axis.text.x=element_text(size=10,angle=90,vjust=1,hjust=1),
            axis.text.y=element_text(size=10),
            axis.title=element_text(size = 10), 
            panel.border = element_blank(),axis.line = element_line(colour = "black",linewidth=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.background = element_blank(),#不显示网格线
            panel.grid.minor = element_blank())+#显示显著性
      ylab(paste0(i," Expression"))+xlab("")+
      annotate(geom='text',size=5,x=feainfo$xpos,
              y=feainfo$y.position,label=feainfo$p.signif,
              color='black') 
ggsave(p,filename = paste0(i,'inCelltype-wilcox.pdf'),width=3,height=3)

#单个基因在某种细胞中的小提琴图加散点图，首先要算差异表达基因哟~手动加pvalue
load('DEG_hipfinal.rda')
load('hipfinalexp.rda')
load('hipmeta.rda')
fea <- c('Zbtb18','Satb2')
fea <- c('Gria1','Kdm5d')
fea <- c('Pde7b','Grik1')
fea <- 'Neurod1'
for (i in fea){
    gene <- i
    Celltype <- 'DG Pyramidal neurons'
    dbox1 <- hipfinalexp[gene,hipmeta[hipmeta$Celltype == Celltype,'Cell']] %>% data.frame(.,check.names=F)
    colnames(dbox1)[1] <- 'expression'
    dbox1$gene <- gene
    dbox1$Type <- hipmeta[rownames(dbox1),'Type']
    print(min(dbox1$expression))
    # dbox1 <- dbox1[dbox1$expression>0,]
    #prepare stat.test
    stat.test <- DEGs[DEGs$gene==i & DEGs$Celltype==Celltype,]
    stat.test <- cbind(data.frame(`.y.`='expression'),stat.test,data.frame(method='Wilcox'))
    stat.test$group1 <- c('mut')
    stat.test$group2 <- c('wt')
    stat.test$p.signif <- ifelse(stat.test$p_val_adj<0.05,
                                ifelse(stat.test$p_val_adj<0.01,
                                    ifelse(stat.test$p_val_adj<0.001,'***','**'),'*'),'ns')
    max <- max(dbox1$expression)
    stat.test$y.position <- max+1
    p <- ggplot(dbox1, aes(x=Type, y=expression,fill=Type)) + 
      geom_violin(width=0.6,trim=TRUE,color="#EED8AE",linewidth=0.6,
                scale = "width",position=position_dodge(width=0.6)) + #小提琴图
      geom_jitter(size=0.6,colour='#EED8AE',
                position=position_jitterdodge(jitter.width = 0.35, 
                                              jitter.height = 0, 
                                              dodge.width = 0.6))+#散点
      stat_summary(fun.data = 'mean_se',
                 size = 0.6, color = "Tomato",
                 position=position_dodge(0.6))+
      # scale_fill_manual(values = c("DarkOrange","SpringGreen","RoyalBlue1","Yellow1","OrangeRed1"))+ #设置填充的颜色
      scale_fill_manual(values=c('#E64B35','#4DBBD5'))+
      theme_bw()+ #背景变为白色
      theme(axis.text=element_text(size=15),
            axis.title=element_text(size = 15), 
            plot.title=element_text(size=15),
            panel.border = element_blank(),axis.line = element_line(colour = "black",linewidth=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.background = element_blank(),#不显示网格线
            panel.grid.minor = element_blank())+#显示显著性
      ylab(paste0(i," Expression"))+xlab("")+ggtitle(paste0(Celltype," ",gene))+
      annotate(geom='text',x=1.5,size=5,
              y=stat.test$y.position,label=stat.test$p.signif,
              color='black')
    ggsave(p,filename=paste0(gene,'_',Celltype,'_vio.pdf'),width = 2.5,height = 2.5)
}

#单个基因整体表达
fea <- 'Mki67'
fea <- 'Prox1'
hipexp <- hipfinal@assays$RNA$data
hipmeta <- hipfinal@meta.data
hipmeta$Cell <- rownames(hipmeta)
for (i in fea){
  gene <- i
  dbox1 <- hipexp[gene,] %>% data.frame(.,check.names=F)
  colnames(dbox1)[1] <- 'expression'
  dbox1$gene <- gene
  dbox1$Type <- hipmeta[rownames(dbox1),'Type']
  print(min(dbox1$expression))
  dbox1$Type <- factor(dbox1$Type,levels=c('wt','mut'))
  p <- ggplot(dbox1, aes(x=Type, y=expression,fill=Type)) + 
  geom_violin(width=0.6,trim=TRUE,color="#EED8AE",linewidth=0.6,scale = "width",position=position_dodge(width=0.6)) + #小提琴图
  geom_jitter(size=0.6,colour='#EED8AE',position=position_jitterdodge(jitter.width = 0.35, jitter.height = 0, dodge.width = 0.6))+#散点
  stat_summary(fun.data = 'mean_se',size = 0.6, color = "Tomato",position=position_dodge(0.6))+
  scale_fill_manual(values=c('#4DBBD5','#E64B35'))+
  theme_bw()+ #背景变为白色
  theme(axis.text=element_text(size=15),axis.title=element_text(size = 15), panel.border = element_blank(),
  axis.line = element_line(colour = "black",linewidth=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
  legend.position = "none",panel.grid.major = element_blank(),panel.background = element_blank(),#不显示网格线
  panel.grid.minor = element_blank())+#显示显著性
  ylab(paste0(i," Expression"))+xlab("")+ggtitle(gene)+
  stat_compare_means(size = 5,label.x=1.4,label = 'p.format',label.y=subset(dbox1,gene==i)[,'expression'] %>% max(.)+0.5)
  ggsave(p,filename=paste0('Hip_',gene,'_vio.pdf'),width = 2.5,height = 3)
}

#共定位分析######
#某个基因大于0就认为表达，本来想取中位数之类的，发现都是0
load('hipfinal.rda')
hipexp <- hipfinal@assays$RNA$data
hipmeta <- hipfinal@meta.data
wtnum <- nrow(subset(hipmeta,Type =='wt'))
mutnum <- nrow(subset(hipmeta,Type =='mut'))
hipmeta$Cell <- rownames(hipmeta)
cutoff <- 0
#pax6,ki67##########
pax6cutoff <- colnames(hipexp[,hipexp['Pax6',]>cutoff])
ki67cutoff <- colnames(hipexp[,hipexp['Mki67',]>cutoff])
pax6ki67 <- intersect(pax6cutoff,ki67cutoff)
dbox1 <- data.frame(Cell=pax6ki67,Type=hipmeta[pax6ki67,'Type'],pax6=hipexp['Pax6',pax6ki67],mki67=hipexp['Mki67',pax6ki67])

#条形图
data <- data.frame(table(dbox1$Type))
colnames(data) <- c('Type','number')
data$Type <- factor(data$Type,levels=c('wt','mut'))
data$Percent <- c(data[1,2]/mutnum*100,data[2,2]/wtnum*100)
data$Percent <- round(data$Percent,2)
p <- ggplot(data = data,aes(x = Type, y = Percent),width=1)+
  geom_bar(aes(fill = Type),color = 'grey21',stat = 'identity', position = 'dodge')+
  scale_fill_manual(values = c('#4DBBD5','#E64B35'))+theme_bw()+
  theme(panel.grid =element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.text=element_text(size=15),axis.title=element_text(size = 15),
  legend.text=element_text(size=15))+
  geom_text(aes(label=Percent), color="white", 
          size=4,position=position_dodge(0.9),vjust=1.5)
# scale_y_continuous(breaks = c(0,0.5,1))
ggsave(p,filename='Pax6+Ki67+cellnumber.pdf',width = 3,height = 3)

#双阳性细胞的表达量
dbox1$Type <- factor(dbox1$Type,levels=c('wt','mut'))
p1 <- ggplot(dbox1, aes(x=Type, y=pax6,fill=Type)) + 
geom_violin(width=0.6,trim=TRUE,color="#EED8AE",linewidth=0.6,scale = "width",position=position_dodge(width=0.6)) + #小提琴图
geom_jitter(size=0.6,colour='#EED8AE',position=position_jitterdodge(jitter.width = 0.35, jitter.height = 0, dodge.width = 0.6))+#散点
stat_summary(fun.data = 'mean_se',size = 0.6, color = "Tomato",position=position_dodge(0.6))+
scale_fill_manual(values=c('#4DBBD5','#E64B35'))+
theme_bw()+ #背景变为白色
theme(axis.text=element_text(size=15),axis.title=element_text(size = 15), panel.border = element_blank(),
axis.line = element_line(colour = "black",linewidth=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
legend.position = "none",panel.grid.major = element_blank(),panel.background = element_blank(),#不显示网格线
panel.grid.minor = element_blank())+#显示显著性
ylab('Expression of Pax6\nin Pax6+/Ki67+ Cells')+xlab("")+
stat_compare_means(size = 6,label.x=1.4,label = 'p.signif',label.y=4)

p2 <- ggplot(dbox1, aes(x=Type, y=mki67,fill=Type)) + 
geom_violin(width=0.6,trim=TRUE,color="#EED8AE",linewidth=0.6,scale = "width",position=position_dodge(width=0.6)) + #小提琴图
geom_jitter(size=0.6,colour='#EED8AE',position=position_jitterdodge(jitter.width = 0.35, jitter.height = 0, dodge.width = 0.6))+#散点
stat_summary(fun.data = 'mean_se',size = 0.6, color = "Tomato",position=position_dodge(0.6))+
scale_fill_manual(values=c('#4DBBD5','#E64B35'))+
theme_bw()+ #背景变为白色
theme(axis.text=element_text(size=15),axis.title=element_text(size = 15), panel.border = element_blank(),
axis.line = element_line(colour = "black",linewidth=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
legend.position = "none",panel.grid.major = element_blank(),panel.background = element_blank(),#不显示网格线
panel.grid.minor = element_blank())+#显示显著性
ylab('Expression of Ki67\nin Pax6+/Ki67+ Cells')+xlab("")+
stat_compare_means(size = 6,label.x=1.4,label = 'p.signif',label.y=4)
p <- p1+p2+patchwork::plot_layout(ncol=2,nrow = 1,guides='collect')
ggsave(p,filename='Pax6_Ki67.pdf',width = 5,height = 3)

#Eomes,ki67##########
Eomescutoff <- colnames(hipexp[,hipexp['Eomes',]>cutoff])
ki67cutoff <- colnames(hipexp[,hipexp['Mki67',]>cutoff])
Eomeski67 <- intersect(Eomescutoff,ki67cutoff)
dbox1 <- data.frame(Cell=Eomeski67,Type=hipmeta[Eomeski67,'Type'],Eomes=hipexp['Eomes',Eomeski67],mki67=hipexp['Mki67',Eomeski67])

#条形图
data <- data.frame(table(dbox1$Type))
colnames(data) <- c('Type','number')
data$Type <- factor(data$Type,levels=c('wt','mut'))
p <- ggplot(data = data,aes(x = Type, y = number),width=1)+
  geom_bar(aes(fill = Type),color = 'grey21',stat = 'identity', position = 'dodge')+
  scale_fill_manual(values = c('#4DBBD5','#E64B35'))+theme_bw()+
  theme(panel.grid =element_blank())+
  theme(axis.title.x = element_blank())+
  theme(axis.text=element_text(size=15),axis.title=element_text(size = 15),
  legend.position='none')+
  geom_text(aes(label=number), color="white", 
          size=4,position=position_dodge(0.9),vjust=1.5)
# scale_y_continuous(breaks = c(0,0.5,1))
ggsave(p,filename='Eomes+Ki67+cellnumber.pdf',width = 2.5,height = 3)

#双阳性细胞的表达量
p1 <- ggplot(dbox1, aes(x=Type, y=Eomes,fill=Type)) + 
geom_violin(width=0.6,trim=TRUE,color="#EED8AE",linewidth=0.6,scale = "width",position=position_dodge(width=0.6)) + #小提琴图
geom_jitter(size=0.6,colour='#EED8AE',position=position_jitterdodge(jitter.width = 0.35, jitter.height = 0, dodge.width = 0.6))+#散点
stat_summary(fun.data = 'mean_se',size = 0.6, color = "Tomato",position=position_dodge(0.6))+
scale_fill_manual(values=c('#E64B35','#4DBBD5'))+
theme_bw()+ #背景变为白色
theme(axis.text=element_text(size=15),axis.title=element_text(size = 15), panel.border = element_blank(),
axis.line = element_line(colour = "black",linewidth=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
legend.position = "none",panel.grid.major = element_blank(),panel.background = element_blank(),#不显示网格线
panel.grid.minor = element_blank())+#显示显著性
ylab('Expression of Eomes\nin Eomes+/Ki67+ Cells')+xlab("")+
stat_compare_means(size = 5,label.x=1.4,label = 'p.format',label.y=subset(dbox1,gene==i)[,'Eomes'] %>% max(.)+0.5)

p2 <- ggplot(dbox1, aes(x=Type, y=mki67,fill=Type)) + 
geom_violin(width=0.6,trim=TRUE,color="#EED8AE",linewidth=0.6,scale = "width",position=position_dodge(width=0.6)) + #小提琴图
geom_jitter(size=0.6,colour='#EED8AE',position=position_jitterdodge(jitter.width = 0.35, jitter.height = 0, dodge.width = 0.6))+#散点
stat_summary(fun.data = 'mean_se',size = 0.6, color = "Tomato",position=position_dodge(0.6))+
scale_fill_manual(values=c('#E64B35','#4DBBD5'))+
theme_bw()+ #背景变为白色
theme(axis.text=element_text(size=15),axis.title=element_text(size = 15), panel.border = element_blank(),
axis.line = element_line(colour = "black",linewidth=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
legend.position = "none",panel.grid.major = element_blank(),panel.background = element_blank(),#不显示网格线
panel.grid.minor = element_blank())+#显示显著性
ylab('Expression of Ki67\nin Eomes+/Ki67+ Cells')+xlab("")+
stat_compare_means(size = 5,label.x=1.4,label = 'p.format',label.y=subset(dbox1,gene==i)[,'mki67'] %>% max(.)+0.5)
p <- p1+p2+patchwork::plot_layout(ncol=2,nrow = 1,guides='collect')
ggsave(p,filename='Eomes_Ki67.pdf',width = 5,height = 3)
