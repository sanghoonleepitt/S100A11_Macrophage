---
title: "Cancer-cell derived S100A11 promotes macrophage recruitment in ER+ breast cancer"
author: "SanghoonLee"
date: "2024-03-28"
output: html_document
---

- This code is tested on R v4.3.2 (2023-10-31) and Seurat_5.0.2 
- Platform: aarch64-apple-darwin20 (64-bit), Apple M2 Pro
- Running under: macOS Sonoma 14.4.1

___If your computer spec is not good, your computer may shut down while learning codes.___

To reproduce the figures published in the paper, as we explained in Section 1, you should contact us and obtain Seurat object file,
"SeuratObject_WuetalBassezetalIntegration.rds" and run the rest of the code. 

## Cancer-cell derived S100A11 promotes macrophage recruitment in ER+ breast cancer

Macrophages play a pivotal role in orchestrating breast tumor development, progression, and therapeutic resistance. Estrogen receptor-positive (ER+) breast tumors exhibit diverse patterns of macrophage infiltration and cancer cell-secreted factors that promote a macrophage-rich breast tumor microenvironment that remains poorly understood. 

In this study, we integrated two publicly available single-cell RNA-sequencing datasets of ER+ breast tumors (n=25) to study cancer cell-macrophage interactions and identify predictors of macrophage infiltration.

This study provides new insights into the relationship between macrophages and cancer cell-derived factors in human ER+ breast tumors. Additionally, it uncovers S100A11 as a targetable paracrine regulator of cancer-macrophage interactions in the pro-tumorigenic macrophage-rich breast tumor microenvironment. 


**Installing necessary libraries**

```{InstallLibrary}
rm(list=ls())
library(BiocManager) # CRAN # install.packages("BiocManager")
library(Seurat)  # CRAN # install.packages("Seurat")
library(stringr) # CRAN
library(ggplot2) # CRAN
library(dplyr) # CRAN
library(tidyr) # CRAN
library(data.table) # CRAN
library(dittoSeq)  # Bioconductor - it requires "nloptr" CRAN package # BiocManager::install("dittoSeq")
library(glmGamPoi) # Bioconuctor package for SCTransform.
library(cowplot) # CRAN, for plot_grid()

library(patchwork)# CRAN. to use plot_annotation()
library(harmony) # CRAN, to use RunHarmony() function, https://github.com/immunogenomics/harmony

```

**Setup working directory.**

```{Setup_WorkingDirectory}
dir<-dirname(rstudioapi::getSourceEditorContext()$path); 
setwd(dir); print(dir)   # Your current working directory # [1]
```


**Raw scRNA-seq data are from of Wu et al. and Bassez et al. scRNA-seq data.**

Ref. Wu et al. Nature Genetics 2021, https://www.nature.com/articles/s41588-021-00911-1

Ref. Bassez et al. Nature Medicine 2021, https://www.nature.com/articles/s41591-021-01323-8  


# Analyze scRNA-seq data and make result figures

##   Section 1. UMAP reduction, clustering cells, and define main cell types. 

**To run this code and reproduce figures in the paper, you should contact us and request Seurat object, "SeuratObject_WuetalBassezetalIntegration.rds." We generated it by merging Wu et al. and Bassez et al. scRNA-seq data and performing SCTransform and Harmony integration.**

```{UMAP_Reduction}

SeuratObject_WuetalBassezetal <- readRDS("SeuratObject_WuetalBassezetalIntegration.rds")

################################################################################
## Step1a.	UMAP Reduction and clustering
################################################################################
SeuratObject_WuetalBassezetal <- SeuratObject_WuetalBassezetal %>% RunUMAP(reduction="harmony", dims=1:30, 
                        verbose=F) %>% FindNeighbors(reduction="harmony", k.param=15, dim=1:30)  
SeuratObject_WuetalBassezetal <- SeuratObject_WuetalBassezetal %>% FindClusters(resolution=1.0) %>% identity()
table(SeuratObject_WuetalBassezetal@active.ident) # 50 clusters
# saveRDS(SeuratObject_WuetalBassezetal, file="Seurat_HarmonyUMAPCluster_Res1.0PC30KP15.rds") 

SeuratObject_WuetalBassezetal$CaseID <- SeuratObject_WuetalBassezetal$orig.ident; 
table(SeuratObject_WuetalBassezetal$CaseID ); table(SeuratObject_WuetalBassezetal@active.ident)

### UMAP by cluster ID
UMAP_ByClusterNumber <- DimPlot(SeuratObject_WuetalBassezetal,reduction="umap", label=TRUE, label.size=8) + 
            plot_annotation(title="Two datasets UMAP_Harmony and Clustering number") +
            theme(axis.text.x=element_text(vjust=0.6, size=25,angle=0), 
            axis.text.y=element_text(vjust=0.6, size=25,angle=0))
ggsave(UMAP_ByClusterNumber, height=8,width=11, dpi=300, 
       filename=paste0("UMAPplot_ByClusterNumber_GSE176078EGA6608_MergeSCTHarmony.pdf"), useDingbats=FALSE)

##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  #####
##### $$$$$$ ###### Supplementary Figure 1A. UMAP plot of integrated datasets by CaseID - 
#                                                           Checking Harmony Integration ##### $$$$$ #####
MyDimplot_SeuratMergeSCTHarmony_ByCaseID <- DimPlot(SeuratObject_WuetalBassezetal,reduction="umap",
          group.by="CaseID", label=FALSE) + plot_annotation(title="Two datasets UMAP_Harmony and CaseID") +
          theme(axis.text.x=element_text(vjust=0.6, size=25,angle=0), 
          axis.text.y=element_text(vjust=0.6, size=25,angle=0))
ggsave(MyDimplot_SeuratMergeSCTHarmony_ByCaseID, height=8,width=11, dpi=300, 
       filename=paste0("FigS1A_OutUMAP_WuetalBassezetal_ByCaseID.pdf"), useDingbats=FALSE)

##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  #####
##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  #####
##### $$$$  Figure 2A. UMAP plot of integrated datasets by study ID - Checking Harmony Integration #### $$$$
SeuratObject_WuetalBassezetal$DatasetID <- ifelse( grepl("CID",SeuratObject_WuetalBassezetal$orig.ident),
                       "Wuetal", ifelse( grepl("BIOKEY",SeuratObject_WuetalBassezetal$orig.ident), 
                       "Bassezetal", "NotAvail") )  
SeuratObject_WuetalBassezetal$DatasetID <- factor(SeuratObject_WuetalBassezetal$DatasetID, 
                                                  levels=c("Wuetal","Bassezetal"))
MyDimplot_ByStudyID <- DimPlot(SeuratObject_WuetalBassezetal,reduction="umap", group.by="DatasetID") + 
                  plot_annotation(title="Two datasets after Harmony integration and clustering")
ggsave(MyDimplot_ByStudyID, height=8,width=9, dpi=300, 
              filename=paste0("Fig2A_OutUMAP_MainCellType_WuetalBassezetal_ByStudyName.pdf"), useDingbats=FALSE)
##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  #####

################################################################################
## Step1b. Define cell types by marker gene expression   ####
################################################################################
MyMainMarker <- c("PECAM1","RAMP2","FLT1","CLDN5",   "EPCAM","KRT19","KRT18","CD24",  
                  "PDGFRB","C1R","DCN","COL1A1",   "ACTA2",  "TPSB2","TPSAB1","CPA3",   
                  "CD68","LYZ","TYROBP",   "CD83","MS4A1","MZB1","CD79A",  
                  "CD2","CD3E","CD3D","CD3G","IL7R") 

MyDotPlot_MainType <- DotPlot(SeuratObject_WuetalBassezetal, features=MyMainMarker)+ 
                theme(axis.text.x=element_text(vjust=0.6, size=15,angle=0), 
                axis.text.y=element_text(vjust=0.6, size=15,angle=0)) +
                scale_colour_gradient2(low = "#1515FA", mid = "#FFFAE2", high = "#ff0000") + 
                coord_flip()    # scale_color_viridis_c()   
ggsave(MyDotPlot_MainType, height=10,width=15, dpi=300, 
       filename=paste0("OutDotplot_ByMainMarker_HarmoneyGSE176EGA6608_Res1.0_40Cluster_MoreMarker.pdf"), 
       useDingbats=FALSE)

## FeaturePlot - To pull gene expression cells to front. When features are burried,  it is useful. 
MyFeaturePlot_SingelGene <- Seurat::FeaturePlot(SeuratObject_WuetalBassezetal, raster=TRUE,features=c("EPCAM"), 
                                                cols=c("lightgrey","blue"), order=TRUE, pt.size=1.5) 
# ggsave(MyFeaturePlot_SingelGene, height=3,width=4, dpi=300,
#    filename=paste0("OutFeatureUMAP_ByMainMarker_HarmoneyGSE176EGA6608_Res1.0PC30KP15_EPCAM.pdf"), 
#     useDingbats=FALSE)

MyFeaturePlot_CD14 <- Seurat::FeaturePlot(SeuratObject_WuetalBassezetal, raster=TRUE,features=c("CD14"), 
                                          cols=c("lightgrey","blue"), order=TRUE, pt.size=1.8) 
# ggsave(MyFeaturePlot_CD14, height=3.5,width=4, dpi=300, 
#   filename=paste0("OutFeatureUMAP_ByMainMarker_HarmoneyGSE176EGA6608_Res1.0PC30KP15_CD14.pdf"), 
#   useDingbats=FALSE)

### Set @active.ident with the new cell type definitions. 
new.cluster.ids_Main<- c("Fibroblast","TNKcell", "Endothelial","CancerEpithelial","TNKcell","TNKcell",   
                       "CancerEpithelial","Bcell","Myeloid","CancerEpithelial","Fibroblast",
                       "CancerEpithelial","Fibroblast","CancerEpithelial","TNKcell","CancerEpithelial",   
                       "CancerEpithelial","Fibroblast","CancerEpithelial", "CancerEpithelial","CancerEpithelial",
                       "TNKcell","CancerEpithelial","Endothelial","CancerEpithelial","TNKcell",    
                       "Myofibroblast","Fibroblast","Bcell","Myeloid","Bcell", 
                       "CancerEpithelial","CancerEpithelial","CancerEpithelial","Bcell","CancerEpithelial",    
                       "Fibroblast","Endothelial","CancerEpithelial","TNKcell", "TNKcell",
                       "Mastcell","Mastcell","CancerEpithelial","CancerEpithelial","CancerEpithelial", 
                       "TNKcell","CancerEpithelial","CancerEpithelial","Myofibroblast")
names(new.cluster.ids_Main) <- levels(SeuratObject_WuetalBassezetal)
SeuratObject_WuetalBassezetal <- RenameIdents(SeuratObject_WuetalBassezetal, new.cluster.ids_Main)
table(SeuratObject_WuetalBassezetal@active.ident); sum(table(SeuratObject_WuetalBassezetal@active.ident)) # 72079

## ==== ##### ====== ##### Add the new cell type definition to Metdata  ## ==== ##### ====== ##### 
SeuratObject_WuetalBassezetal <- AddMetaData(object=SeuratObject_WuetalBassezetal, 
          metadata=c(SeuratObject_WuetalBassezetal@active.ident), col.name=c("CellTypeByMarker_GSE176EGA6608"))
saveRDS(SeuratObject_WuetalBassezetal, file="SeuratObject_Harmony_WuetalBassezetal_CellTypeCluster.rds")    

##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  #####
##### $$$$$ ###### Figure 2B. UMAP plot of integrated datasets by cell type definition ##### $$$$$ #####
MyDimplot_CellType <- DimPlot(SeuratObject_WuetalBassezetal, reduction="umap", 
                        group.by="CellTypeByMarker_GSE176EGA6608", 
                        pt.size=.1, label=FALSE, label.size=8) + NoLegend(); # label=TRUE,
ggsave(MyDimplot_CellType, height=7,width=7, dpi=300, 
       filename=paste0("Fig2B_OutUMAP_MainCellType_WuetalBassezetal_NoLabel.pdf"), useDingbats=FALSE)

##### $$$$$  ###### Figure 2C. UMAP Feature plot by marker gene expression.  ###### $$$$$ #####
MyFeaturePlot <- FeaturePlot(SeuratObject_WuetalBassezetal,raster=TRUE, cols=c("lightgrey","darkblue"), 
                  order=TRUE, pt.size=1.8, features=c("EPCAM","KRT19","KRT18","CD14", "PDGFRB","ACTA2", 
                      "CD3D","CD68","MS4A1", "PECAM1","IL7R",'CD8A'))   #  raster.dpi = c(512, 512)   
ggsave(MyFeaturePlot, height=12.1,width=18, dpi=300, 
  filename=paste0("Fig2C_OutFeatureUMAP_WuetalBassezetal_ByMainMarkergeneExpression.pdf"), useDingbats=FALSE)
##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  #####

## Metadata has the cluster IDs already. The last column, "seurat_clusters" has the cluster IDS
MyMainMetadata <- SeuratObject_WuetalBassezetal@meta.data %>% 
                  tibble::rownames_to_column("CellID_GSE176EGA6608") %>% 
                  dplyr::select(-c(orig.ident,nCount_RNA,nFeature_RNA,percent.mt,SCT_snn_res.1));
dim(MyMainMetadata) 
table(MyMainMetadata$seurat_clusters)
table(MyMainMetadata$CellTypeByMarker_GSE176EGA6608)
MyMainMetadata[1:2,]

### Violin plot by marker genes. 
MyViolinplot_Main <- VlnPlot(object=SeuratObject_WuetalBassezetal, features=c("CPA3",
                     "PECAM1","EPCAM","KRT19","PDGFRB","ACTA2","CD68","MS4A1","CD3D"), stack=TRUE,flip=TRUE )
ggsave(MyViolinplot_Main, height=6,width=9, dpi=300, 
       filename=paste0("OutViolin_MainCellTypeMarker_HarmoneyGSE176EGA6608_Res1.0PC30KP15_ERp.pdf"), 
       useDingbats=FALSE)

```

##  Section 2. Subset Myeloid cells and find macrophage cells by marker genes expression

```{MyeloidSubsuet}
################################################################################
## Step2a. Subset Myeloid and define cell subtypes
################################################################################
# Subset by active.ident
SeuratObject_WuetalBassezetal_Myeloid <- subset(x=SeuratObject_WuetalBassezetal, idents =c("Myeloid")) 
length(SeuratObject_WuetalBassezetal_Myeloid@active.ident)  # 3041

HighVarGene <- SeuratObject_WuetalBassezetal_Myeloid@assays$SCT@var.features;length(HighVarGene) 
SeuratObject_WuetalBassezetal_Myeloid <- ScaleData(SeuratObject_WuetalBassezetal_Myeloid, features=HighVarGene) 
SeuratObject_WuetalBassezetal_Myeloid<- Seurat::RunPCA(SeuratObject_WuetalBassezetal_Myeloid, verbose = FALSE)   
SeuratObject_WuetalBassezetal_Myeloid <- FindNeighbors(SeuratObject_WuetalBassezetal_Myeloid, reduction = "pca", 
                                                       dims = 1:30)
SeuratObject_WuetalBassezetal_Myeloid <- FindClusters(SeuratObject_WuetalBassezetal_Myeloid, resolution = 1.2) 
table(SeuratObject_WuetalBassezetal_Myeloid@active.ident)

# saveRDS(SeuratObject_WuetalBassezetal_Myeloid, file="SeuratObj_WuetalBassezetal_Myeloid.rds")

## By TSNE
SeuratObject_WuetalBassezetal_Myeloid <- RunTSNE(SeuratObject_WuetalBassezetal_Myeloid, dims = 1:30)  
# lower dims numbers will give less number of clusters image (less white space)
MyDimplot_ERpos_Myeloid_TSNE <- DimPlot(SeuratObject_WuetalBassezetal_Myeloid, pt.size=1, 
                                        label.size=10,reduction="tsne", label=TRUE)
ggsave(MyDimplot_ERpos_Myeloid_TSNE, height=8,width=9, dpi=300, 
      filename=paste0("OutTSNEplot_GSE176078_MacrophageDC_Res1.2_21clst_FromWholeRes1.0PC30KP15_ERp.pdf"), 
      useDingbats=FALSE)

################################################################################
## Step2b. Dotplot by Myeloid markers
################################################################################
MyFeature_Myeloid <- c( "IL1B","CSF1R","S100A9","FCGR3A",   "CD14","SIGLEC1","CXCL10","EGR1","CD68","FCGR1A",
                        "FABP5","APOE",  "CD80","CD86","CD163","MRC1","MSR1",    "CLEC10A","THBD","CD1C",
                        "ITGAX","HLA-DRB1","LAMP3","CLEC9A","GZMB","FLT3",  "IL3RA","SELL","IRF7","TYROBP", 
                        "LILRA4","CXCR3","CLEC4C",   "IL10","CD40",  "CCR2","CCL2","CCL18","MMP9","CX3CR1",
                        "MT1G","SLC2A1","LYVE1","LYZ",   "ACE", "ADGRE1",  "TEK") 
                          # "ADGRE1","TEK", don't have expression. 
MyFeature_Myeloid_Cheng <- c("LST1","LILRB2","FCN1","CD14", "PTPRC","CD93","S100A8","S100A9","CD86",    
                             "C1QA","C1QC","SEPP1","PLTP","EREG","NLRP3","CCL4","CD68","FCGR3A","MRC1",
                             "IDO1","CCR7","FSCN1","SELL","IRF7")
                        
MyDotPlot_Myeloid <- DotPlot(SeuratObject_WuetalBassezetal_Myeloid, features=MyFeature_Myeloid_Cheng)+ 
              theme(axis.text.x=element_text(vjust=0.6, angle=90)) + 
              scale_colour_gradient2(low = "#1515FA", mid = "#FFFAE2", high = "#ff0000") + coord_flip() 

ggsave(MyDotPlot_Myeloid, height=8,width=8, dpi=300, 
       filename=paste0("OutDotplot_ByMyeloidMarker_GSE176EGA6608_Res1.0_21clst_FromWholecellRes1.0.pdf"),
       useDingbats=FALSE)

## Metadata has the cluster IDs already. The last column, "seurat_clusters" has the cluster IDS
MyMetadata <- SeuratObject_WuetalBassezetal_Myeloid@meta.data %>% tibble::rownames_to_column("CellID_") %>% 
              dplyr::select(-c(orig.ident,nCount_RNA,nFeature_RNA,percent.mt,SCT_snn_res.1,SCT_snn_res.1.2));
dim(MyMetadata) # 3658   14
table(MyMetadata$seurat_clusters)
# 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17 
# 382 309 289 282 271 260 227 179 157 134 130  93  84  80  58  51  31  24 

## Macrophage assignment by Wu et al. paper.
MyMetadata_Mac <- MyMetadata %>% dplyr::filter(CellTypeMinor=="Macrophage", ); dim(MyMetadata_Mac)# 1578 27 
table(MyMetadata_Mac$seurat_clusters)  # 0, 7, 8, 10, 16 are macrophage, maybe 9 is not macrophage

## from whole cell res 1.5  - Best
new.cluster.ids_Myeloid <- c("Macrophage","Macrophage","Macrophage","Macrophage","Macrophage","Macrophage",   
                             "Monocyte","Macrophage","Macrophage","Monocyte","Macrophage",   
                             "DC","Macrophage","Macrophage","Macrophage","Monocyte",             
                             "DC","Monocyte")  
names(new.cluster.ids_Myeloid) <- levels(SeuratObject_WuetalBassezetal_Myeloid)
SeuratObject_WuetalBassezetal_Myeloid <- RenameIdents(SeuratObject_WuetalBassezetal_Myeloid, 
                                          new.cluster.ids_Myeloid)
table(SeuratObject_WuetalBassezetal_Myeloid@active.ident)
# Macrophage   Monocyte         DC 
# 2481        436        124

MeyloidSubsetMetadata <- SeuratObject_WuetalBassezetal_Myeloid@meta.data[, c("CaseID","seurat_clusters" )]
MeyloidSubtypeDF <- data.frame(SeuratObject_WuetalBassezetal_Myeloid@active.ident); 
colnames(MeyloidSubtypeDF)[1]<-"MyeloidSubtype"
MeyloidClusterNumbSubtype <- cbind(MeyloidSubsetMetadata,MeyloidSubtypeDF )
# fwrite(MeyloidClusterNumbSubtype, file="MeyloidClusterNumbSubtype_ERp.txt",
#                           row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

## By TSNE
SeuratObject_WuetalBassezetal_Myeloid <- RunTSNE(SeuratObject_WuetalBassezetal_Myeloid, dims = 1:30) 
# lower dims numbers will give less number of clusters image (less white space)
MyDimplot_ERpos_Myeloid_Rename_TSNE <- DimPlot(SeuratObject_WuetalBassezetal_Myeloid, pt.size=1, 
                                               label.size=8, reduction="tsne", label=TRUE)
ggsave(MyDimplot_ERpos_Myeloid_Rename_TSNE, height=8,width=9, dpi=300, 
       filename=paste0("OutTSNE_ByMyeloidMarker_GSE176EGA6608_Res1.2_21clst_FromWholeRes1.0PC30KP15_ERp.pdf"), 
       useDingbats=FALSE)

## Cell type annotation by Wu et al. and Bassez et al. 
SeuratObject_WuetalBassezetal_Myeloid$CellType_WuBassez<-
                            ifelse(is.na(SeuratObject_WuetalBassezetal_Myeloid$cellType),
                            paste0(SeuratObject_WuetalBassezetal_Myeloid$CellTypeMinor, "_Wu"), 
                            paste0(SeuratObject_WuetalBassezetal_Myeloid$cellType, "_Bassez" ) ) #
MyDimplot_ERpos_Myeloid_Rename_TSNE_Wu <- DimPlot(SeuratObject_WuetalBassezetal_Myeloid, pt.size=1, 
                            label.size=8, reduction="tsne", label=TRUE, group.by="CellType_WuBassez")
#ggsave(MyDimplot_ERpos_Myeloid_Rename_TSNE_Wu, height=8,width=9, dpi=300, 
#    filename=paste0("OutTSNE_ByMyeloidMarker_GSE176EGA6608_WuBassezCellType.pdf"), useDingbats=FALSE)
table(SeuratObject_WuetalBassezetal_Myeloid@active.ident,SeuratObject_WuetalBassezetal_Myeloid$CellType_WuBassez)

MyFeaturePlot_Myeloid <- FeaturePlot(SeuratObject_WuetalBassezetal_Myeloid,raster=FALSE,reduction="tsne",
                  pt.size=0.5,features=c("CD68","FCGR3A","CD14","CD86","MRC1","LAMP3","CLEC9A","SELL","IRF7"))  
                  #  raster.dpi = c(512, 512)   
# ggsave(MyFeaturePlot_Myeloid, height=10,width=10, dpi=300, 
#                    filename=paste0("OutFeatureUMAP_Myeloid_Res1.0PC30KP15_Res1.2_ERp.pdf"), useDingbats=FALSE)

##### $$$$$  ###### Supplementary Figure 1B. Feature violin plot by marker gene expression for Myeloid 
MyViolinplot <- VlnPlot(object = SeuratObject_WuetalBassezetal_Myeloid, 
                        features = c("CD68","FCGR3A","CD14","CD86","MRC1","CD1C","LAMP3","CLEC9A","SELL","IRF7"),
                        pt.size=2, stack=TRUE, flip=TRUE )
ggsave(MyViolinplot, height=5,width=5, dpi=300, filename=paste0("FigS1B_OutViolinplot_MyeloidMarkerExp.pdf"), 
                      useDingbats=FALSE)
#####  ============== $$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  #####  ============== $$$$$$$$$$  

################################################################################
## Step2c. Store the cell type identity as a data.frame
################################################################################
### Store the cell type identity as a data.frame
CellTypeMyeloidDataFrame <- data.frame((SeuratObject_WuetalBassezetal_Myeloid@active.ident)); 
dim(CellTypeMyeloidDataFrame) # 3041 1
colnames(CellTypeMyeloidDataFrame)[1]<-"MyeloidType"; 
CellTypeMyeloidDataFrame_Proc <- CellTypeMyeloidDataFrame %>% tibble::rownames_to_column("CellID_GSE176EGA6608"); 
head(CellTypeMyeloidDataFrame_Proc); table(CellTypeMyeloidDataFrame_Proc$MyeloidType)
#               CellID_GSE176EGA6608 MyeloidType
# 1 GSE176_CID3586_AACCATGCAGGTCGTC  Macrophage
# 2 GSE176_CID3586_AACTTTCGTGACCAAG  Macrophage

```

##   Section 3. Subset NKT cells and find cell sub types by marker genes expression 
```{NKT}
################################################################################
## Step3a. Subset NKT cells 
################################################################################
SeuratObject_WuetalBassezetal_NKT <- subset(x=SeuratObject_WuetalBassezetal, idents =c("TNKcell")) 
length(SeuratObject_WuetalBassezetal_NKT@active.ident)  # 16956

# Normalize => Find Variable Features
HighVarGene <- SeuratObject_WuetalBassezetal_NKT@assays$SCT@var.features;length(HighVarGene)  # 2000
SeuratObject_WuetalBassezetal_NKT <- ScaleData(SeuratObject_WuetalBassezetal_NKT, features=HighVarGene) 
SeuratObject_WuetalBassezetal_NKT<- Seurat::RunPCA(SeuratObject_WuetalBassezetal_NKT, verbose = FALSE) 

SeuratObject_WuetalBassezetal_NKT <- FindNeighbors(SeuratObject_WuetalBassezetal_NKT, reduction="pca",
                                                   dims = 1:30)
SeuratObject_WuetalBassezetal_NKT <- FindClusters(SeuratObject_WuetalBassezetal_NKT, resolution=1.0)  
table(SeuratObject_WuetalBassezetal_NKT@active.ident)

# saveRDS(SeuratObject_WuetalBassezetal_NKT, 
#          file="SeuratObj_HarmonyGSE176EGA6608s_Tcell_PCA_AfQC_FromWholeRes1.0PC30KP15_Res1.0_ERp.rds")

SeuratObject_WuetalBassezetal_NKT <- RunTSNE(SeuratObject_WuetalBassezetal_NKT, dims = 1:30)  
# lower dims numbers will give less number of clusters image (less white space)
MyDimplot_ERpos_Tcell_TSNE <- DimPlot(SeuratObject_WuetalBassezetal_NKT, pt.size=1, 
                                      label.size=8,reduction="tsne", label=TRUE)
# ggsave(MyDimplot_ERpos_Tcell_TSNE, height=8,width=8, dpi=300, 
#         filename=paste0("OutTSNEplot_GSE176EGA6608_Tcell_Res1.0_31clst_FromWholeRes1.0PC30KP15_ERp.pdf"), 
#         useDingbats=FALSE)

################################################################################
## Step3b. Dotplot by NK, CD4+T, CD8+T cell markers
################################################################################
Marker_Tcell <- c("FCGR3A","IL2RB","KLRB1","KLRC1","KLRD1","KLRF1","KLRK1","NCR3",
                  "NCR1","FGFBP2","AREG","XCL1","KIR2DL4","NCAM1",
                  "CD8A","CD8B","GZMA","GZMB",  "LAG3","PDCD1","CTLA4",   
                  "TNFRSF4","BTLA","CD40LG","STAT4","STAT1","FOXP3","CD4","CCR7",   
                  "IL7R","IL2RA","PTPRC","CD3D","CD27")   # , 

MyDotPlot_Tcell <- DotPlot(SeuratObject_WuetalBassezetal_NKT, features=Marker_Tcell)+ 
  theme(axis.text.x=element_text(vjust=0.6, angle=90)) +
  scale_colour_gradient2(low = "#1515FA", mid = "#FFFAE2", high = "#ff0000") + coord_flip()
# ggsave(MyDotPlot_Tcell, height=8,width=8, dpi=300, 
#     filename=paste0("OutDotplot_ByTcellMarker_GSE176EGA6608_Res0.8_17clst_FromWholecellRes1.5_ERp.pdf"), 
#     useDingbats=FALSE)

## Metadata has the cluster IDs already. The last column, "seurat_clusters" has the cluster IDS
MyMetadata <- SeuratObject_WuetalBassezetal_NKT@meta.data %>% tibble::rownames_to_column("CellID_") %>% 
                    dplyr::select(-c(orig.ident,nCount_RNA,nFeature_RNA,percent.mt,SCT_snn_res.1));
dim(MyMetadata) # 16575    14
table(MyMetadata$seurat_clusters) # 30 clusters

#### =========== New cell type assignment   ============= ########### 
new.cluster.ids_Tcell <- c("NKcell","CD8Tcell","CD4Tcell","CD8Tcell","CD8Tcell","NKTcell",    
                           "CD4Tcell","CD4Tcell","Treg","NKTcell","NotAvail",
                           "CD8Tcell","CD4Tcell","NKTcell","NKTcell","NotAvail",         
                           "CD4Tcell","NKTcell","CD4Tcell","CD4Tcell","CD8Tcell",
                           "NKcell","NKcell","CD4Tcell","CD8Tcell","NotAvail",       
                           "CD4Tcell","CD8Tcell","CD4Tcell","NKcell")
length(new.cluster.ids_Tcell)
names(new.cluster.ids_Tcell) <- levels(SeuratObject_WuetalBassezetal_NKT)
SeuratObject_WuetalBassezetal_NKT <- RenameIdents(SeuratObject_WuetalBassezetal_NKT, new.cluster.ids_Tcell)
table(SeuratObject_WuetalBassezetal_NKT@active.ident); 
sum(table(SeuratObject_WuetalBassezetal_NKT@active.ident))  # 17103
# NKcell CD8Tcell CD4Tcell  NKTcell     Treg NotAvail 
# 2141     5084     4877     3040      732     1229 
# saveRDS(SeuratObject_WuetalBassezetal_NKT, "SeuratObject_WuetalBassezetal_NKT.rds")

# SeuratObject_WuetalBassezetal_NKT <- subset(x=SeuratObject_WuetalBassezetal_NKT, idents = c("NotAvail"),
#                                   invert=TRUE) # invert=TRUE willl exclude the ident cells.
table(SeuratObject_WuetalBassezetal_NKT@active.ident); 
sum(table(SeuratObject_WuetalBassezetal_NKT@active.ident)) # 17103
# NKcell CD8Tcell CD4Tcell  NKTcell     Treg NotAvail 
# 2141     5084     4877     3040      732     1229 

SeuratObject_WuetalBassezetal_NKT <- RunTSNE(SeuratObject_WuetalBassezetal_NKT, dims = 1:30)  
# lower dims numbers will give less number of clusters image (less white space)
MyDimplot_ERpos_Tcell_Rename_TSNE <- DimPlot(SeuratObject_WuetalBassezetal_NKT,label.size=8, 
                                             pt.size=1, reduction="tsne", label=FALSE)
# ggsave(MyDimplot_ERpos_Tcell_Rename_TSNE, height=8,width=8.5, dpi=300, 
# filename=paste0("OutTSNE_ByTcellMarker_GSE176EGA6608_Res1.0_31clst_FromWholeRes1.5PC30KP15_ERp.pdf"),
# useDingbats=FALSE)

##### $$$$$$ ###### Supplementary Figure 1C. Feature violin plot by marker gene expression for NK/T cell   
SeuratObject_WuetalBassezetal_NKT@active.ident <- factor(SeuratObject_WuetalBassezetal_NKT@active.ident, 
                                          levels=c("NKcell","NKTcell","CD4Tcell","CD8Tcell","Treg","NotAvail"))
MyViolinplot_NKTcell <- VlnPlot(object = SeuratObject_WuetalBassezetal_NKT, 
                                features = c("AREG","KLRB1", "XCL1", "IL7R","IL2RB","CD4","CD8A","CD8B","FOXP3"), 
                                stack=TRUE, flip=TRUE)
ggsave(MyViolinplot_NKTcell, height=4.5,width=5, dpi=300, 
                filename=paste0("FigS1C_OutViolin_NKTcellMarkerExp.pdf"), useDingbats=FALSE)
#####  ============== $$$$$$$$$$$$$  ##### ============== $$$$$$$$$$$$$  ##### $$$$$$$$$$$$$  

#################################################################################
## Step3c. Store the cell type identity as a data.frame
#################################################################################
CellTypeTcellDataFrame <- data.frame((SeuratObject_WuetalBassezetal_NKT@active.ident)); 
dim(CellTypeTcellDataFrame) # 17103   1
colnames(CellTypeTcellDataFrame)[1]<-"TcellType"; 
CellTypeTcellDataFrame_Proc <- CellTypeTcellDataFrame %>% tibble::rownames_to_column("CellID_GSE176EGA6608");
head(CellTypeTcellDataFrame_Proc)
table(CellTypeTcellDataFrame_Proc$TcellType); dim(CellTypeTcellDataFrame_Proc)  #  17103     2
# NKcell  NKTcell CD4Tcell CD8Tcell     Treg NotAvail 
# 2141     3040     4877     5084      732     1229 
```


##Section 4.  Annotate Myeloid cell subtypes and NKT cell subtypes, and calculate Macrophage fraction per samples. 

```{CellTypeAnnotation}
#################################################################################
## Step4a. Annotate Myeloid cell subtypes and NKT cell subtypes  
#################################################################################
CellTypeClusterDataFrame_Proc <- MyMainMetadata[, c("CellID_GSE176EGA6608","CellTypeByMarker_GSE176EGA6608")]
CellTypeClusterMyeloidTcell <- dplyr::left_join(CellTypeClusterDataFrame_Proc, CellTypeMyeloidDataFrame_Proc) %>% 
          dplyr::left_join(CellTypeTcellDataFrame_Proc) # join_by(CellID_GSE176EGA6608)
dim(CellTypeClusterMyeloidTcell) # [1] 72079  5
CellTypeClusterMyeloidTcell[1:2,]

table(CellTypeClusterMyeloidTcell$CellTypeByMarker_GSE176EGA6608); 
sum(table(CellTypeClusterMyeloidTcell$CellTypeByMarker_GSE176EGA6608)) # 72079
# Fibroblast          TNKcell      Endothelial CancerEpithelial     Bcell      Myeloid    Myofibroblast     Mastcell 
# 12503            17103             6421            27014         4460        3041           987         550 

CellTypeClusterMyeloidTcell$ByMainMarker_MyeloidTcell <- 
                            ifelse(CellTypeClusterMyeloidTcell$CellTypeByMarker_GSE176EGA6608 == "Myeloid",
                            as.vector(CellTypeClusterMyeloidTcell$MyeloidType),
                            ifelse(CellTypeClusterMyeloidTcell$CellTypeByMarker_GSE176EGA6608 == "TNKcell",
                            as.vector(CellTypeClusterMyeloidTcell$TcellType),
                            as.vector(CellTypeClusterMyeloidTcell$CellTypeByMarker_GSE176EGA6608)))
table(CellTypeClusterMyeloidTcell$ByMainMarker_MyeloidTcell); 
sum(table(CellTypeClusterMyeloidTcell$ByMainMarker_MyeloidTcell))


################################################################################
## Step4b. Add new metadata of new cell subtypes to Seurat object
################################################################################
SeuratObject_WuetalBassezetal <- AddMetaData(object=SeuratObject_WuetalBassezetal, 
                                            metadata=c(CellTypeClusterMyeloidTcell$ByMainMarker_MyeloidTcell), 
                                            col.name=c("CellTypeMacroTcell_GSE176EGA6608"))
table(SeuratObject_WuetalBassezetal$CellTypeMacroTcell_GSE176EGA6608) ## New celltype

# saveRDS(SeuratObject_WuetalBassezetal, file="SeuratObj_Harmony_WuetalBassezetal_CelltypeDefined.rds")

SeuratObject_WuetalBassezetal <- SetIdent(SeuratObject_WuetalBassezetal, 
                                  value=SeuratObject_WuetalBassezetal$CellTypeMacroTcell_GSE176EGA6608)
table(SeuratObject_WuetalBassezetal@active.ident)

################################################################################
## Step4c. Make violin plots by main marker gene expression 
################################################################################
SeuratObject_WuetalBassezetal <- subset(x=SeuratObject_WuetalBassezetal, idents = c("NotAvail"), invert=TRUE) 
# invert=TRUE willl exclude the ident cells.  #  I alrady removed "NotAvail" cells  
table(SeuratObject_WuetalBassezetal@active.ident)
table(SeuratObject_WuetalBassezetal$CellTypeMacroTcell_GSE176EGA6608);
sum(table(SeuratObject_WuetalBassezetal$CellTypeMacroTcell_GSE176EGA6608)) # 70850

CellTypeClusterMyeloidTcell_Rmv <- CellTypeClusterMyeloidTcell %>% 
                        dplyr::filter(!ByMainMarker_MyeloidTcell %in% c("NotAvail"))

table(CellTypeClusterMyeloidTcell_Rmv$ByMainMarker_MyeloidTcell)  # This has NO "NotAvail)
saveRDS(SeuratObject_WuetalBassezetal, file="SeuratObj_Harmony_WuetalBassezetal_CelltypeDefined_RmvNotAvail.rds")

### Violin plot by marker genes per main cell types 
MyViolinplot_MainSub <- VlnPlot(object=SeuratObject_WuetalBassezetal, 
                        features=c("FOXP3", "CPA3","PECAM1","EPCAM","KRT19","PDGFRB","ACTA2",
                        "FCGR1A","CD163","CD68","CD14","CD86", "MS4A1","CD3D","CD8A","IL7R","IRF7","IL2RB"), 
                        group.by="CellTypeByMarker_GSE176EGA6608", stack=TRUE,flip=TRUE )
ggsave(MyViolinplot_MainSub, height=6,width=7, dpi=300, 
       filename=paste0("OutViolin_MainSubCellTypeMarker_HarmoneyGSE176EGA6608_Res1.5PC30KP15_20230819.pdf"), 
       useDingbats=FALSE)

#Violin plot by marker genes per sub celltypes <= I should load "SeuratObject_WuetalBassezetal"first at line 678
SeuratObject_WuetalBassezetal$CellTypeMacroTcell_GSE176EGA6608 <- 
                 factor(SeuratObject_WuetalBassezetal$CellTypeMacroTcell_GSE176EGA6608, 
                 levels=c("Bcell","DC","Mastcell", "Macrophage", "Monocyte", "NKcell","NKTcell",
                 "CD4Tcell","CD8Tcell","Treg","CancerEpithelial","Endothelial","Fibroblast","Myofibroblast"))
MyViolinplot_MainSub <- VlnPlot(object=SeuratObject_WuetalBassezetal, 
                    features=c("FOXP3", "CPA3","PECAM1","EPCAM","KRT19","PDGFRB","ACTA2","FCGR1A","CD68","CD14", 
                               "MS4A1","CD3D","CD8A", "IL7R","IRF7","IL2RB", "XCL1","KLRD1","KLRB1", "AREG"), 
                      group.by="CellTypeMacroTcell_GSE176EGA6608", stack=TRUE,flip=TRUE )
ggsave(MyViolinplot_MainSub, height=6.5,width=6.5, dpi=300, 
       filename=paste0("OutViolin_MainSubCellTypeMarker_HarmoneyGSE176EGA6608_Res1.5PC30KP15_20240222.pdf"), 
       useDingbats=FALSE)

## Dotplot visualizing averaged expression of cannonical markers in cell clusters. 
## Line513 may cause an error,'Error in `.rowNamesDF<-`(x, value = value):duplicate 'row.names' are not allowed'. 
## According to your Seurat package version, you may pass or not pass this step. 
## If you get an error, email me "Sanghoon Lee, sal170@pitt.edu" and ask for 
## "SeuratObj_Harmony_WuetalBassezetal_CelltypeDefined_RmvNotAvail.rds" file at line 678
## When you read this rds file and run the following lines, you can get the dotplot without an error. 
MyMainMarker <- c("PECAM1","RAMP2","FLT1","CLDN5", "EPCAM","KRT19","KRT18","CD24", "PDGFRB","C1R","DCN","COL1A1",
                "ACTA2",  "TPSB2","TPSAB1","CPA3", "CD68","LYZ","TYROBP",  "LILRA4","LAMP3",  
                "CLEC9A","STAT1","IRF7","CD83","MS4A1","MZB1","CD79A", "CD8A","CD2","CD3E","CD3D","CD3G","IL7R") 
MyDotPlot_MainType <- DotPlot(SeuratObject_WuetalBassezetal, 
                       features=MyMainMarker, group.by="CellTypeMacroTcell_GSE176EGA6608")+ 
                       theme(axis.text.x=element_text(vjust=0.6, size=15,angle=0), 
                       axis.text.y=element_text(vjust=0.6, size=15,angle=0)) +
                       scale_colour_gradient2(low = "#1515FA", mid = "#FFFAE2", high = "#ff0000") + coord_flip()    
# ggsave(MyDotPlot_MainType, height=10,width=15, dpi=300, 
#       filename=paste0("OutDotplot_ByMainMarker_HarmoneyGSE176EGA6608_Res1.0_40Cluster_MoreMarker_20240222.pdf"),
#       useDingbats=FALSE)

```

# E. Calculate cell tpy fraction and make necessary figures

##   Section 5.  Calculate cell type fraction per samples. 

Make dotplot of Macrophage fraction per sample. Make fraction barplot for all samples 

```{CellTypeFraction}
################################################################################
## Step5a. Cell type fraction per samples. Make dotplot for macrophage fraction. 
## number of cells per cluster: https://github.com/satijalab/seurat/issues/2825
################################################################################
CellTypeFraction<- with(SeuratObject_WuetalBassezetal@meta.data, 
                        table(orig.ident, CellTypeMacroTcell_GSE176EGA6608)); 
CellTypeFraction[1:3,]  # class is table.  # 14 10

####   To make input data for cell type fraction correlation heatmap  ##################
CellTypeFraction_DF <- CellTypeFraction %>% matrix(ncol=ncol(CellTypeFraction)) %>% data.frame
colnames(CellTypeFraction_DF) <- colnames(CellTypeFraction)
rownames(CellTypeFraction_DF) <- rownames(CellTypeFraction)
head(CellTypeFraction_DF)

ProportionCalculation <- CellTypeFraction_DF/rowSums(CellTypeFraction_DF);dim(ProportionCalculation) 
saveRDS(ProportionCalculation, "ProportionCalculation.rds")

library(ie2misc) # CRAN, for madstat() ##   error: X11 library is missing: install XQuartz from www.xquartz.org
cell_type <-"Macrophage" #  "Monocyte"
MacroFractMean <- mean(ProportionCalculation[,cell_type]*100); print(MacroFractMean)   # 4.917426
MacroFractMedian <- median(ProportionCalculation[,cell_type]*100);  print(MacroFractMedian)  # 4.217754
MacroFractSTD <- sqrt(var(ProportionCalculation[,cell_type]*100));  print(MacroFractSTD)  # 2.725528
MacroFractMAD <- sqrt(madstat(ProportionCalculation[,cell_type]*100));  print(MacroFractMAD) # 1.46275

dplyr::ntile(sort(ProportionCalculation[,cell_type])*100, 4)  
sort(ProportionCalculation[,cell_type])

##### -----------   ##### ----------- ============== $$$$$$$$$$$$$  ##### ----- $$$$$$$$$$$$$  
#### $$$$ #### Supplementary Figure 1D. Dotplot to represent Macrophage infiltration in individual samples 
pdf(file="FigS1D_OutDotplot_MacrophageInfiltration_IndividualSample.pdf", width=5, height=5)
    plot(sort(ProportionCalculation[,cell_type])*100, pch=16, cex=1.3)
dev.off()
##### ----------- ============== $$$$$$$$$$$$$  ##### ----------- ============== $$$$$$$$$$$$$  

#################################################################################################
### Step5b. Make fraction barplot that represent each cell type fraction in different samples
#################################################################################################
NewIndex<-c()
#for (EachCase in names(sort(ProportionCalculation[,cell_type]))) { 
for (EachCase in rownames(ProportionCalculation)[order(ProportionCalculation[,cell_type])] ) { 
  Index<-print(grep(paste0(EachCase, "$"), rownames(ProportionCalculation) ) )
  NewIndex<-c(NewIndex, Index)
}

table(SeuratObject_WuetalBassezetal$CellTypeMacroTcell_GSE176EGA6608)


##### -----------  ============== $$$$$$$$$$$$$  ##### ---------- ============== $$$$$$$$$$$$$  
##### $$$ ####  Figure 2D. Barplot to represent cell type infiltration in individual samples #### $$$$$ ####
FractionBarplot <- dittoBarPlot(object = SeuratObject_WuetalBassezetal, var="CellTypeMacroTcell_GSE176EGA6608", 
              group.by="orig.ident", x.reorder=NewIndex, var.labels.reorder=c(8,1,2,3,4,5,6,7,9,10,11,12,13, 14),
              color.panel=c("#999999", "green","orange", "magenta","#ff3300",    
              "yellow","#66ffff","#ff6633", "#0099ff",  "black","#9966ff","#339900" ,"#99FF99" ,"blue"))  
               #  x.reorder = seq(1,20,1)) x.reorder=c(5,3,4,1,2)
ggsave(FractionBarplot, height=5,width=8, dpi=300, 
       filename=paste0("Fig2D_OutFractionBarplot_CellTypeInfiltration.pdf"), useDingbats=FALSE)
##### ----------- ============== $$$$$$$$$$$$$  ##### -----------  ----------- ============== $$$$$$$$$$$$$  

#################################################################################################
### Step5c. Make a table to record cell type fraction in different samples. 
#################################################################################################
ProportionCalculationProcess <- ProportionCalculation*100
ProportionCalculationProcess$MacroPoorMidRich <- ifelse(ProportionCalculationProcess[, "Macrophage"] >= 4, 
                               "MacroRich", ifelse(ProportionCalculationProcess[, "Macrophage"] <= 2, 
                               "MacroPoor","MacroMid"))
ProportionCalculationMacroRatio <- ProportionCalculationProcess %>% tibble::rownames_to_column("CaseID"); 
colnames(ProportionCalculationMacroRatio)[ncol(ProportionCalculationMacroRatio)] <- "MacrophageRatio"
table(ProportionCalculationMacroRatio$MacrophageRatio)
ProportionCalculationMacroRatio$CaseID[ProportionCalculationMacroRatio$MacrophageRatio=="MacroMid"] 
fwrite(ProportionCalculationMacroRatio, file="MacrophageInfiltrarion_IndividualSample.txt", 
           col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE) 
          ## Use this for celltype fraction correlation heatmap.

#################################################################################################
### Step5d. Cell type correlation pyramid plot.
#################################################################################################
library(ComplexHeatmap)   #BiocManager::install("ComplexHeatmap")
library(circlize)  # cran  for colorRamp2 function
library(corrplot)
options(rgl.useNULL = TRUE)  # https://stackoverflow.com/a/66127391/2554330

Correlation.corplot_SH<-function(corMatrix,fontsize=0.8){
        rownames(corMatrix)=gsub("_"," ",rownames(corMatrix))
        colnames(corMatrix)=gsub("_"," ",colnames(corMatrix))
        mycolor=colorRampPalette(c("darkblue","blue", "white","orange", "red","red"))(n=100)
        corrplot(corMatrix, method="circle",type="lower",tl.col="black",col=mycolor,tl.cex=fontsize,tl.srt=30)
}

CellTypeProportionData <- ProportionCalculationMacroRatio %>% tibble::column_to_rownames("CaseID") %>% 
                  dplyr::select(-c(MacrophageRatio ))

## Calculate Spearman correlation.   === Calculate correlation with continuous values, not count data. 
CellTypeProportionData_Correlation <- CellTypeProportionData [,c(1,  6,7,8,9,10,   2,3,4,5,11,12,13,14)]
colnames(CellTypeProportionData_Correlation)

CellTypeProportionData_Cor <- cor(CellTypeProportionData_Correlation, method="spearman"); 
dim(CellTypeProportionData_Cor); CellTypeProportionData_Cor[1:2,] 

pdf(file="Fig2E_PyramidCorrelation_ByCellTypeFraction.pdf",width=6, height=6)
         Correlation.plot=Correlation.corplot_SH(corMatrix=CellTypeProportionData_Cor)
graphics.off()
```


## Section 6. Scaled by DefaultAssay RNA - this is used to calculate 
##            correlation between gene expression and Macrophage fraction. 

```{Scaling}
#################################################################################################
## Step6a. Scaled by DefaultAssay RNA 
#################################################################################################
library(patchwork) # plot_annotation function
library(gridExtra)  # for 'grid.arrange' function
library(ggrepel) # for geom_text_repel function. 

### Data scaling.
DefaultAssay(SeuratObject_WuetalBassezetal) <- "RNA"
TotalGeneNumb<-nrow(SeuratObject_WuetalBassezetal[["RNA"]]@features); print(TotalGeneNumb)  # 35929
SeuratObject_WuetalBassezetal_Normal <- SeuratObject_WuetalBassezetal %>% NormalizeData %>% 
          FindVariableFeatures(selection.method="vst", nfeatures=TotalGeneNumb)  
          # nfeatures 20,000 is bettern than 5000
GeneSymb <- rownames(SeuratObject_WuetalBassezetal_Normal@assays$RNA@features)
SeuratObject_WuetalBassezetal_Scaled <- ScaleData(SeuratObject_WuetalBassezetal_Normal, features=GeneSymb) 
# This step requires a lot of memory. It may not run in your local computer. You will need a supercomputer.

## saveRDS(SeuratObject_WuetalBassezetal_Scaled, file="SeuratObjject_WuetalBassezetal_RNAScaled.rds")    
## SeuratObject_WuetalBassezetal_Scaled <- readRDS("SeuratObject_WuetalBassezetal_RNAScaled.rds")

#################################################################################################
## Step6b. calculate correlation between gene expression and cell type fraction in cancer epithelial cells
#################################################################################################

MetaData_WuetalBassezetal <- SeuratObject_WuetalBassezetal_Scaled@meta.data; dim(MetaData_WuetalBassezetal) 
MetaData_WuetalBassezetal[,which(colnames(MetaData_WuetalBassezetal)=="Cell_ID")] <- NULL   
## Delete "Cell_ID" column because it has some NAs.  "CellID" column also has some NAs. Fix it. 
MetaData_WuetalBassezetal$CellID <- rownames(MetaData_WuetalBassezetal)   
MetaData_WuetalBassezetal$CaseID <- MetaData_WuetalBassezetal$orig.ident
MetaData_WuetalBassezetal[1:2,]        # I have "CaseID" column already

##  inner_join between MetaData_WuetalBassezetal and ProportionCalculationMacroRatio
MetaData_WuetalBassezetal_CellTypeFraction <- dplyr::inner_join(MetaData_WuetalBassezetal,
                                 ProportionCalculationMacroRatio) %>%   # by "CaseID"
                                 dplyr::select(-c(orig.ident, nCount_RNA, nFeature_RNA, percent.mt, 
                                 nCount_SCT, nFeature_SCT, CellTypeMajor, CellTypeMinor,CellTypeSubset, 
                                 seurat_clusters,SCT_snn_res.1, CellID,  BC_type, cellType, 
                                 CellTypeByMarker_GSE176EGA6608, DatasetID)); 
dim(MetaData_WuetalBassezetal_CellTypeFraction) # 70850     18
table(MetaData_WuetalBassezetal_CellTypeFraction$CellTypeMacroTcell_GSE176EGA6608)

## calculate cell fraction for each cell type. 
## count the cells per cell type and per patients 
CellTypeCount <- table(MetaData_WuetalBassezetal_CellTypeFraction$CaseID, 
                       MetaData_WuetalBassezetal_CellTypeFraction$CellTypeMacroTcell_GSE176EGA6608) 

CellTypeCount_DF <- data.frame(matrix(CellTypeCount, nrow=nrow(CellTypeCount)))
rownames(CellTypeCount_DF) <- rownames(CellTypeCount); colnames(CellTypeCount_DF) <- colnames(CellTypeCount)
dim(CellTypeCount); CellTypeCount   # 14  10


CellTypeRatio <- CellTypeCount/rowSums(CellTypeCount); class(CellTypeRatio); rowSums(CellTypeRatio); 
class(CellTypeRatio)
## Make data.frame format. 
CellTypeRatio_CaseID <- data.frame(rbind(CellTypeRatio)) %>% tibble::rownames_to_column("CaseID")
colnames(CellTypeRatio_CaseID) <- gsub("\\.","_", colnames(CellTypeRatio_CaseID))

################################################################################
## Extract gene scaled expression data from seurat object 
################################################################################
ScaledExp <- SeuratObject_WuetalBassezetal_Scaled[["RNA"]]$scale.data; ScaledExp[1:2,1:5] %>% 
  data.frame; dim(ScaledExp) # 33186 90532  #  @scale.data is old version of Seurat package. 

#### Filter only protein coding genes because All Genes are too many, so I can't run for loop below ##
HumanGeneInfoFile <- "Homo_sapiens.gene_info.txt"  ## You can download this file from NCBI FTP
HumanGeneInfo <- fread(HumanGeneInfoFile, header=TRUE, stringsAsFactors=FALSE); dim(HumanGeneInfo);
HumanGeneInfo[1:3,1:10]

length(HumanGeneInfo$Symbol[HumanGeneInfo$type_of_gene=="protein-coding"]) # 19992
ScaledExp_AllGene <- ScaledExp[rownames(ScaledExp) %in% 
                                 (HumanGeneInfo$Symbol[HumanGeneInfo$type_of_gene=="protein-coding"]), ]; 
dim(ScaledExp_AllGene); ScaledExp_AllGene[1:3,1:3]  # 17284 71610

## ========= ##======== Filter out genes of low variance ====== ## =========== 
vars <- apply(ScaledExp_AllGene, 1, var); table(vars==0); table(vars<0.01) # FALSE: 16660 TRUE 988  
saveRDS(vars, "CalculateVariance_AcrossAllCell.rds")

ScaledExp_HighVarGene <- ScaledExp_AllGene[vars >= 0.01,]; dim(ScaledExp_HighVarGene) # 16660 90532
saveRDS(ScaledExp_HighVarGene, file="ScaledExp_HighVarGene_12942g.rds")
## ========= ######### ============== ############# =======  ######## ========

ScaledExp_Process <- ScaledExp_HighVarGene %>% t %>% data.frame 
saveRDS(ScaledExp_Process, file="ScaledExp_Process.rds")
ScaledExp_Process[1:3,1:3]  

## Let's remove genes that have low expression across 70,850 cells
ScaledExp_RmvLowExpGene <- ScaledExp_Process %>% tibble::rownames_to_column("CellID") %>% 
  dplyr::inner_join(MetaData_WuetalBassezetal[,c("CellID","CaseID","CellTypeMacroTcell_GSE176EGA6608")]) %>% 
  dplyr::select(-CellID)
table(ScaledExp_RmvLowExpGene$CaseID); dim(ScaledExp_RmvLowExpGene);   # 70850 16402
# saveRDS(ScaledExp_RmvLowExpGene, file="ScaledExp_RmvLowExpGene.rds")
ScaledExp_RmvLowExpGene <- readRDS(file="ScaledExp_RmvLowExpGene.rds")

```
 
Calculate mean of gene expression per caseID and calculate correlation 
 
```{Step6c}
#################################################################################################
## Step6c. Calculate mean of gene expression per caseID and calculate correlation 
## between gene expression and Macrophage fraction in different cell types. 
#################################################################################################
MyCellTypeSubset <- unique(ScaledExp_RmvLowExpGene$CellTypeMacroTcell_GSE176EGA6608); print(MyCellTypeSubset)

CorrTestSummary_AllCellType <- data.frame(matrix(ncol=1, nrow=ncol(ScaledExp_RmvLowExpGene))); CellTypeCount=0;
for(EachCellTypeSubset in MyCellTypeSubset) {
    # EachCellTypeSubset <- MyCellTypeSubset[6]
    CellTypeCount=CellTypeCount+1;
    print(paste0("currnt cell type: ", EachCellTypeSubset))
    ScaledExp_RmvLowExpGene_Subset <- ScaledExp_RmvLowExpGene %>% 
      dplyr::filter(CellTypeMacroTcell_GSE176EGA6608==EachCellTypeSubset) %>% 
      dplyr::select(-CellTypeMacroTcell_GSE176EGA6608); dim(ScaledExp_RmvLowExpGene_Subset) # 3014 245
    
    ## Calculate mean of all gene expression per case ID. 
    GeneExp_MeanInAllCase <- data.frame(); 
    for (EachCase in unique(ScaledExp_RmvLowExpGene_Subset$CaseID)) {
          # EachCase <- "CID4461"
          ScaledExp_AllGene_SubsetByCase <- ScaledExp_RmvLowExpGene_Subset %>% filter(CaseID==EachCase) %>% 
                    dplyr::select(-CaseID)  
          dim(ScaledExp_AllGene_SubsetByCase)  # 3006   264
          ScaledExp_MeanGeneExp <- data.frame(t(colMeans(ScaledExp_AllGene_SubsetByCase)))
          rownames(ScaledExp_MeanGeneExp) <- EachCase
          GeneExp_MeanInAllCase <- rbind(GeneExp_MeanInAllCase, ScaledExp_MeanGeneExp)
    }
    dim(GeneExp_MeanInAllCase); GeneExp_MeanInAllCase[1:8,1:5]; # 25 16400 

    GeneExp_MeanInAllCase_CaseID <- GeneExp_MeanInAllCase %>% tibble::rownames_to_column("CaseID"); 
    GeneExp_MeanInAllCase_CaseID[1:3,1:3]
    colnames(GeneExp_MeanInAllCase_CaseID) <- gsub("\\.","_",colnames(GeneExp_MeanInAllCase_CaseID))
    
    ################################################################################
    ## inner_join between ProportionCalculationMacroRatio, CellTypeRatio_CaseID,  
    ## and  GeneExp_MeanInAllCase_CaseID
    ################################################################################
    CellTypeRatio_GeneMeanExp<-dplyr::inner_join(ProportionCalculationMacroRatio[,c("CaseID","MacrophageRatio")], 
                      CellTypeRatio_CaseID) %>% dplyr::inner_join(GeneExp_MeanInAllCase_CaseID) %>% 
      data.frame; dim(CellTypeRatio_GeneMeanExp)
    colnames(CellTypeRatio_GeneMeanExp) <- gsub("\\.","_",colnames(CellTypeRatio_GeneMeanExp))
    CellTypeRatio_GeneMeanExp[1:2,1:18]

    #####################       #####################       #####################       #####################     
    ## == Calculate the correlation between Macrophage fraction and all gene expression in different cell types. 
    #####################       #####################       #####################       #####################    
    MyCellTypeFraction <- colnames(CellTypeRatio_CaseID)[-1]  
    MyGene <- colnames(GeneExp_MeanInAllCase_CaseID)[-1]; length(MyGene)  # 16660

    print(paste0("cell type count: ", CellTypeCount))
    CorrX <- CellTypeRatio_GeneMeanExp[, "Macrophage"];
    All_CorrTestSummary_Cyto <- data.frame(); SubLoopCount=0;
    for(EachGene in MyGene) {
          # EachGene <- "S100A11" # MyGene[706]
          print(paste0("MyGene: ", EachGene))
          SubLoopCount <- SubLoopCount + 1; 
          print(paste0("current subloopnumb: ", SubLoopCount))
          CorrY <- CellTypeRatio_GeneMeanExp[, EachGene];
          CorrTest <- cor.test(CorrX, CorrY, method=c( "spearman"), exact=FALSE)
          CorrTestSummary_Cyto <- data.frame(t(c(CorrTest$estimate, CorrTest$p.value)))
          colnames(CorrTestSummary_Cyto) <- c(paste0("Macrophage_Rho_In",EachCellTypeSubset),
                                              paste0("Macrophage_pval_In",EachCellTypeSubset))
          rownames(CorrTestSummary_Cyto) <- EachGene
          
          All_CorrTestSummary_Cyto <- rbind(All_CorrTestSummary_Cyto, CorrTestSummary_Cyto)
    }
    if(CellTypeCount==1) {
          CorrTestSummary_AllCellType <- All_CorrTestSummary_Cyto
    } else if (CellTypeCount>=1) {
          CorrTestSummary_AllCellType<-cbind(CorrTestSummary_AllCellType,  All_CorrTestSummary_Cyto)
    }
}

### When I get correlation file between Macrophage infiltration and all 
## gene exp in all cell type, several rows (genes) have NA for all Cell types.
dim(CorrTestSummary_AllCellType) # 16400 28
table(is.na(rowMeans(CorrTestSummary_AllCellType)))  # FALSE 16341  TRUE 461
CorrTestSummary_AllCellType[is.na(CorrTestSummary_AllCellType)] <- ""
fwrite(CorrTestSummary_AllCellType, file="SpearCorrel_MacrophageFraction_AllGeneExp.txt", 
       row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

```

Linedotplot of correlation - Figure 3A

```{Step6d}
#################################################################################################
## Step6d. Linedotplot of correlation - Figure 3A
#################################################################################################
colnames(CorrTestSummary_AllCellType)[1] <- "GeneSymb"; dim(CorrTestSummary_AllCellType); 
CorrTestSummary_AllCellType[1:2,]
AllCellType<-unique(gsub("(.*)_","",colnames(CorrTestSummary_AllCellType)[2:ncol(CorrTestSummary_AllCellType)]))

RhoThres=0.5; pvalThres=0.05
# for (EachCellType in AllCellType) {
    EachCellType <- AllCellType[6]; print(EachCellType) # InCancerEpithelial
    CellType_Rho <- paste0("Macrophage_Rho_",EachCellType); CellType_pval <- 
                paste0("Macrophage_pval_",EachCellType); print(CellType_Rho)
    CorrTestSummary_AllCellType_Macro <- CorrTestSummary_AllCellType %>% data.frame %>% 
        dplyr::select_if(colnames(CorrTestSummary_AllCellType) %in% c("GeneSymb", CellType_Rho,  CellType_pval)) 
    
    CorrTestSummary_AllCellType_Macro$SigGene <-ifelse(CorrTestSummary_AllCellType_Macro$GeneSymb=="S100A11"|
                               CorrTestSummary_AllCellType_Macro[, CellType_Rho] >= 0.60 | 
                               CorrTestSummary_AllCellType_Macro[, CellType_Rho] <= -0.65 ,
                               CorrTestSummary_AllCellType_Macro$GeneSymb, ""  )
    CorrTestSummary_AllCellType_Macro$PosNeg<-ifelse(CorrTestSummary_AllCellType_Macro[,CellType_Rho]>=RhoThres & 
                              (CorrTestSummary_AllCellType_Macro[, CellType_pval]) < pvalThres, "Pos", 
                              ifelse(CorrTestSummary_AllCellType_Macro[, CellType_Rho] <= -RhoThres & 
                              (CorrTestSummary_AllCellType_Macro[, CellType_pval]) < pvalThres, "Neg","Stable"))
    
    table(CorrTestSummary_AllCellType_Macro$PosNeg ); length(table(CorrTestSummary_AllCellType_Macro$SigGene))  
    ## Sort by Macrophage_Rho ascending order
    CorrTestSummary_AllCellType_Macro_RhoSort<-
            CorrTestSummary_AllCellType_Macro[with(CorrTestSummary_AllCellType_Macro,
            order(CorrTestSummary_AllCellType_Macro[, CellType_Rho], decreasing=FALSE)),] %>%
            dplyr::mutate(Index=seq(1:nrow(CorrTestSummary_AllCellType_Macro))) %>% data.frame
    ######################################################################
    ##  Make dot plot
    ######################################################################
    CorrelType <- length(table(CorrTestSummary_AllCellType_Macro$PosNeg )); print(CorrelType)
    if(CorrelType==3) {
          MyColorAnnot <- c("blue", "red", "black")
    } else if (CorrelType==2) {
          if( any(names(table(CorrTestSummary_AllCellType_Macro$PosNeg ))== c("Pos", "Stable")))  {
              MyColorAnnot <- c("red", "black")
          } else if ( any(names(table(CorrTestSummary_AllCellType_Macro$PosNeg ))== c("Neg", "Stable"))) {
              MyColorAnnot <- c("blue", "black")
          }
    }
    
    colnames(CorrTestSummary_AllCellType_Macro_RhoSort)[2] <- "CorrelRho"
    Dotplot<-ggplot(CorrTestSummary_AllCellType_Macro_RhoSort,aes(Index,CorrelRho,colour=PosNeg,label=SigGene))+ 
      geom_point() + scale_color_manual(values=MyColorAnnot)+  ylim(c(-1.0, 1.0)) +  # "blue",
      theme(legend.position="none", panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black")) # theme_bw(base_size = 10)
    # Dotplot
    MyPlot<-Dotplot + geom_text_repel(box.padding=0.5, size=2,  max.overlaps = Inf)
    MyPlot
    
    OutFile <-paste0("Fig3A_OutLineDotplot_CorrelBtwGeneExpMacroFract_InCanEpi.pdf");print(OutFile)
    ggsave(MyPlot, file=OutFile, width=8,height=8)
# }  # end of for loop  

```

Outgrid correlation dotplot - Figure 3E  

```{Step6e}
    
#################################################################################################
## Step6e. Outgrid correlation dotplot - Figure 3E  
##          Correlation between Macrophage infiltration and S100A11 gene expression in different cell types.  
#################################################################################################
library(gridExtra)  # for 'grid.arrange' function      
      
# MyCellTypeSubset <- unique(SeuratObject_WuetalBassezetal@active.ident); 
# MyCellTypeSubset <- MyCellTypeSubset[!is.na(MyCellTypeSubset)]; print(MyCellTypeSubset)
MyInterestGene <- "S100A11"; CellTypeNumb<-0; MyPlot_AllCellType <- list(); 

CellTypeCount<- 0;
for(EachCellTypeSubset in MyCellTypeSubset) {
      # EachCellTypeSubset <- MyCellTypeSubset[6]
      CellTypeCount=CellTypeCount+1;
      print(paste0("currnt cell type: ", EachCellTypeSubset))
      ScaledExp_RmvLowExpGene_Subset <- ScaledExp_RmvLowExpGene %>% 
        dplyr::filter(CellTypeMacroTcell_GSE176EGA6608==EachCellTypeSubset) %>% 
        dplyr::select(-CellTypeMacroTcell_GSE176EGA6608); dim(ScaledExp_RmvLowExpGene_Subset)  # 3014 245
      
      ## Calculate mean of all gene expression per case ID. 
      GeneExp_MeanInAllCase <- data.frame(); 
      for (EachCase in unique(ScaledExp_RmvLowExpGene_Subset$CaseID)) {
            # EachCase <- "CID4461"
            ScaledExp_AllGene_SubsetByCase <- ScaledExp_RmvLowExpGene_Subset %>% filter(CaseID==EachCase) %>% 
              dplyr::select(-CaseID)  
            dim(ScaledExp_AllGene_SubsetByCase)  # 3006   264
            ScaledExp_MeanGeneExp <- data.frame(t(colMeans(ScaledExp_AllGene_SubsetByCase)))
            rownames(ScaledExp_MeanGeneExp) <- EachCase
            GeneExp_MeanInAllCase <- rbind(GeneExp_MeanInAllCase, ScaledExp_MeanGeneExp)
      }
      dim(GeneExp_MeanInAllCase); GeneExp_MeanInAllCase[1:8,1:5]; # 25 16400 
      #             TNFRSF18    TNFRSF4    TNFRSF14   TNFRSF25    TNFRSF9
      # CID3586  -0.2614056 -0.2704576  0.18663825 -0.2035496 -0.1529137
      # CID4066  -0.1950509 -0.2352626 -0.16339175 -0.2150040 -0.1621458
      GeneExp_MeanInAllCase_CaseID <- GeneExp_MeanInAllCase %>% tibble::rownames_to_column("CaseID"); 
      GeneExp_MeanInAllCase_CaseID[1:3,1:3]
      colnames(GeneExp_MeanInAllCase_CaseID) <- gsub("\\.","_",colnames(GeneExp_MeanInAllCase_CaseID))
      
      ################################################################################
      ## inner_join between ProportionCalculationMacroRatio  and  
      ## GeneExp_MeanInAllCase_CaseID, and select only S100A11 exp column
      ################################################################################
      CellTypeRatio_GeneMeanExp <- dplyr::inner_join(ProportionCalculationMacroRatio, 
                         GeneExp_MeanInAllCase_CaseID[,c("CaseID",MyInterestGene)]) %>% 
                         data.frame; dim(CellTypeRatio_GeneMeanExp) # 25 16402
      CellTypeRatio_GeneMeanExp[1:2, ]
      CorrX <- CellTypeRatio_GeneMeanExp[, "S100A11"]; 
      CorrY <- CellTypeRatio_GeneMeanExp[, "Macrophage"]; 
      
      ## Corr test
      CorrTest <- cor.test(CorrX, CorrY, method=c("spearman"), exact=FALSE); print(CorrTest)
      # plot(CorrX, CorrY)
      
      if(is.na(CorrTest$estimate)) {
        CorrTest$estimate <- 0
        CorrTest$p.value <- 1
      }
      
      IPcol <- c("blue", "#66B2FF","#CCFFFF","#E0E0E0","#FFFFCC", "#feb24c", "#fd8d3c","#FF0000",
                 "red", "#CC0000")#, "#bd0026")
      
      XMin<- floor(min(CellTypeRatio_GeneMeanExp$S100A11)) ; print(paste0("XMin: ", XMin)) # -1
      XMax<- max(CellTypeRatio_GeneMeanExp$S100A11); print(paste0("XMax: ", XMax)) # 
      
      YMin<- floor(min(CellTypeRatio_GeneMeanExp$Macrophage)) ; print(paste0("YMin: ", XMin)) # -1
      YMax<- ceiling(max(CellTypeRatio_GeneMeanExp$Macrophage)); print(paste0("YMax: ", YMax)) # -1
      
      MidMac <-  mean(CellTypeRatio_GeneMeanExp$Macrophage); print(MidMac) # 3.284987
      
      if(EachCellTypeSubset == "Macrophage") {
          MyPlotColor <- ggplot(CellTypeRatio_GeneMeanExp, aes(S100A11, Macrophage)) + 
            geom_point(aes(color = Macrophage), alpha=1, size=5, show.legend=T)      
      
      } else if(EachCellTypeSubset != "Macrophage") {
          colnames(CellTypeRatio_GeneMeanExp)[colnames(CellTypeRatio_GeneMeanExp)==EachCellTypeSubset]<-
                  "CellTypeSelect"
          MyPlotColor <- ggplot(CellTypeRatio_GeneMeanExp, aes(S100A11, Macrophage)) + 
            geom_point(aes(color = CellTypeSelect), alpha=1, size=5, show.legend=T)      
      }
      
      
      MyPlot <- MyPlotColor + scale_colour_gradientn(colours = IPcol) + xlim(c(-0.9, 1.6)) + 
        ylim(c(-0.5, YMax+1)) +labs(x=paste0("S100A11 expression (scaled)"), y=paste0("TAM infiltration"),
             title=paste0(EachCellTypeSubset, " Sp Rho=", print(round(CorrTest$estimate,3)), 
                          ", pval=", print(round(CorrTest$p.value,5)) ) ) + theme_bw() +  
        theme(axis.text.x=element_text(colour="Black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),  
              axis.text.y=element_text(colour="Black",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
              plot.title = element_text(size=14),
              axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
              axis.title.y=element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),
              axis.line.y=element_line(color="black", linewidth=0.7),
              axis.line.x=element_line(color="black", linewidth =0.7),
              axis.ticks=element_line(colour="black",linewidth=1),
              axis.ticks.length=unit(.22, "cm"), text=element_text(size=22),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        stat_smooth(method = "glm", formula = y ~ x, geom = "smooth", color="grey", se=FALSE)
      MyPlot
      
      MyPlot_AllCellType <- c(MyPlot_AllCellType, list(MyPlot))

}      
      
length(MyPlot_AllCellType)  # 14

OutBoxplot <- paste0("Fig3E_OutGridDotplot_CorrelS100A11TAMInfiltration.pdf"); print(OutBoxplot)
pdf(OutBoxplot, width=30, height=20)
grid.arrange(MyPlot_AllCellType[[1]], MyPlot_AllCellType[[2]], MyPlot_AllCellType[[3]],MyPlot_AllCellType[[4]], 
         MyPlot_AllCellType[[5]], MyPlot_AllCellType[[6]], MyPlot_AllCellType[[7]],MyPlot_AllCellType[[8]], 
         MyPlot_AllCellType[[9]], MyPlot_AllCellType[[10]], MyPlot_AllCellType[[11]],MyPlot_AllCellType[[12]], 
         MyPlot_AllCellType[[13]],MyPlot_AllCellType[[14]], nrow=4, top="" )
dev.off()   


```
    
## Section 7.  S100A11 expression in different samples or in different cell types. 

```{S100A11}
#################################################################################################
## Step7a. Violin plot of S100A11 gene expression in different CaseID, when the samples
## are sorted by macrophage infiltration ascending order. 
#################################################################################################
ProfortionCal_SortByMacrophage <- ProportionCalculationMacroRatio[with(ProportionCalculationMacroRatio, 
                                                                       order(Macrophage, decreasing=FALSE)),]
SeuratObject_WuetalBassezetal_Scaled$CaseID <- factor(SeuratObject_WuetalBassezetal_Scaled$CaseID, 
                                                      levels=ProfortionCal_SortByMacrophage$CaseID)
MyViolinplot_ByCaseID <- VlnPlot(object = SeuratObject_WuetalBassezetal_Scaled, 
                                 features = c("S100A11"), group.by="CaseID" ) 
ggsave(MyViolinplot_ByCaseID, height=8,width=12, dpi=300, 
       filename="Fig3C_OutVlnPlot_S100A11Exp_InDifferentCaseID.pdf", useDingbats=FALSE)

#################################################################################################
## Step7b. Violin plot of S100A11 gene expression in different cell types. 
#################################################################################################
SeuratObject_WuetalBassezetal_Scaled$CellTypeMacroTcell_GSE176EGA6608 <- 
                            factor(SeuratObject_WuetalBassezetal_Scaled$CellTypeMacroTcell_GSE176EGA6608, 
                            levels=c("Bcell","NKcell","NKTcell","CD4Tcell","CD8Tcell","Treg","DC","Mastcell",
                                  "Macrophage","Monocyte","CancerEpithelial","Endothelial","Fibroblast","Myofibroblast"))
MyViolinplot_ByCellType <- VlnPlot(object = SeuratObject_WuetalBassezetal_Scaled, features = c("S100A11"), 
                                   group.by="CellTypeMacroTcell_GSE176EGA6608" ) 
ggsave(MyViolinplot_ByCellType, height=8,width=12, dpi=300,
       filename="Fig3D_OutVlnPlot_S100A11Exp_InDifferentCellType.pdf", useDingbats=FALSE)
```


**The End**
