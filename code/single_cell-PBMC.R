setwd("/PATH/FILE")

library(dplyr)
library(Seurat)
library(ggplot2)
library(MAST)

#========Read10X matrix, save rds file==========
data_dir <- './raw_mtx'
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data <- Read10X(data.dir = data_dir)
meta  <- read.csv("PBMC_meta.csv",row.names = "cell",check.names=FALSE)
sc_raw <- CreateSeuratObject(counts=data, meta.data=meta)
Idents(sc_raw) <- "source_name"
sc<-subset(sc_raw,idents=c("nCoV_PBMC","Normal_PBMC"))

#==========reload rds file, scRNA-seq analysis==========
sc <- PercentageFeatureSet(sc, pattern = "^MT-",col.name = "percent.mito")
sc.list <- SplitObject(sc, split.by = "pat_id")
for (i in 1:length(sc.list)) { print(i)
  sc.list[[i]] <- NormalizeData(sc.list[[i]], verbose=FALSE)
  sc.list[[i]] <- FindVariableFeatures(sc.list[[i]], selection.method="vst", nfeatures=2000, verbose=FALSE)
}
sc.anchors <- FindIntegrationAnchors(object.list=sc.list, dims=1:30)
sc.integrated <- IntegrateData(anchorset=sc.anchors, dims=1:30)
sc.integrated <- ScaleData(sc.integrated)

#===========PCA & tsne===========
sc.integrated <- RunPCA(sc.integrated, npcs=15, verbose=FALSE,seed.use = 42)
sc.integrated <- RunTSNE(sc.integrated, dims=1:15, verbose=FALSE,k.seed=1)
#Find Clusters
sc.integrated <- FindNeighbors(sc.integrated, dims=1:15, verbose=FALSE)
sc.integrated <- FindClusters(sc.integrated, resolution = 0.5, verbose=FALSE)

#=========DefaultAssay with all genes===========
# <*> #
DefaultAssay(sc.integrated) <- "RNA"
saveRDS(sc.integrated, file = "./PBMC_noFlu_DE.rds")
# <*> #

# ==== Markers from paper ====
markers <- c('CD3E', 'CD4', 'CCR7', 'CD8A',
             'NCAM1', 'CD14', 'FCGR3A', 'NR4A1',
             'CD19', 'FCER1A', 'PPBP', 'HBB')

pdf("./image/Markers_FeaturePlot_p.pdf",width=20, height=15)
print(FeaturePlot(sc.integrated, features = markers)+ NoLegend()+ NoAxes())
dev.off()

png("./image/Markers_FeaturePlot_p.png",width=32, height=20,units="cm",res=800)
plot(FeaturePlot(sc.integrated, features = markers)+ NoLegend())
dev.off()

#========tSNE plot======
cols=c("#AC92EC","#8067B7",
       "#4FC1E9","#4A89DC", "#F08080","#DA4453",
       "#E6E9ED","#CCD1D9","#AAB2BD",
       "#FC6E51",
       "#48CFAD","#6B8E23","#8CC152",
       "#FFCE54",
       "#BAA286",
       '#EC87C0' )

pdf("./image/tSNE_map_cell_legend.pdf", height=6, width=8)
print(DimPlot(sc.integrated, label=F, pt.size=0.5, label.size=6,group.by="celltype_sub",cols=cols)+ NoAxes())
dev.off()

png("./image/tSNE_map_cell_legend.png",width=16, height=12,units="cm",res=800)
print(DimPlot(sc.integrated, label=F, pt.size=0.5, label.size=6,group.by="celltype_sub",cols=cols)+ NoAxes())
dev.off()

