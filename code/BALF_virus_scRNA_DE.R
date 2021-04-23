setwd("/hdd/yuan/project/CELL_zhang_Virus")
library(dplyr)
library(Seurat)
library(MAST)

sc <- readRDS(file.path("./","./Virus_plus_BALF.rds"))
Idents(sc) <- "cellType"
print(table(Idents(sc)))
cell_list=c('Macrophage',
            'Plasma',
            'Secretory',
            'Neutrophil',
            'NK',
            'Ciliated',
            'T',
            'Squamous')

# ==== sex DE critical====
for (c in cell_list) {
  print(c)
  #m vs f
  markers <- FindMarkers(sc, ident.1="M", ident.2="F", group.by="Sex", subset.ident=c, logfc.threshold=0, min.pct=0.01,test.use = "MAST")
  write.table(markers, file=sprintf('./BALF_virus_DE/BALF_virus_DE_%s.tsv', gsub("/"," ",c)), quote=FALSE, sep="\t")
  print(dim(markers))
}


