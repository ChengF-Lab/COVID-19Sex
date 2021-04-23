library(Seurat)
library(edgeR)

setwd("/hdd/yuan/project/2_covid19/PBMC_check_0415")
#load .rds file, a output file from single_cell-PBMC.R
sc <- readRDS(file.path("./","PBMC_noFlu_DE.rds"))
severity_list=c('Normal','COV_Mild','COV_Severe')
Idents(sc) <- "Disease_type"
print(table(Idents(sc)))
for (s in severity_list) {
  print(s)
  sub_severity<-subset(sc,idents=c(s))
  saveRDS(sub_severity, file=sprintf('./PBMC_%s.rds', gsub("/"," ",s)))
}

#==========Reload rds. files for each COVID-19 severity=========
#cell_list may vary because cell number was lower than 10 in some celltypes.
#PBMC_DE, make a folder before running code.
cell_list=c('B cell lgG-',
            'B cell lgG+',
            'CD4-T cell EM',
            'CD4-T cell nEM',
            'CD8-T cell EM',
            'CD8-T cell nEM',
            'DC',
            'Monocyte',
            'Monocyte-i',
            'Monocyte-n',
            'NK',
            'Platelet',
            'RBC')
sub_severity <- readRDS(file.path("./","PBMC_Normal.rds"))
Idents(sub_severity) <- "celltype_sub"
print(table(Idents(sub_severity)))
for (c in cell_list) {
  print(c)
  sub<-subset(sub_severity,idents=c(c))
  #===========DE==========
  count_ <- sub$RNA@counts
  meta_ <- sub@meta.data
  y <- DGEList(counts=count_, group=meta_$gender)
  y <- calcNormFactors(y)
  design <- model.matrix(~ meta_$gender)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit)
  print(qlf$comparison)
  write.table(topTags(qlf, n=Inf), file=sprintf('./PBMC_DE/Normal_%s_DE_inf.tsv',c), quote=FALSE, sep="\t")
  DE_file<-topTags(qlf,n=Inf)
  DEGs_q05_FC_0 <- DE_file$table$FDR < 0.05 & abs(DE_file$table$logFC) > 0
  DEGs_q05_FC_0<-DE_file[DEGs_q05_FC_0,]
  write.table(DEGs_q05_FC_0, file=sprintf('./PBMC_DE/Normal_%s_DE_05.tsv',c), quote=FALSE, sep="\t")
    
}
