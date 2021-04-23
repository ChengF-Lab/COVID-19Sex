library(Seurat)
library(edgeR)

setwd("/PATH/FILE")
#'covid_nbt_main.rds' download from https://doi.org/10.6084/m9.figshare.12436517.
#Save sub-files of COVID-19 severity
sc <- readRDS(file.path("./","covid_nbt_main.rds"))
severity_list=c("critical","moderate","control")
Idents(sc) <- "severity"
print(table(Idents(sc)))
for (s in severity_list) {
  print(s)
  sub_severity<-subset(sc,idents=c(s))
  saveRDS(sub_severity, file=sprintf('./NS_%s.rds', gsub("/"," ",s)))
}
  
#==========Reload rds. files for each COVID-19 severity=========
#cell_list may vary because cell number was lower than 10 in some celltypes.
cell_list=c('Basal',
            'Ciliated',
            'Ciliated-diff',
            'FOXN4',
            'Ionocyte',
            'IRC',
            'Secretory',
            'Secretory-diff',
            'Squamous',
            'B cell',
            'CTL',
            'MC',
            'MoD-Ma',
            'moDC',
            'Neu',
            'NK',
            'NKT',
            'NKT-p',
            'nrMa',
            'pDC',
            'rMa',
            'Treg')

sub_severity <- readRDS(file.path("./","NS_critical.rds"))
Idents(sub_severity) <- "celltype"
print(table(Idents(sub_severity)))
for (c in cell_list) {
  print(c)
  sub<-subset(sub_severity,idents=c(c))
  #================DE==========
  count_ <- sub$RNA@counts
  meta_ <- sub@meta.data
  y <- DGEList(counts=count_, group=meta_$sex)
  y <- calcNormFactors(y)
  design <- model.matrix(~ meta_$sex)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit)
  print(qlf$comparison)
  write.table(topTags(qlf, n=Inf), file=sprintf('./NS_DE/critical_%s_DE_inf.tsv',c), quote=FALSE, sep="\t")
  DE_file<-topTags(qlf,n=Inf)
  DEGs_q05_FC_0 <- DE_file$table$FDR < 0.05 & abs(DE_file$table$logFC) > 0
  DEGs_q05_FC_0<-DE_file[DEGs_q05_FC_0,]
  write.table(DEGs_q05_FC_0, file=sprintf('./NS_DE/critical_%s_DE_05.tsv',c), quote=FALSE, sep="\t")
}


