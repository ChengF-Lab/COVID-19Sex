library(dplyr)
library(Seurat)
setwd("/hdd/yuan/project/2_covid19/PBMC_check_0415")
#PBMC single cell samples
sc <- readRDS(file.path("./","PBMC_noFlu_DE.rds"))
severity_list=c('Normal','COV_Mild','COV_Severe')
Idents(sc) <- "Disease_type"
print(table(Idents(sc)))
for (c in severity_list){
  print(c)
  sub<-subset(sc,idents=c)
  
  Idents(sub) <- "gender"
  print(table(Idents(sub)))
  
  sub_m<-subset(sub,idents=c("male"))
  # Take raw data and normalize it.q90
  count_raw <- sub_m$RNA@counts
  count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
  write.table(count_norm,file=sprintf('./CPDB/%s_male_cpdb_count.txt',c), sep='\t', quote=F)
  # Generating metadata file.
  meta_data <- cbind(rownames(sub_m@meta.data), sub_m@meta.data[,'celltype_sub', drop=F])
  write.table(meta_data, file=sprintf('./CPDB/%s_male_cpdb_meta.txt',c), sep='\t', quote=F,row.names=F)
  
  sub_f<-subset(sub,idents=c("female"))
  # Take raw data and normalize it.
  count_raw_f <- sub_f$RNA@counts
  count_norm_f <- apply(count_raw_f, 2, function(x) (x/sum(x))*10000)
  write.table(count_norm_f,file=sprintf('./CPDB/%s_female_cpdb_count.txt',c), sep='\t', quote=F)
  # Generating metadata file.
  meta_data_f <- cbind(rownames(sub_f@meta.data), sub_f@meta.data[,'celltype_sub', drop=F])
  write.table(meta_data_f,file=sprintf('./CPDB/%s_female_cpdb_meta.txt',c), sep='\t', quote=F,row.names=F)
}