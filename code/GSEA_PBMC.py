# -*- coding: utf-8 -*-
import os
import pandas as pd 
import numpy as np 
import re
import gseapy as gp
import matplotlib.pyplot as plt
#PBMC_DE generate from PBMC_scRNA_DE.R
#Keep immunne_pathway_kegg22.gmt in folder PBMC_DE; make a folder in PBMC_DE, named GSEA_gender
os.chdir("/PATH/FILE")
listname = os.listdir("/PATH/FILE")
#print(listname)

for name in listname:
    if ".tsv" in name:
        f = open (name)
        a=f.readline()
        f_w = open ("name_"+name,"w")
        f_w.write("gene\t"+a)
        f_w.write(f.read())
        f.close()
        f_w.close()

listname_new = os.listdir("/PATH/FILE")
#print(listname_new)
for f in listname_new:
    if "name_" in f:
        print(f)
        df_DE=pd.read_csv(f,sep="\t")
        if df_DE.shape[0]>50:
            print(df_DE.shape)
            df_DE=df_DE.sort_values(by='logFC',ascending=False)
            print(df_DE.head())
            gene = df_DE["gene"].tolist()
            rank_value = df_DE['logFC'].tolist()
            gene_value_list=[]
            for i in zip(gene,rank_value):
                gene_value_list.append(list(i))
            df_prerank=pd.DataFrame(gene_value_list)
            print(df_prerank)
            print('===========')
            #GSEA-pre-rank 
            #note: multiprocessing may not work on windows
            pre_res = gp.prerank(rnk=df_prerank, gene_sets='./immunne_pathway_kegg22.gmt',
                                 processes=4,
                                 min_size=1,
                                 permutation_num=1000, # reduce number to speed up testing
                                 outdir='./GSEA_png/prerank_zscore_report_GSEA'+str(name), format='png')
            
            print(pre_res.res2d.head())
            final_name=f.strip(".tsv")
            pre_res.res2d.to_csv("./GSEA_gender/GSEA-%s.csv" % final_name)

