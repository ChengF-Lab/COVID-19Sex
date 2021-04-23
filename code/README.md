# COVID-19Sex

paper "Multimodal Single-Cell Omics Analysis Identifies Epithelium-Immune Cell Interactions and Immune Vulnerability Associated with Sex Differences in COVID-19"

Data set
* `data/PBMC_meta.csv` meta informations of PBMC single cell data.
* `data/immunne_pathway_kegg22.gmt` is the immune pathways for GSEA analysis.
* .rds files are available at https://figshare.com/s/1e9dc06d2b80b7361f99

Code
* `code/single_cell-PBMC.R`, is single cell analysis code for raw matrix.
* `code/PBMC_scRNA_DE.R`, code of differential expression analysis for each cell type of PBMCs by sex
* `code/CPDB_PBMC_sex.R`, generate CellPhoneDB input files from Seurat.
* `code/cellphoneDB_PBMC.sh`, CellPhoneDB analysis, more,details please refer to 2.	Efremova, M. et al. Nat Protoc. 2020.
* `code/GSEA_PBMC.py`, GSEA analysis.
* `code/KM_PS_analysis.R`, code cumulative hazard analysis based on propensity score match method.
* `code/NS_scRNA-DE.R`, code of differential expression analysis for each cell type of nasal samples by sex
* `code/BALF_virus_scRNA_DE.R`, code of differential expression analysis for each cell of BALF samples type by sex

Requirements
* Python 3.7. The required dependencies for pandas, numpy, scipy, scikit-learn, and gseapy.
* R 4.0.3. The required dependencies for ggplot2, dplyr, MatchIt,Matching,survival, survminer, Seurat, edgeR and MAST.

Reference:
Efremova, M. et al. CellPhoneDB: inferring cell-cell communication from combined expression of multi-subunit ligand-receptor complexes. Nat Protoc. 2020 Apr;15(4):1484-1506.

