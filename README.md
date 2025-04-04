# Chadarevian_2024_Neprilysin
Seurat analysis files associated with the single cell RNA-seq data from Chadarevian et al. (2025) Harnessing human iPSC-microglia for CNS-wide delivery of disease-modifying proteins

This analysis utilizes the following datasets:
  - Microglia promote anti-tumour immunity and suppress breast cancer brain metastasis 
  - Evans et al. 2023 PMID: 37957324
  - GEO: GSE147948
    - Samples:
      - Control
      - Met1
      - Met2

  - Development of a Chimeric Model to Study and Manipulate Human Microglia In Vivo 
  - Hasselmann et al. 2019 PMID: 31375314
  - GEO: GSE133433
  - Samples:
    - WT_Male
    - WT_Female
    - 5X_Male
    - 5X_Female

  - Cuprizone Demyelination
  - This Paper
  - GEO:
  - Samples:
    - Sublibrary1
    - Sublibrary2

In this publication, the Hasselmann et al. and Evans et al. 10X Genomics datasets were aligned to the Ensembl GRCh38 release 111 transcriptome, to retain only protein coding genes in the analysis, using CellRanger v6.1.2. The Cuprizone demyelination Parse dataset was aligned to the Ensembl GRCh38 release 111 transcriptome using Parse's split-pipe v1.2.0. Subsequent analysis for all datasets was performed in Seurat v5.0.1 using R v4.3.3.
  -   sctransform
  -   dplyr
  -   cowplot
  -   ggplot2
  -   jsonlite

The analysis steps that were taken include:
  - Dataset quality control was performed on each individual sample using the process outlined in "2024_Chadarevian_NEP_Seurat_v5_QC.R" and the analysis parameters included in "XXXXXXXXX"
  - Sample merging and integration within each dataset using the relevant R scripts:
      - "2024_Chadarevian_NEP_5XfAD_Seurat_v5_Integrated.R" for the Hasselmann et al. dataset
      - "2024_Chadarevian_NEP_Cancer_Seurat_v5_Integrated.R" for the Evans et al. dataset
      - "2024_Chadarevian_NEP_Cuprizone_Seurat_v5_Integrated.R" for the Cuprizone demyelination dataset from this paper
