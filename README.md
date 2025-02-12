# DESeq2-RNASeq-Pipeline
This repository contains an R pipeline for differential gene expression analysis using DESeq2. It processes gene expression quantification results generated after the Dragen RNA-seq analysis pipeline.

## Requirements

- **R version**: 4.4.2 (2024-10-31 ucrt)
- **RStudio version**: 2024.12.0 Build 467
- **Required R packages**:
  - `DESeq2` (version 1.46.0)
  - `tximport` (version 1.34.0)
  - `dplyr` (version 1.1.4)
  - `rtracklayer` (version 1.66.0)

## Pipeline Overview
The pipeline performs the following steps:
- Import transcript quantification data from `.quant.sf` files (Salmon output).
- Create a `DESeqDataSet` and filter low-expression genes.
- Perform differential expression analysis (DEG) using DESeq2.
- Generate PCA plots to visualize sample clustering.
- Annotate genes using a reference GTF file.

## Input & Output Files
- .quant.sf files – transcript quantification results from Salmon
- .gtf file – gene annotations for transcript-to-gene mapping
- Differential expression results (.csv) – contains log2 fold changes and p-values
- Annotated results (_annotated.csv) – includes gene names from GTF
- Filtered mitochondrial genes (_mt.csv)
- PCA plots for visualization

## License & Citation
If you use this pipeline in your research, please cite the appropriate sources.
