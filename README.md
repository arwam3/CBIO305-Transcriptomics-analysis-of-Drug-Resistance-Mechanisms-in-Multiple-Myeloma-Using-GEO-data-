# 🔬CBIO305-Transcriptomics-analysis-of-Drug-Resistance-Mechanisms-in-Multiple-Myeloma-Using-GEO-data-
This project accompanies the study:

**📄 Citation** 

Riz I, Hawley TS, Hawley RG. KLF4-SQSTM1/p62-associated prosurvival autophagy contributes to carfilzomib resistance in multiple myeloma models. Oncotarget 2015 Jun 20;6(17):14814-31. PMID: 26109433
## 📚 Table of Contents

- [Overview](#-overview)
- [Workflow](#-Workflow)
- [Repository Structure](#-repository-structure)
- [Data Sources](#-data-sources)
- [Setup and Requirements](#️-setup-and-requirements)


## 🧪 Overview

This project investigates transcriptomic alterations underlying **carfilzomib resistance** in **multiple myeloma (MM)** using microarray data from the GEO database (**GSE69078**). The analysis focuses on identifying **differentially expressed genes (DEGs)** and associated **autophagy-related survival mechanisms**, including those regulated by **KLF4** and **SQSTM1/p62**.

![Graphical Abstract](https://github.com/user-attachments/assets/b2784e7e-89b4-4c88-ac17-7646b328f63b)

### 🔬Workflow

1. **Microarray preprocessing**: QC, normalization, and expression matrix construction using Bioconductor.
2. **DEA**: Differential expression using the `limma` package.
3. **Enrichment**: GO/KEGG pathway analysis and PPI network construction.
4. **Visualization**: Volcano plots, UMAP, PCA, and heatmaps


## 📁 Repository Structure

.
├── data/
│ ├── raw/ # Raw CEL files (GSE69078)
│ └── processed/ # Normalized matrix, DEG lists
├── R scripts/
│ ├── DEA and visualization.R # Differential expression analysis and plots
│ └── Enrichment Analysis.R # GO/KEGG/PPI enrichment
├── figures/ # Plots: volcano, PCA, UMAP, etc.
├── 
---


## 📁 Repository Structure

.
├── Workflow_diagram.png # Pipeline diagram (optional)
├── Data/
│ ├── raw/ # Raw CEL files (GSE69078)
│ └── processed/ # Normalized matrix, DEG lists
├── scripts/
│ ├── DEA and Visualization.R # Differential expression
│ └── Enrichment analysis.R # GO/KEGG/PPI enrichment
├── Figures/ # Plots: volcano, PCA, UMAP, etc.
├── README.md # Project documentation
└──Rrequirements.txt # R package dependencies

---

## 📂 Data Sources

### 🔸 Microarray Expression Data
- **Source**: GEO dataset [GSE69078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69078)
- **Platform**: Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)
- **Samples**: KMS-11 and KMS-34 MM cell lines (parental and carfilzomib-resistant variants)

---

## ⚙️ Setup and Requirements

Install R and required packages:

```r
install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery",
  "limma",
  "clusterProfiler",
  "pheatmap",
  "plotly",
  "openxlsx",
  "umap",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "enrichplot"
))
