# :dna: CBIO305 - Transcriptomics Analysis of Drug Resistance Mechanisms in Multiple Myeloma Using GEO Data

This project accompanies the study:

**📄 Citation**  
Riz I, Hawley TS, Hawley RG. KLF4-SQSTM1/p62-associated prosurvival autophagy contributes to carfilzomib resistance in multiple myeloma models. Oncotarget 2015 Jun 20;6(17):14814-31. PMID: 26109433

---

## 📚 Table of Contents

- [Overview](#-overview)
- [Workflow](#-workflow)
- [Repository Structure](#-repository-structure)
- [Data Sources](#-data-sources)
- [Setup and Requirements](#-setup-and-requirements)
- [Visulaization and interpretation](#visualization_and_interpretation)
- [Top 5 Differentially Expressed Genes (DEGs)](#top_5_differentially_expressed_genes_)

---

## :pushpin: Overview

This project investigates transcriptomic alterations underlying **carfilzomib resistance** in **multiple myeloma (MM)** using microarray data from the GEO database (**GSE69078**). The analysis focuses on identifying **differentially expressed genes (DEGs)** and associated **autophagy-related survival mechanisms**, including those regulated by **KLF4** and **SQSTM1/p62**.

![Graphical Abstract](https://github.com/user-attachments/assets/b2784e7e-89b4-4c88-ac17-7646b328f63b)  


---

## :clipboard: Workflow

1. **Microarray preprocessing:** QC, normalization, and expression matrix construction using Bioconductor.  
2. **DEA:** Differential expression using the `limma` package.  
3. **Enrichment:** GO/KEGG pathway analysis and PPI network construction.  
4. **Visualization:** Volcano plots, UMAP, PCA, and heatmaps.

---

## 📁 Repository Structure

```
.

├── R Script/
│   ├── DEA_and_Visualization.R       # Differential expression and visualization
│   └── Enrichment_Analysis.R         # GO/KEGG/PPI enrichment
├── Figures/                          # Plots: volcano, PCA, UMAP, etc.
    ├── Figure 2.png
    ├── Figure 3.png
    ├── Figure 4.png
    └── Figure 5.png
├──  Significant_DEGs_filtered.xlxs   # Normalized matrix, DEG lists
└── README.md                         # Project documentation

```

---

## 📂 Data Sources

### 🔸 Microarray Expression Data  
- **Source:** GEO dataset [GSE69078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69078)  
- **Platform:** Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)  
- **Samples:** KMS-11 and KMS-34 MM cell lines (parental and carfilzomib-resistant variants)  

---

## ⚙️ Setup and Requirements

Install R and the required packages with:

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
```

---
---

## :bar_chart: Visualization and Interpretation

### Figure 2: Box Plot & Heatmap
- **A.** The box plot shows that the overall distribution of gene expression levels (log2 transformed) is remarkably similar across all four sample groups: KMS 11 control, KMS 11 resistant, KMS 34 control, and KMS 34 resistant. The medians, interquartile ranges, and overall spread of the data points within each group appear consistent.
- **B.**  The heatmap visually represents the expression levels of genes identified as differentially expressed (DEGs).
- **C.** he histogram shows the distribution of adjusted p-values from the differential gene expression analysis. A "good" p-value distribution for a differential expression analysis typically shows a peak near 0 and a relatively flat distribution for higher p-values (approaching 1)

### Figure 3 : MD Plot

- **A. KMS-11 Resistant vs Control:**  
Most genes show minimal change, but some are significantly up- or downregulated (logFC > ±1), indicating selective gene expression changes associated with resistance.

- **B. KMS-34 Resistant vs Control:**  
Stronger transcriptional shifts are observed, with many genes showing logFC > ±2, suggesting a more robust resistance-associated response.

- **C. KMS-11 vs KMS-34 Resistant:**  
Major expression differences (logFC > ±4) reveal distinct resistance mechanisms between the two cell lines.

---
 
### Figure 4 : Volcano Plot

- **A. KMS-11 Resistant vs Control:**  
Moderate differential expression is observed, with significant genes showing log2 fold change > ±1 and -log10 adjusted p-value > 2.

- **B. KMS-34 Resistant vs Control:**  
Stronger expression changes are seen, with many genes exceeding log2 fold change ±2.5 and -log10 p-value > 6.

- **C. KMS-11 vs KMS-34 Resistant:**  
Extensive expression differences, with numerous genes showing log2 fold change > ±6 and high significance (-log10 p-value > 10).

---

### Figure 5 : PCA & UMAP

- The UMAP shows clear separation between KMS-11 and KMS-34 samples, and between control and resistant groups, highlighting distinct expression profiles across both cell lines and treatment conditions.  
- PCA reveals strong clustering by both cell line and resistance status.  
- This distinct separation supports robust differences in gene expression driven by both cell type and resistance development.

## 🔝 Top 5 Differentially Expressed Genes (DEGs)

| **Gene**   | **Biological Role**                                                                                              | **Associated Pathways**                    | **Supporting Reference**   |
|------------|-----------------------------------------------------------------------------------------------------------------|--------------------------------------------|----------------------------|
| **GJA1 (Connexin 43)** | Connexin 43 is a major component of gap junctions enabling intercellular communication and ion exchange. It regulates cell cycle arrest, apoptosis, and differentiation, essential for tissue homeostasis. | Cell cycle regulation, apoptosis, tissue development | Zhang et al., 2014         |
| **S100A4**  | Calcium-binding protein involved in cell migration, invasion, and metastasis. Plays a critical role in cytoskeletal remodeling and immune regulation in cancer progression. | Cancer metastasis, immune response modulation          | Lim et al., 2021           |
| **CD74**   | Acts as a chaperone for MHC class II molecules and receptor for macrophage migration inhibitory factor (MIF). Regulates antigen presentation and immune cell signaling, overexpressed in hematologic malignancies. | Antigen processing and presentation, immune system regulation | Stein et al., 2004         |
| **EPCAM**  | Transmembrane glycoprotein involved in epithelial cell adhesion, Wnt signaling, and tumorigenesis. Overexpressed in various cancers contributing to proliferation and stemness. | Epithelial cell signaling, cancer progression             | Trzpis et al., 2021        |

