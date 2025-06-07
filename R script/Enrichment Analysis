# =========================================================================================
# STEP 1: Load Required Libraries
# =========================================================================================
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(c("GEOquery", "limma", "openxlsx", "clusterProfiler", "org.Hs.eg.db", "enrichplot"), ask = FALSE)
install.packages(c(
"openxlsx",   # for Excel export
"ggplot2",    # base for plotting
"plotly",     # makes plots interactive
"htmlwidgets",
"STRINGdb", 
"pathview"# to save plotly as HTML
))
# Load libraries
library(GEOquery)         # For downloading GEO data
library(limma)            # For differential expression analysis
library(openxlsx)         # For exporting Excel files
library(clusterProfiler)  # For enrichment analysis
library(org.Hs.eg.db)     # Gene annotation for humans
library(enrichplot)       # For visualizing enrichment results
library(umap)             # Uniform Manifold Approximation and Projection (UMAP) for dimensionality reduction
library(ggplot2)          # Powerful and flexible plotting package based on the grammar of graphics
library(plotly)           # For creating interactive plots (can convert ggplot2 plots to interactive)
library(htmlwidgets)      # Needed to save interactive widgets like plotly plots as HTML
library(pheatmap)         # For creating enhanced heatmaps
# =========================================================================================
# STEP 2: Download GEO Dataset
# =========================================================================================
gset_list <- getGEO("GSE69078", GSEMatrix = TRUE, AnnotGPL = TRUE)
gset <- gset_list[[1]]

# Make feature names syntactically valid for R
fvarLabels(gset) <- make.names(fvarLabels(gset))
# Group membership string for each sample (each character is a group code)
gsms <- "222333000111"

# Split the string into a vector of single characters representing group IDs
sml <- strsplit(gsms, split="")[[1]]

# Extract expression matrix from the ExpressionSet (genes x samples)
ex <- exprs(gset)

# Calculate quantiles of expression values to assess data distribution
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=TRUE))

# Check if log2 transformation is needed:
# - Condition 1: 99th percentile > 100 (very high values likely not log-transformed)
# - OR Condition 2: range > 50 and 25th percentile > 0 (large spread, likely raw counts)
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

if (LogC) {
  # Replace any zero or negative values with NaN to avoid log2 issues
  ex[which(ex <= 0)] <- NaN
  
  # Apply log2 transformation to the expression matrix
  exprs(gset) <- log2(ex)
}
# --- Assign samples to groups ---
gs <- factor(sml)  # Convert group membership vector to a factor

# Define meaningful group names and make them syntactically valid variable names
groups <- make.names(c("KMS-11 ctrl","KMS-34 ctrl","KMS-11 resistant","KMS-34 resistant"))

# Assign these group names as factor levels
levels(gs) <- groups

# Store the group factor in the ExpressionSet phenotype data
gset$group <- gs

# --- Create design matrix for linear modeling ---
# Model without intercept (~group + 0) so each group has its own coefficient
design <- model.matrix(~group + 0, gset)

# Rename columns of design matrix to match group names exactly
colnames(design) <- levels(gs)

# --- Clean the dataset ---
# Remove genes with missing values to avoid errors in linear modeling
gset <- gset[complete.cases(exprs(gset)), ]

# --- Fit linear model to expression data ---
fit <- lmFit(gset, design)  # Fit gene-wise linear models

# --- Set up contrasts to compare groups ---
# Define pairwise comparisons of interest using makeContrasts
cont.matrix <- makeContrasts(
  KMS11ctrl_vs_KMS34ctrl = `KMS.11.ctrl` - `KMS.34.ctrl`,
  KMS11ctrl_vs_KMS11res = `KMS.11.ctrl` - `KMS.11.resistant`,
  KMS11ctrl_vs_KMS34res = `KMS.11.ctrl` - `KMS.34.resistant`,
  KMS34ctrl_vs_KMS11res = `KMS.34.ctrl` - `KMS.11.resistant`,
  KMS34ctrl_vs_KMS34res = `KMS.34.ctrl` - `KMS.34.resistant`,
  KMS11res_vs_KMS34res = `KMS.11.resistant` - `KMS.34.resistant`,
  levels = design
)

# Apply contrast matrix to the fitted model
fit2 <- contrasts.fit(fit, cont.matrix)

# --- Empirical Bayes moderation to improve statistics ---
fit2 <- eBayes(fit2, 0.01)  # Shrinkage of variance estimates

# --- Prepare to extract top genes from each contrast ---
base_cols <- c("P.Value", "adj.P.Val", "F", "logFC")  # Stats of interest

contrast_tables <- list()  # To store results for each contrast

# Loop through each contrast column in the model fit
for (i in 1:ncol(fit2)) {
  contrast_name <- colnames(fit2)[i]
  
  # Extract top table for this contrast, sorted by B-statistic
  tT_contrast <- topTable(fit2, coef=i, adjust="fdr", sort.by="B", number=Inf)
  
  # Identify columns present in this table that we want to rename
  existing_cols <- intersect(base_cols, colnames(tT_contrast))
  
  # Rename these columns to append contrast name suffix (to avoid collisions)
  renamed_cols <- paste0(existing_cols, "_", contrast_name)
  colnames(tT_contrast)[match(existing_cols, colnames(tT_contrast))] <- renamed_cols
  
  # Select only relevant columns (gene IDs and stats)
  cols_to_keep <- c("ID", "Gene.symbol", "Gene.title", renamed_cols)
  cols_to_keep <- intersect(cols_to_keep, colnames(tT_contrast))  # Safety check
  
  # Store in list
  contrast_tables[[contrast_name]] <- tT_contrast[, cols_to_keep]
}

# --- Merge all contrasts by common gene identifiers ---
final_tT <- Reduce(function(x, y) merge(x, y, by=c("ID", "Gene.symbol", "Gene.title"), all=TRUE), contrast_tables)

# --- Output the complete table of DE results ---
write.table(final_tT, file=stdout(), row.names=FALSE, sep="\t")

# --- Create simplified table with only IDs and logFC + adj.P.Val columns ---
cols_of_interest <- grep("logFC_|adj.P.Val_|^ID$|^Gene.symbol$|^Gene.title$", colnames(final_tT), value=TRUE)
simplified_tT <- final_tT[, cols_of_interest]

# Save simplified differential expression table to disk
write.table(simplified_tT, file="simplified_DEG_table.tsv", row.names=FALSE, sep="\t")
###################################################################################################
# Filter for significant DEGs: genes with |log2 fold change| > 1 and adjusted p-value < 0.05 in any contrast

deg_list <- list()  # Initialize empty list to store filtered DEGs from each contrast

# Extract contrast names from column names starting with "logFC_"
contrast_names <- gsub("logFC_", "", grep("^logFC_", colnames(simplified_tT), value = TRUE))

# Loop through each contrast to filter significant DEGs
for (contrast in contrast_names) {
  logFC_col <- paste0("logFC_", contrast)       # Construct column name for log fold change
  adjP_col <- paste0("adj.P.Val_", contrast)   # Construct column name for adjusted p-value
  
  # Check if both columns exist (safety check)
  if (all(c(logFC_col, adjP_col) %in% colnames(simplified_tT))) {
    
    # Subset genes meeting criteria: abs(logFC) > 1 and adj.P.Val < 0.05
    filtered <- subset(simplified_tT,
                       abs(simplified_tT[[logFC_col]]) > 1 &
                         simplified_tT[[adjP_col]] < 0.05)
    
    # If any genes pass filtering, add a column to identify contrast and save to list
    if (nrow(filtered) > 0) {
      filtered$Contrast <- contrast
      deg_list[[contrast]] <- filtered
    }
  }
}

# Combine all filtered DEGs across contrasts into one data.frame
significant_DEGs <- do.call(rbind, deg_list)

# Order significant DEGs by Contrast and descending absolute logFC for clarity
significant_DEGs <- significant_DEGs[order(significant_DEGs$Contrast,
                                           -abs(significant_DEGs[[grep("^logFC_", colnames(significant_DEGs))[1]]])), ]


###################################################################################################
# Exporting simplified_tT and significant_DEGs tables to Excel files

# Define output folder path for saving files
output_path <- "E:/NU/5th semester/Omics I/Project/"

# Save the full simplified differential expression table (all genes, all contrasts)
write.xlsx(simplified_tT, file = paste0(output_path, "all_DEGs_table.xlsx"), rowNames = FALSE)

# Save the filtered significant DEGs table (only genes passing criteria)
write.xlsx(significant_DEGs, file = paste0(output_path, "significant_DEGs_filtered.xlsx"), rowNames = FALSE)


########################################################################################################
# Filter for significant genes using only specific contrasts and export final subset

significant_genes <- subset(final_tT,
                            (abs(logFC_KMS11ctrl_vs_KMS11res) > 1 & adj.P.Val_KMS11ctrl_vs_KMS11res < 0.05) |
                              (abs(logFC_KMS34ctrl_vs_KMS34res) > 1 & adj.P.Val_KMS34ctrl_vs_KMS34res < 0.05) |
                              (abs(logFC_KMS11res_vs_KMS34res) > 1 & adj.P.Val_KMS11res_vs_KMS34res < 0.05)
)

# Select columns for gene IDs, symbols, titles, and stats of the selected contrasts
significant_genes_subset <- significant_genes[, c(
  "ID", "Gene.symbol", "Gene.title",
  "logFC_KMS11ctrl_vs_KMS11res", "adj.P.Val_KMS11ctrl_vs_KMS11res", "P.Value_KMS11ctrl_vs_KMS11res",
  "logFC_KMS34ctrl_vs_KMS34res", "adj.P.Val_KMS34ctrl_vs_KMS34res", "P.Value_KMS34ctrl_vs_KMS34res",
  "logFC_KMS11res_vs_KMS34res", "adj.P.Val_KMS11res_vs_KMS34res", "P.Value_KMS11res_vs_KMS34res"
)]

# Add helper column for max absolute logFC across these contrasts to sort by effect size
significant_genes_subset$max_abs_logFC <- apply(significant_genes_subset[, grep("^logFC_", names(significant_genes_subset))], 1, max.abs <- function(x) max(abs(x)))

# Order descending by maximum absolute logFC (strongest changes at top)
significant_genes_subset <- significant_genes_subset[order(-significant_genes_subset$max_abs_logFC), ]

# Remove helper column before exporting
significant_genes_subset$max_abs_logFC <- NULL

# Export final filtered and ordered DEG table to Excel
library(openxlsx)
write.xlsx(
  significant_genes_subset,
  file = "E:/NU/5th semester/Omics I/Project/significant_DEGs_filtered.xlsx",
  rowNames = FALSE
)

# ======================================================================================
# Extract gene symbols from your significant DEGs table

gene_symbols <- significant_genes_subset$Gene.symbol

# ======================================================================================
# Convert gene symbols to Entrez IDs
# ======================================================================================
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", 
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# ======================================================================================
# Define GO enrichment function
# ======================================================================================
run_go_enrichment <- function(entrez_gene_list, ontology) {
  enrichGO(
    gene          = entrez_gene_list,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
}

# Run enrichment for BP, MF, CC
go_BP <- run_go_enrichment(entrez_ids$ENTREZID, "BP")
go_MF <- run_go_enrichment(entrez_ids$ENTREZID, "MF")
go_CC <- run_go_enrichment(entrez_ids$ENTREZID, "CC")
# ======================================================================================
# Visualizations: Barplot and Dotplot
# ======================================================================================
bp_plot <- dotplot(go_BP, showCategory = 10, title = "GO BP Enrichment")
mf_plot <- dotplot(go_MF, showCategory = 10, title = "GO MF Enrichment")
cc_plot <- dotplot(go_CC, showCategory = 10, title = "GO CC Enrichment")

# Make interactive
bp_plotly <- ggplotly(bp_plot)
mf_plotly <- ggplotly(mf_plot)
cc_plotly <- ggplotly(cc_plot)

# Save HTML widgets
saveWidget(bp_plotly, "GO_BP_Enrichment.html")
saveWidget(mf_plotly, "GO_MF_Enrichment.html")
saveWidget(cc_plotly, "GO_CC_Enrichment.html")

# ======================================================================================
# Tree Plots (emapplot)
# ======================================================================================
go_BP <- pairwise_termsim(go_BP)
go_MF <- pairwise_termsim(go_MF)
go_CC <- pairwise_termsim(go_CC)

emapplot(go_BP, showCategory = 10, title = "GO BP Tree")
emapplot(go_MF, showCategory = 10, title = "GO MF Tree")
emapplot(go_CC, showCategory = 10, title = "GO CC Tree")

# ======================================================================================
# Gene-Concept Network (cnetplot)
# ======================================================================================
cnetplot(go_BP, categorySize = "pvalue", foldChange = NULL, showCategory = 5)
cnetplot(go_MF, categorySize = "pvalue", foldChange = NULL, showCategory = 5)
cnetplot(go_CC, categorySize = "pvalue", foldChange = NULL, showCategory = 5)

# ======================================================================================
# STRING PPI Network
# ======================================================================================
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)
mapped <- string_db$map(significant_genes_subset, "Gene.symbol", removeUnmappedRows = TRUE)
string_db$plot_network(mapped$STRING_id)

# ======================================================================================
# KEGG Pathway Enrichment
# ======================================================================================
kegg <- enrichKEGG(gene = entrez_ids$ENTREZID,
                   organism = "hsa",
                   pvalueCutoff = 0.05)

dotplot(kegg, showCategory = 10, title = "KEGG Pathway Enrichment")
saveWidget(ggplotly(dotplot(kegg)), "KEGG_Enrichment.html")

# =============================
# KEGG Pathway Diagram (Pathview)
# =============================
# Provide dummy FC values (optional: use logFC from DEGs)
fc <- setNames(rep(1, length(entrez_ids$ENTREZID)), entrez_ids$ENTREZID)
if (length(kegg@result$ID) > 0) {
  pathview(gene.data = fc, pathway.id = kegg@result$ID[1], species = "hsa")
}
