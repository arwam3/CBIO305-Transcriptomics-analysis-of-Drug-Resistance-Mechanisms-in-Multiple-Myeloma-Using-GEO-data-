# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
#############################################################################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)
library(openxlsx) 
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(pheatmap)

# Load series and platform data from GEO
gset <- getGEO(filename = "D:/UNIVERSITY/Junior/semester 2/OMICS/project/GSE69078_series_matrix.txt.gz", AnnotGPL = TRUE)

# Fix feature variable names
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Group membership
gsms <- "222333000111"
sml <- strsplit(gsms, split="")[[1]]

# Log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)
}

# Assign groups and design matrix
gs <- factor(sml)
groups <- make.names(c("KMS-11 ctrl","KMS-34 ctrl","KMS-11 resistant","KMS-34 resistant"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

# Remove rows with missing values
gset <- gset[complete.cases(exprs(gset)), ]

# Fit linear model
fit <- lmFit(gset, design)

# CONTRASTS: resistant - control (to get upregulated in carfilzomib)
cont.matrix <- makeContrasts(
  KMS11res_vs_KMS11ctrl = `KMS.11.resistant` - `KMS.11.ctrl`,
  KMS34res_vs_KMS34ctrl = `KMS.34.resistant` - `KMS.34.ctrl`,
  KMS11res_vs_KMS34res = `KMS.11.resistant` - `KMS.34.resistant`,
  KMS11ctrl_vs_KMS34ctrl = `KMS.11.ctrl` - `KMS.34.ctrl`,
  KMS11ctrl_vs_KMS34res = `KMS.11.ctrl` - `KMS.34.resistant`,
  KMS34ctrl_vs_KMS11res = `KMS.34.ctrl` - `KMS.11.resistant`,
  levels = design
)
fit2 <- contrasts.fit(fit, cont.matrix)

# Compute statistics
fit2 <- eBayes(fit2, 0.01)

# Extract results per contrast
base_cols <- c("P.Value", "adj.P.Val", "F", "logFC")
contrast_tables <- list()

for (i in 1:ncol(fit2)) {
  contrast_name <- colnames(fit2)[i]
  tT_contrast <- topTable(fit2, coef=i, adjust="fdr", sort.by="B", number=Inf)
  existing_cols <- intersect(base_cols, colnames(tT_contrast))
  renamed_cols <- paste0(existing_cols, "_", contrast_name)
  colnames(tT_contrast)[match(existing_cols, colnames(tT_contrast))] <- renamed_cols
  cols_to_keep <- c("ID", "Gene.symbol", "Gene.title", renamed_cols)
  cols_to_keep <- intersect(cols_to_keep, colnames(tT_contrast))
  contrast_tables[[contrast_name]] <- tT_contrast[, cols_to_keep]
}

# Merge all contrasts
final_tT <- Reduce(function(x, y) merge(x, y, by=c("ID", "Gene.symbol", "Gene.title"), all=TRUE), contrast_tables)

# Output the full table
write.table(final_tT, file=stdout(), row.names=FALSE, sep="\t")

# Simplified table
cols_of_interest <- grep("logFC_|adj.P.Val_|^ID$|^Gene.symbol$|^Gene.title$", colnames(final_tT), value=TRUE)
simplified_tT <- final_tT[, cols_of_interest]
write.table(simplified_tT, file="simplified_DEG_table.tsv", row.names=FALSE, sep="\t")

# Filter for significant DEGs: abs(logFC) > 1 and adj.P.Val < 0.05
deg_list <- list()
contrast_names <- gsub("logFC_", "", grep("^logFC_", colnames(simplified_tT), value = TRUE))

for (contrast in contrast_names) {
  logFC_col <- paste0("logFC_", contrast)
  adjP_col <- paste0("adj.P.Val_", contrast)
  if (all(c(logFC_col, adjP_col) %in% colnames(simplified_tT))) {
    filtered <- subset(simplified_tT,
                       abs(simplified_tT[[logFC_col]]) > 1 &
                         simplified_tT[[adjP_col]] < 0.05)
    if (nrow(filtered) > 0) {
      filtered$Contrast <- contrast
      deg_list[[contrast]] <- filtered
    }
  }
}

significant_DEGs <- do.call(rbind, deg_list)
significant_DEGs <- significant_DEGs[order(significant_DEGs$Contrast,
                                           -abs(significant_DEGs[[grep("^logFC_", colnames(significant_DEGs))[1]]])), ]

# Export to Excel
output_path <- "D:/UNIVERSITY/Junior/semester 2/OMICS/project/"
write.xlsx(simplified_tT, file = paste0(output_path, "all_DEGs_table.xlsx"), rowNames = FALSE)
write.xlsx(significant_DEGs, file = paste0(output_path, "significant_DEGs_filtered.xlsx"), rowNames = FALSE)

# Keep only key contrasts (carfilzomib resistant vs control)
significant_genes <- subset(final_tT,
                            (abs(logFC_KMS11res_vs_KMS11ctrl) > 1 & adj.P.Val_KMS11res_vs_KMS11ctrl < 0.05) |
                              (abs(logFC_KMS34res_vs_KMS34ctrl) > 1 & adj.P.Val_KMS34res_vs_KMS34ctrl < 0.05) |
                              (abs(logFC_KMS11res_vs_KMS34res) > 1 & adj.P.Val_KMS11res_vs_KMS34res < 0.05)
)

significant_genes_subset <- significant_genes[, c(
  "ID", "Gene.symbol", "Gene.title",
  "logFC_KMS11res_vs_KMS11ctrl", "adj.P.Val_KMS11res_vs_KMS11ctrl", "P.Value_KMS11res_vs_KMS11ctrl",
  "logFC_KMS34res_vs_KMS34ctrl", "adj.P.Val_KMS34res_vs_KMS34ctrl", "P.Value_KMS34res_vs_KMS34ctrl",
  "logFC_KMS11res_vs_KMS34res", "adj.P.Val_KMS11res_vs_KMS34res", "P.Value_KMS11res_vs_KMS34res"
)]

significant_genes_subset$max_abs_logFC <- apply(significant_genes_subset[, grep("^logFC_", names(significant_genes_subset))], 1, function(x) max(abs(x)))
significant_genes_subset <- significant_genes_subset[order(-significant_genes_subset$max_abs_logFC), ]
significant_genes_subset$max_abs_logFC <- NULL

write.xlsx(
  significant_genes_subset,
  file = "D:/UNIVERSITY/Junior/semester 2/OMICS/project/significant_DEGs_filtered.xlsx",
  rowNames = FALSE
)
######################################################################################################################################################################
#Visualization#
# Histogram of adjusted p-values for your 3 main contrasts
adj_pvals <- unlist(final_tT[, c(
  "adj.P.Val_KMS11res_vs_KMS11ctrl",
  "adj.P.Val_KMS34res_vs_KMS34ctrl",
  "adj.P.Val_KMS11res_vs_KMS34res"
)], use.names = FALSE)

# Remove NAs
adj_pvals <- adj_pvals[!is.na(adj_pvals)]

# Plot histogram
hist(adj_pvals, col = "grey", border = "white",
     xlab = "Adjusted P-Value", ylab = "Number of Genes",
     main = "Adjusted P-Value Distribution (Selected Contrasts)",
     breaks = 30)
###################################################################################################################
#Static volcano
# Volcano plot function
plot_volcano <- function(logFC_col, adjP_col, contrast_name) {
  df <- data.frame(
    logFC = final_tT[[logFC_col]],
    adjP = final_tT[[adjP_col]]
  )
  df <- df[complete.cases(df), ]
  
  df$Regulation <- "No Diff"
  df$Regulation[df$adjP < 0.05 & df$logFC > 1] <- "Up"
  df$Regulation[df$adjP < 0.05 & df$logFC < -1] <- "Down"
  
  df$Regulation <- factor(df$Regulation, levels = c("Up", "Down", "No Diff"))
  
  ggplot(df, aes(x = logFC, y = -log10(adjP), color = Regulation)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "No Diff" = "black")) +
    labs(title = paste("Volcano Plot:", contrast_name),
         x = "Log2 Fold Change", y = "-log10(Adjusted P-value)") +
    theme_minimal() +
    theme(legend.title = element_blank())
}

# Example usage:
plot_volcano("logFC_KMS11res_vs_KMS11ctrl", "adj.P.Val_KMS11res_vs_KMS11ctrl", "KMS-11 Resistant vs Control")
plot_volcano("logFC_KMS34res_vs_KMS34ctrl", "adj.P.Val_KMS34res_vs_KMS34ctrl", "KMS-34 Resistant vs Control")
plot_volcano("logFC_KMS11res_vs_KMS34res", "adj.P.Val_KMS11res_vs_KMS34res", "KMS-11 Resistant vs KMS-34 Resistant")
########################################################################################################################################################
# Interactive volcano plot function
create_volcano_plot <- function(data, lfc_col, pval_col, contrast_name) {
  volcano_data <- data[, c("ID", "Gene.symbol", lfc_col, pval_col)]
  colnames(volcano_data) <- c("ID", "Gene", "logFC", "adj.P.Val")
  
  volcano_data <- volcano_data[complete.cases(volcano_data), ]
  volcano_data$negLog10Pval <- -log10(volcano_data$adj.P.Val)
  
  volcano_data$Color <- "Not Significant"
  volcano_data$Color[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.05] <- "Upregulated"
  volcano_data$Color[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.05] <- "Downregulated"
  volcano_data$Color <- factor(volcano_data$Color, levels = c("Upregulated", "Downregulated", "Not Significant"))
  
  p <- ggplot(volcano_data, aes(x = logFC, y = negLog10Pval, text = Gene)) +
    geom_point(aes(color = Color), alpha = 0.7) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black")) +
    labs(
      title = paste("Interactive Volcano Plot:", contrast_name),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-Value",
      color = "Expression"
    ) +
    theme_minimal()
  
  interactive_plot <- ggplotly(p, tooltip = "text")
  return(interactive_plot)
}

# Generate and display plots
plot1 <- create_volcano_plot(final_tT, "logFC_KMS11res_vs_KMS11ctrl", "adj.P.Val_KMS11res_vs_KMS11ctrl", "KMS-11 Resistant vs Control")
plot2 <- create_volcano_plot(final_tT, "logFC_KMS34res_vs_KMS34ctrl", "adj.P.Val_KMS34res_vs_KMS34ctrl", "KMS-34 Resistant vs Control")
plot3 <- create_volcano_plot(final_tT, "logFC_KMS11res_vs_KMS34res", "adj.P.Val_KMS11res_vs_KMS34res", "KMS-11 Resistant vs KMS-34 Resistant")

# Display in RStudio viewer or save
plot1
plot2
plot3
########################################################################################################################################################
# Interactive MD
#Interactive MD (corrected)
# Function to create interactive MD plot using Gene.symbol from final_tT
plot_md_interactive <- function(data, lfc_col, adjp_col, contrast_name) {
  # Extract needed columns
  md_data <- data[, c("ID", "Gene.symbol", lfc_col, adjp_col)]
  colnames(md_data) <- c("ID", "Gene", "logFC", "adj.P.Val")
  
  # Add mean expression from limma fit (Amean is in same row order as final_tT)
  md_data$meanExpr <- fit2$Amean
  
  # Remove NAs
  md_data <- md_data[complete.cases(md_data), ]
  
  # Assign colors based on thresholds
  md_data$Color <- "No Diff"
  md_data$Color[md_data$adj.P.Val < 0.05 & md_data$logFC > 1] <- "Up"
  md_data$Color[md_data$adj.P.Val < 0.05 & md_data$logFC < -1] <- "Down"
  
  color_map <- c("Up" = "red", "Down" = "blue", "No Diff" = "black")
  
  # Create interactive plot
  plot <- plot_ly(
    data = md_data,
    x = ~meanExpr,
    y = ~logFC,
    text = ~paste("Gene:", Gene,
                  "<br>Avg Log Expr:", round(meanExpr, 2),
                  "<br>LogFC:", round(logFC, 2),
                  "<br>Adj.P.Val:", signif(adj.P.Val, 3)),
    type = "scatter",
    mode = "markers",
    color = ~Color,
    colors = color_map,
    marker = list(size = 5),
    hoverinfo = "text"
  ) %>%
    layout(
      title = paste("Interactive MD Plot:", contrast_name),
      xaxis = list(title = "Average Log Expression"),
      yaxis = list(title = "Log Fold Change"),
      shapes = list(list(
        type = "line",
        x0 = min(md_data$meanExpr),
        x1 = max(md_data$meanExpr),
        y0 = 0, y1 = 0,
        line = list(color = "gray", dash = "dash")
      ))
    )
  
  return(plot)
}

# ✅ Generate and display plots
md1 <- plot_md_interactive(final_tT, "logFC_KMS11res_vs_KMS11ctrl", "adj.P.Val_KMS11res_vs_KMS11ctrl", "KMS11res vs KMS11ctrl")
md2 <- plot_md_interactive(final_tT, "logFC_KMS34res_vs_KMS34ctrl", "adj.P.Val_KMS34res_vs_KMS34ctrl", "KMS34res vs KMS34ctrl")
md3 <- plot_md_interactive(final_tT, "logFC_KMS11res_vs_KMS34res", "adj.P.Val_KMS11res_vs_KMS34res", "KMS11res vs KMS34res")


# Display one plot in RStudio
md1
md2
md3

# ✅ Save to HTML
saveWidget(md1, "D:/UNIVERSITY/Junior/semester 2/OMICS/project/KMS11res_vs_KMS11ctrl_MD.html")
saveWidget(md2, "D:/UNIVERSITY/Junior/semester 2/OMICS/project/KMS34res_vs_KMS34ctrl_MD.html")
saveWidget(md3, "D:/UNIVERSITY/Junior/semester 2/OMICS/project/KMS11res_vs_KMS34res_MD.html")
#################################################################################################################################################3
# Venn Diagram:Showing overlap of DEGs (significant genes) for your 3 contrasts. 
library(limma)

# Decide tests on selected contrasts
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=1)

# Subset only your contrasts of interest:
dT_sub <- dT[, c("KMS11res_vs_KMS11ctrl", "KMS34res_vs_KMS34ctrl", "KMS11res_vs_KMS34res")]

# Plot Venn diagram
vennDiagram(dT_sub, circle.col=c("red","blue","green"), main="Venn Diagram of DEGs")
######################################################################################################################
#Boxplot: Showing expression distribution for all samples grouped by your 4 groups
# Order samples by group for plotting
ord <- order(gset$group)

palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02"))  # same colors as UMAP

boxplot(exprs(gset)[, ord], 
        boxwex=0.6, notch=TRUE, outline=FALSE, las=2, col=gset$group[ord],
        main=paste("Expression Boxplot -", annotation(gset)),
        ylab="Expression (log2)",
        xlab="Samples ordered by group")
legend("topleft", legend=levels(gset$group), fill=palette(), bty="n")
########################################################################################################################
# UMAP plot (dimensionality reduction) (exactly like geo2r script)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 5, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=5", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE69078")

#############################################################################################################################
#Heatmap for Significant DEGs Across Samples

# 1. Extract expression matrix from GEO object
exprs_data <- exprs(gset)

# 2. Get probe IDs from your filtered significant genes table
sig_gene_ids <- significant_genes_subset$ID  # OR rownames(significant_genes_subset) if IDs are row names

# 3. Subset expression data
heatmap_data <- exprs_data[sig_gene_ids, ]

# 4. Scale expression data per gene (row)
scaled_data <- t(scale(t(heatmap_data)))

# 5. Prepare sample annotation (make sure `gs` is your group info vector)
annotation_col <- data.frame(Group = gs)
rownames(annotation_col) <- colnames(exprs_data)

# 6. Plot the heatmap
pheatmap(
  scaled_data,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_col = 8,
  main = "Heatmap of Significant DEGs",
  color = colorRampPalette(c("blue", "white", "red"))(100)
)
##############################################################################################################################
# PCA plot
# For PCA, you need actual expression values per sample — not statistics like logFC.
# So, the correct approach is:
# Use exprs(gset) and filter its rows (genes) using the gene IDs or symbols from significant_genes_subset.

# Step 1: Get gene IDs or symbols
deg_ids <- significant_genes_subset$ID  # Or Gene.symbol if applicable

# Step 2: Match them to rownames of expression matrix
expr_data <- exprs(gset)
deg_expr <- expr_data[rownames(expr_data) %in% deg_ids, ]

# Step 3: Do PCA on the transposed expression matrix
pca <- prcomp(t(deg_expr), scale. = TRUE)

# Step 4: Plot
pca_df <- data.frame(pca$x[, 1:2])
pca_df$Group <- gset$group
pca_df$Sample <- colnames(deg_expr)

library(ggplot2)
library(plotly)

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group,
                        text = paste("Sample:", Sample, "<br>Group:", Group))) +
  geom_point(size = 3) +
  labs(title = "PCA of Significant DEGs", x = "PC1", y = "PC2") +
  theme_minimal()

interactive_pca <- ggplotly(p, tooltip = "text")
interactive_pca

#Save as HTML
saveWidget(interactive_pca, 
           "D:/UNIVERSITY/Junior/semester 2/OMICS/project/PCA_Significant_DEGs.html")






