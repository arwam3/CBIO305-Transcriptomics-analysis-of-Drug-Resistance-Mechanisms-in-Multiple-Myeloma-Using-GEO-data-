# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma

library(GEOquery)
library(limma)
library(umap)
library(openxlsx) 
library(ggplot2)
library(umap)
library(plotly)
library(htmlwidgets)
library(pheatmap)

# load series and platform data from GEO

#gset <- getGEO("GSE69078", GSEMatrix =TRUE, AnnotGPL=TRUE)
#if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
#gset <- gset[[idx]]

# Adjust the path based on where you saved the file
gset <- getGEO(filename = "D:/UNIVERSITY/Junior/semester 2/OMICS/project/GSE69078_series_matrix.txt.gz", AnnotGPL = TRUE)



# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "222333000111"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("KMS-11 ctrl","KMS-34 ctrl","KMS-11 resistant","KMS-34 resistant"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
#cts <- paste(groups, c(tail(groups, -1), head(groups, 1)), sep="-")
#cont.matrix <- makeContrasts(contrasts=cts, levels=design)
#fit2 <- contrasts.fit(fit, cont.matrix)

# set up all pairwise contrasts of interest and recalculate model coefficients
cont.matrix <- makeContrasts(
  KMS11ctrl_vs_KMS34ctrl = `KMS.11.ctrl` - `KMS.34.ctrl`,
  KMS11ctrl_vs_KMS11res = `KMS.11.ctrl` - `KMS.11.resistant`,
  KMS11ctrl_vs_KMS34res = `KMS.11.ctrl` - `KMS.34.resistant`,
  KMS34ctrl_vs_KMS11res = `KMS.34.ctrl` - `KMS.11.resistant`,
  KMS34ctrl_vs_KMS34res = `KMS.34.ctrl` - `KMS.34.resistant`,
  KMS11res_vs_KMS34res = `KMS.11.resistant` - `KMS.34.resistant`,
  levels = design
)
fit2 <- contrasts.fit(fit, cont.matrix)

####################################################################################################

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)

# Define the columns we want to extract and rename from each topTable
base_cols <- c("P.Value", "adj.P.Val", "F", "logFC")

# Create a list to store each contrast's table
contrast_tables <- list()

# Loop over each contrast
for (i in 1:ncol(fit2)) {
  contrast_name <- colnames(fit2)[i]
  
  # Get topTable for the i-th contrast
  tT_contrast <- topTable(fit2, coef=i, adjust="fdr", sort.by="B", number=Inf)
  
  # Ensure required columns are present before renaming
  existing_cols <- intersect(base_cols, colnames(tT_contrast))
  renamed_cols <- paste0(existing_cols, "_", contrast_name)
  colnames(tT_contrast)[match(existing_cols, colnames(tT_contrast))] <- renamed_cols
  
  # Keep only needed columns
  cols_to_keep <- c("ID", "Gene.symbol", "Gene.title", renamed_cols)
  cols_to_keep <- intersect(cols_to_keep, colnames(tT_contrast))  # safeguard
  contrast_tables[[contrast_name]] <- tT_contrast[, cols_to_keep]
}

# Merge all contrast tables by ID
final_tT <- Reduce(function(x, y) merge(x, y, by=c("ID", "Gene.symbol", "Gene.title"), all=TRUE), contrast_tables)

# Output the full table
write.table(final_tT, file=stdout(), row.names=FALSE, sep="\t")

#It only renames columns if they exist.
#It merges all contrasts correctly by shared gene IDs.
#It gives you separate logFC, P.Value, adj.P.Val, and F values per contrast.
###################################################################################################
#Simplified table
# Extract only gene identifiers + logFC and adj.P.Val columns
cols_of_interest <- grep("logFC_|adj.P.Val_|^ID$|^Gene.symbol$|^Gene.title$", colnames(final_tT), value=TRUE)

# Create simplified table
simplified_tT <- final_tT[, cols_of_interest]

# Output the simplified version
write.table(simplified_tT, file="simplified_DEG_table.tsv", row.names=FALSE, sep="\t")

###################################################################################################
#Filter for significant DEGs: (Where absolute log2 fold change > 1 and adjusted p-value < 0.05 in any contrast)

# Create an empty list to collect filtered DEGs from each contrast
deg_list <- list()

# Get all contrast names from the logFC columns
contrast_names <- gsub("logFC_", "", grep("^logFC_", colnames(simplified_tT), value = TRUE))

# Filter significant DEGs per contrast
for (contrast in contrast_names) {
  logFC_col <- paste0("logFC_", contrast)
  adjP_col <- paste0("adj.P.Val_", contrast)
  
  # Check if both columns exist (safety)
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

# Combine all filtered DEGs into one table
significant_DEGs <- do.call(rbind, deg_list)

# Order by absolute logFC (highest first) within each contrast
significant_DEGs <- significant_DEGs[order(significant_DEGs$Contrast,
                                           -abs(significant_DEGs[[grep("^logFC_", colnames(significant_DEGs))[1]]])), ]

#NOTES
# It keeps all genes that are significant in any contrast.
#Genes not significant in any contrast get excluded.
#For genes that are significant in some contrasts but not others, the entire row is kept, showing all contrasts' values.
###################################################################################################
#Exporting simplified_tT & significant_DEGs to excel
# Define your output path
output_path <- "D:/UNIVERSITY/Junior/semester 2/OMICS/project/"

# Save full simplified table
write.xlsx(simplified_tT, file = paste0(output_path, "all_DEGs_table.xlsx"), rowNames = FALSE)

# Save significant filtered DEGs
write.xlsx(significant_DEGs, file = paste0(output_path, "significant_DEGs_filtered.xlsx"), rowNames = FALSE)

###################################################################################################
# #ONLY keep needed contrasts (NO pvalue)
# # Define which contrasts to keep
# #selected_contrasts <- c(
#   "logFC_KMS11ctrl_vs_KMS11res",
#   "adj.P.Val_KMS11ctrl_vs_KMS11res",
#   "logFC_KMS34ctrl_vs_KMS34res",
#   "adj.P.Val_KMS34ctrl_vs_KMS34res",
#   "logFC_KMS11res_vs_KMS34res",
#   "adj.P.Val_KMS11res_vs_KMS34res"
# )
# 
# # keep gene identification columns
# #id_columns <- c("ID", "Gene.symbol", "Gene.title")
# 
# # Subset the merged table
# #selected_DEGs <- final_tT[, c(id_columns, selected_contrasts)]
# 
# # OPTIONAL: Filter for significance (|logFC| > 1 & adj.P.Val < 0.05) in *any* contrast
# #sig_DEGs <- subset(selected_DEGs,
#                    #(abs(logFC_KMS11ctrl_vs_KMS11res) > 1 & adj.P.Val_KMS11ctrl_vs_KMS11res < 0.05) |
#                      #(abs(logFC_KMS34ctrl_vs_KMS34res) > 1 & adj.P.Val_KMS34ctrl_vs_KMS34res < 0.05) |
#                      #(abs(logFC_KMS11res_vs_KMS34res) > 1 & adj.P.Val_KMS11res_vs_KMS34res < 0.05)
# #)
# 
# # Export to Excel if desired
# #library(openxlsx)
# #write.xlsx(sig_DEGs, file = "D:/UNIVERSITY/Junior/semester 2/OMICS/project/selected_significant_DEGs.xlsx", row.names = FALSE)
########################################################################################################
#ONLY keep needed contrasts (With pvalue)
# Filter for significant genes in any of the selected contrasts
significant_genes <- subset(final_tT,
                            (abs(logFC_KMS11ctrl_vs_KMS11res) > 1 & adj.P.Val_KMS11ctrl_vs_KMS11res < 0.05) |
                              (abs(logFC_KMS34ctrl_vs_KMS34res) > 1 & adj.P.Val_KMS34ctrl_vs_KMS34res < 0.05) |
                              (abs(logFC_KMS11res_vs_KMS34res) > 1 & adj.P.Val_KMS11res_vs_KMS34res < 0.05)
)

# Select only the necessary columns
significant_genes_subset <- significant_genes[, c(
  "ID", "Gene.symbol", "Gene.title",
  "logFC_KMS11ctrl_vs_KMS11res", "adj.P.Val_KMS11ctrl_vs_KMS11res", "P.Value_KMS11ctrl_vs_KMS11res",
  "logFC_KMS34ctrl_vs_KMS34res", "adj.P.Val_KMS34ctrl_vs_KMS34res", "P.Value_KMS34ctrl_vs_KMS34res",
  "logFC_KMS11res_vs_KMS34res", "adj.P.Val_KMS11res_vs_KMS34res", "P.Value_KMS11res_vs_KMS34res"
)]

# Order by highest absolute logFC across the 3 contrasts
significant_genes_subset$max_abs_logFC <- apply(significant_genes_subset[, grep("^logFC_", names(significant_genes_subset))], 1, function(x) max(abs(x)))
significant_genes_subset <- significant_genes_subset[order(-significant_genes_subset$max_abs_logFC), ]
significant_genes_subset$max_abs_logFC <- NULL  # Remove helper column

# Export to Excel
library(openxlsx)
write.xlsx(
  significant_genes_subset,
  file = "D:/UNIVERSITY/Junior/semester 2/OMICS/project/significant_DEGs_filtered.xlsx",
  rowNames = FALSE
)

########################################################################################################
# Visualize and quality control test results.

#Histogram of adjusted p-values for your 3 contrasts
# Check and convert to numeric, removing NAs
# Select adjusted p-value columns for your 3 contrasts
adj_pvals <- unlist(final_tT[, c(
  "adj.P.Val_KMS11ctrl_vs_KMS11res",
  "adj.P.Val_KMS34ctrl_vs_KMS34res",
  "adj.P.Val_KMS11res_vs_KMS34res"
)])

# Remove NAs
adj_pvals <- adj_pvals[!is.na(adj_pvals)]

hist(adj_pvals, col = "grey", border = "white",
     xlab = "Adjusted P-Value", ylab = "Number of genes",
     main = "Adjusted P-Value Distribution (Selected Contrasts)",
     breaks = 30)
#############################################################################################################
# #Volcano plot for each of the 3 contrasts
# library(ggplot2)
# 
# plot_volcano <- function(logFC_col, adjP_col, contrast_name) {
#   df <- data.frame(
#     logFC = final_tT[[logFC_col]],
#     adjP = final_tT[[adjP_col]]
#   )
#   df <- df[complete.cases(df), ]
#   
#   df$Regulation <- "No Diff"
#   df$Regulation[df$adjP < 0.05 & df$logFC > 1] <- "Up"
#   df$Regulation[df$adjP < 0.05 & df$logFC < -1] <- "Down"
#   
#   df$Regulation <- factor(df$Regulation, levels = c("Up", "Down", "No Diff"))
#   
#   ggplot(df, aes(x = logFC, y = -log10(adjP), color = Regulation)) +
#     geom_point(alpha = 0.6, size = 1.5) +
#     scale_color_manual(values = c("Up" = "red", "Down" = "blue", "No Diff" = "black")) +
#     labs(title = paste("Volcano Plot:", contrast_name),
#          x = "Log Fold Change", y = "-log10(Adjusted P-value)") +
#     theme_minimal() +
#     theme(legend.title = element_blank())
# }
# 
# # Example usage:
# plot_volcano("logFC_KMS11ctrl_vs_KMS11res", "adj.P.Val_KMS11ctrl_vs_KMS11res", "KMS11ctrl vs KMS11res")
# plot_volcano("logFC_KMS34ctrl_vs_KMS34res", "adj.P.Val_KMS34ctrl_vs_KMS34res", "KMS34ctrl vs KMS34res")
# plot_volcano("logFC_KMS11res_vs_KMS34res", "adj.P.Val_KMS11res_vs_KMS34res", "KMS11res vs KMS34res")
####################################################################################################################
#Interactive Volcano
#Helper
create_volcano_plot <- function(data, lfc_col, pval_col, contrast_name) {
  volcano_data <- data[, c("ID", "Gene.symbol", lfc_col, pval_col)]
  colnames(volcano_data) <- c("ID", "Gene", "logFC", "adj.P.Val")
  
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
#Generate and display:
# KMS11 ctrl vs res
plot1 <- create_volcano_plot(final_tT, "logFC_KMS11ctrl_vs_KMS11res", "adj.P.Val_KMS11ctrl_vs_KMS11res", "KMS11ctrl vs KMS11res")
plot1

# KMS34 ctrl vs res
plot2 <- create_volcano_plot(final_tT, "logFC_KMS34ctrl_vs_KMS34res", "adj.P.Val_KMS34ctrl_vs_KMS34res", "KMS34ctrl vs KMS34res")
plot2

# KMS11res vs KMS34res
plot3 <- create_volcano_plot(final_tT, "logFC_KMS11res_vs_KMS34res", "adj.P.Val_KMS11res_vs_KMS34res", "KMS11res vs KMS34res")
plot3

#Save as interactive HTML
# Save all volcano plots to your OMICS project folder
saveWidget(plot1, "D:/UNIVERSITY/Junior/semester 2/OMICS/project/KMS11ctrl_vs_KMS11res_volcano.html")
saveWidget(plot2, "D:/UNIVERSITY/Junior/semester 2/OMICS/project/KMS34ctrl_vs_KMS34res_volcano.html")
saveWidget(plot3, "D:/UNIVERSITY/Junior/semester 2/OMICS/project/KMS11res_vs_KMS34res_volcano.html")

#############################################################################################################
# #MD plot
# plot_md <- function(logFC_col, adjP_col, contrast_name) {
#   logFC <- final_tT[[logFC_col]]
#   adjP <- final_tT[[adjP_col]]
#   meanExpr <- fit2$Amean  # mean log-expression from limma fit
#   
#   valid <- complete.cases(logFC, adjP, meanExpr)
#   logFC <- logFC[valid]
#   adjP <- adjP[valid]
#   meanExpr <- meanExpr[valid]
#   
#   # Define colors
#   colors <- rep("black", length(logFC))
#   colors[adjP < 0.05 & logFC > 1] <- "red"
#   colors[adjP < 0.05 & logFC < -1] <- "blue"
#   
#   plot(meanExpr, logFC,
#        pch = 20, col = colors,
#        xlab = "Average Log Expression", ylab = "Log Fold Change",
#        main = paste("MD Plot:", contrast_name))
#   abline(h = 0, col = "gray", lty = 2)
#   
#   legend("topright", legend = c("Up", "Down", "No Diff"),
#          col = c("red", "blue", "black"), pch = 20)
# }
# 
# # Example usage:
# plot_md("logFC_KMS11ctrl_vs_KMS11res", "adj.P.Val_KMS11ctrl_vs_KMS11res", "KMS11ctrl vs KMS11res")
# plot_md("logFC_KMS34ctrl_vs_KMS34res", "adj.P.Val_KMS34ctrl_vs_KMS34res", "KMS34ctrl vs KMS34res")
# plot_md("logFC_KMS11res_vs_KMS34res", "adj.P.Val_KMS11res_vs_KMS34res", "KMS11res vs KMS34res")
################################################################################################################################
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
md1 <- plot_md_interactive(final_tT, "logFC_KMS11ctrl_vs_KMS11res", "adj.P.Val_KMS11ctrl_vs_KMS11res", "KMS11ctrl vs KMS11res")
md2 <- plot_md_interactive(final_tT, "logFC_KMS34ctrl_vs_KMS34res", "adj.P.Val_KMS34ctrl_vs_KMS34res", "KMS34ctrl vs KMS34res")
md3 <- plot_md_interactive(final_tT, "logFC_KMS11res_vs_KMS34res", "adj.P.Val_KMS11res_vs_KMS34res", "KMS11res vs KMS34res")

# Display one plot in RStudio
md1
md2
md3

# ✅ Save to HTML
saveWidget(md1, "D:/UNIVERSITY/Junior/semester 2/OMICS/project/KMS11ctrl_vs_KMS11res_MD.html")
saveWidget(md2, "D:/UNIVERSITY/Junior/semester 2/OMICS/project/KMS34ctrl_vs_KMS34res_MD.html")
saveWidget(md3, "D:/UNIVERSITY/Junior/semester 2/OMICS/project/KMS11res_vs_KMS34res_MD.html")
 
########################################################################################################
# Venn Diagram:Showing overlap of DEGs (significant genes) for your 3 contrasts. 
library(limma)

# Decide tests on selected contrasts
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=1)

# Subset only your contrasts of interest:
dT_sub <- dT[, c("KMS11ctrl_vs_KMS11res", "KMS34ctrl_vs_KMS34res", "KMS11res_vs_KMS34res")]

# Plot Venn diagram
vennDiagram(dT_sub, circle.col=c("red","blue","green"), main="Venn Diagram of DEGs")
#######################################################################################################
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

#######################################################################################################
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

#################################################################################################################
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
