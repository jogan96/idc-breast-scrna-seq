# -----------------------------------------------------------------------------------------
# scRNA-seq analysis: human invasive ductal carcinoma (10x genome)
# Author: Jessica Ogan
# Date: 18/10/2026
# Input: raw_data/filtered.feature_bc_matrix
# Outputs: plots, cluster_markers_all.csv, cluster_markers_top10, DEG.csv


# Install Packages
# BiocManager::install(c("clusterProfiler", "enrichplot", "org.Hs.eg.db"), update = FALSE, ask = FALSE)

# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SingleR)
library(celldex)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)

#Allow for reproducibility
set.seed(1234)


# Extract zipped files
untar("raw_data/Breast_Cancer_3p_filtered_feature_bc_matrix.tar.gz", exdir = "raw_data")
list.files("raw_data", recursive = TRUE)

#Check files are in the raw_data folder
list.files("raw_data")


### Step 1: Load files
raw_data <- Read10X(data.dir = "raw_data/filtered_feature_bc_matrix")

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = raw_data,
  project = "IDC_breast",
  min.cells = 3,       # Keep genes detected in at least 3 cells
  min.features = 200   # Keep cells with at least 200 genes
)

# Add metadata (manually specify tumour type and sample)
seurat_obj$tumour_type <- "IDC"
seurat_obj$sample <- "IDC_sample1"

#Basic summary
seurat_obj

# View metadata
head(seurat_obj@meta.data)


### Step 2: Quality Control ##
# QC performed to remove technical noise and ensure only biologically meaningful cells are used

#2a. Calculate the percentage of mitochondrial genes per cell (human: "^MT-")
# High MT% often indicate apoptotic or low quality cells
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Cross check if MT genes are in the data
grep("^MT-", rownames(seurat_obj), value = TRUE)

#2b. Visualise QC metrics
# nFeature_RNA: no. of genes detected per cell
# nCount_RNA: total RNA counts per cell
# percent.mt: mitochondrial gene percentage per cell
qc_plot <- VlnPlot(seurat_obj,
                   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   pt.size = 0.1,
                   ncol = 3)
qc_plot

ggsave("plots/qc_metrics.png", plot = qc_plot, width = 12, height = 5, dpi = 300)

#2c. Filter out low-quality cells based on:
# - cells expressing fewer than 200 genes (likely dead)
# - cells expressing more than 9000 genes (potential doublets)
# - cells with >15% mitochondrial gene content (stressed/dead)
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > 200 &
                       nFeature_RNA < 9000 &
                       percent.mt < 15)


### Step 3: Normalisation and Feature Selection ##
# Raw UMI counts are biased by sequencing depth 
# Normalise so cells with different sequencing depth can be compared

#3a. Normalise the data using log-normalisation:
# Transforms raw UMI counts into normalised values to account for sequencing depth
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

#3b. Selection of variable genes
# Only highly variable genes are informative for clustering and downstream analysis
# Identify 2000 most variable genes across all cells
# These genes will be used for downstream dimensionality reduction
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

#3c. Visualise the variable features
# The top 20 variable genes are labelled for interpretability
var_plot <- VariableFeaturePlot(seurat_obj)
LabelPoints(plot = var_plot, points = head(VariableFeatures(seurat_obj), 20), repel = TRUE)


### Step 4: Dimensionality Reduction & Clustering ### 
#4a. Scale and centre the expression values of each gene across all cells
# Scaling because PCA assumes normally distributed features and centres gene expression
# Essential step for PCA and downstream linear modelling
seurat_obj <- ScaleData(seurat_obj)

#4b. Run PCA to reduce dimensionality and identify major sources of variation
# PCA reduces noise and retains the most biological signal in fewer dimensions
# The number of components is left at the Seurat default (50)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

#4c. Visualise Elbow Plot
# Visualise the standard deviation of each PC to decide how many PCs to retain
# A sharp drop in the plot helps determine the optimal number of PCs
ElbowPlot(seurat_obj, ndims = 50)

# Based on the elbow plot, we select the first 15 principal components
pcs <- 15 #this can be tuned depending on the dataset

#4d. Identify the k-nearest neighbours for each cell based on the selected PCs
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pcs)

# Perform graph-based clustering of the cells
# Clustering groups cells with similar expression patterns (possible cell types)
# Resolution affects the granularity: higher = more clusters
seurat_obj <- FindClusters(seurat_obj, resolution = 0.6)

#4e. Run UMAP for 2D visualisation of the clusters 
# For easy & interpretable 2D visualisation of complex transcriptomics differences
# UMAP works best for low dimensional PCs
seurat_obj <- RunUMAP(seurat_obj, dims = 1:pcs)

# Plot UMAP clustering result
DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE) + ggtitle("UMAP: Clustered Cells")

# Cross check output
ncol(seurat_obj) # number of cells remaining after QC
summary(seurat_obj$percent.mt) # mitochondrial % should all be below the 15% cutoff
max(seurat_obj$nFeature_RNA) # highest gene count should be below the 9,000 cutoff
table(Idents(seurat_obj)) # number of cells in each cluster (should sum to the total)


### Step 5: Marker Identification (data-driven, before naming any cluster)
# Find genes that distinguish each cluster from all the others
Idents(seurat_obj) <- "seurat_clusters"

markers <- FindAllMarkers(seurat_obj,
                          only.pos = TRUE,      # keep only genes higher in the cluster
                          min.pct = 0.25,       # expressed in at least 25% of cells in the cluster
                          logfc.threshold = 0.25)

# Check that the table is from this analysis
unique(markers$cluster)
head(markers$gene)

# Top 10 specific markers per cluster:
# pct.1 = fraction of cells in the cluster expressing the gene, pct.2 = fraction in all other cells
top_markers <- markers %>%
  filter(pct.1 > 0.5, pct.2 < 0.3) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

print(top_markers[, c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2")], n = 110)

# # Save the full marker table and the top 10 per cluster
write.csv(markers, "cluster_markers_all.csv", row.names = FALSE)
write.csv(top_markers, "cluster_markers_top10.csv", row.names = FALSE)

# Dot plot of the top 3 data-derived markers per cluster
# Note: pct.2 < 0.3 leaves clusters 3 and 6 with few or no markers (they share genes with other epithelial clusters)
dot_genes <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC) %>%
  pull(gene) %>%
  unique()

DotPlot(seurat_obj, features = dot_genes, group.by = "seurat_clusters") + RotatedAxis()


### Step 6: Cell Type Annotation with SingleR (human reference)
# Compare each cell to a reference of known cell types
ref <- celldex::HumanPrimaryCellAtlasData()

annotations <- SingleR(
  test = GetAssayData(seurat_obj, layer = "data"),
  ref = ref,
  labels = ref$label.main
)

# Raw per-cell SingleR labels
seurat_obj$SingleR_label <- annotations$labels

DimPlot(seurat_obj, group.by = "SingleR_label", label = TRUE, repel = TRUE) +
  ggtitle("UMAP: SingleR Labels")

# How the SingleR labels fall across clusters (rows = clusters)
label_table <- table(seurat_obj$seurat_clusters, seurat_obj$SingleR_label)
print(label_table)

# Majority SingleR label for each cluster
apply(label_table, 1, function(x) names(x)[which.max(x)])

# Compare this majority label with the top markers from Step 5.
# Look the marker genes up in a database (e.g. CellMarker 2.0, PanglaoDB)


### Step 7: Final Cell Type Labels and Differential Expression
# 7a. Final cluster-level labels. SingleR's majority label, checked against the markers.
# Check that these still match your tables if you change any setting or the seed.
cluster_labels <- c("0" = "Epithelial", "1" = "Epithelial", "2" = "Epithelial",
                    "3" = "Epithelial", "4" = "Macrophage/Monocyte", "5" = "Epithelial",
                    "6" = "Epithelial", "7" = "Epithelial", "8" = "Fibroblast",
                    "9" = "Endothelial",
                    "10" = "Pericyte") #SingleR: Tissue_stem_cells; top markers (COX4I2, HIGD1B, NDUFA4L2, ACTA2) indicate pericytes

seurat_obj$celltype_final <- unname(cluster_labels[as.character(seurat_obj$seurat_clusters)])
Idents(seurat_obj) <- "celltype_final"

table(seurat_obj$celltype_final)

p_final <- DimPlot(seurat_obj, group.by = "celltype_final", label = TRUE, repel = TRUE) +
  ggtitle("UMAP: Final Cell Types")
p_final

dir.create("plots", showWarnings = FALSE)
ggsave("plots/umap_final_celltypes.png", plot = p_final, width = 8, height = 6, dpi = 300)

# 7b. Differential expression between two cell types
# ident.1 / ident.2 can be changed to any labels in celltype_final
de_genes <- FindMarkers(seurat_obj,
                        ident.1 = "Epithelial",
                        ident.2 = "Macrophage/Monocyte",
                        min.pct = 0.25,
                        logfc.threshold = 0.25)

# Positive avg_log2FC = higher in Epithelial; negative = higher in Macrophage/Monocyte
de_genes$gene <- rownames(de_genes)
de_genes$p_val_adj[de_genes$p_val_adj == 0] <- 1e-300   # avoid infinite values on the log scale
de_genes$significant <- ifelse(de_genes$p_val_adj < 0.05 & abs(de_genes$avg_log2FC) > 0.5,
                               "Yes", "No")

head(de_genes[order(de_genes$avg_log2FC, decreasing = TRUE), ], 10)   # top epithelial genes
head(de_genes[order(de_genes$avg_log2FC), ], 10)                      # top macrophage genes

write.csv(de_genes, "Epithelial_vs_Macrophage_DEGs.csv", row.names = FALSE)

# 7c. Volcano plot
volcano <- ggplot(de_genes, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = significant)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Epithelial vs Macrophage/Monocyte",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")
volcano

ggsave("plots/volcano_epithelial_vs_macrophage.png", plot = volcano, width = 8, height = 6, dpi = 300)


### 8.Label the genes with the 10 biggest fold changes on each side,
# i.e. genes that are significant and expressed in >50% of cells in at least one group
sig <- de_genes %>%
  filter(significant == "Yes", pmax(pct.1, pct.2) > 0.5)

top_epi <- sig %>% slice_max(n = 10, order_by = avg_log2FC)   # higher in Epithelial
top_mac <- sig %>% slice_min(n = 10, order_by = avg_log2FC)   # higher in Macrophage/Monocyte

de_genes$label <- ifelse(de_genes$gene %in% c(top_epi$gene, top_mac$gene),
                         de_genes$gene, NA)

# Label volcano plot
volcano <- ggplot(de_genes, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = significant)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_text_repel(aes(label = label),
                  colour = "black", size = 3,
                  max.overlaps = Inf, na.rm = TRUE) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Epithelial vs Macrophage/Monocyte",
       subtitle = "Right = higher in Epithelial, Left = higher in Macrophage/Monocyte",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")
volcano

ggsave("plots/volcano_labelled.png", plot = volcano, width = 9, height = 7, dpi = 300)


### 9. Enrichment Analysis ###

# Genes higher in epithelial cells (significant, log2FC > 1)
epi_genes <- de_genes %>%
  filter(significant == "Yes", avg_log2FC > 1) %>%
  pull(gene)

ego <- enrichGO(gene = epi_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
p_epi <- dotplot(ego, showCategory = 10) +
  ggtitle("GO enrichment: genes higher in epithelial cells")
p_epi

ggsave("plots/go_epithelial.png", plot = p_epi, width = 8, height = 6, dpi = 300)

# # Genes higher in macrophage cells (significant, log2FC > 1)
mac_genes <- de_genes %>%
  filter(significant == "Yes", avg_log2FC < -1) %>%
  pull(gene)

ego_mac <- enrichGO(gene = mac_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
p_mac <- dotplot(ego_mac, showCategory = 10) +
  ggtitle("GO enrichment: genes higher in Macrophage/Monocyte cells")
p_mac

ggsave("plots/go_macrophage_monocyte.png", plot = p_mac, width = 8, height = 6, dpi = 300)


# ------------------------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "session_info.txt")
