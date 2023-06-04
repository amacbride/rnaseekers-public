library(Seurat)
library(SeuratDisk)

## Loading the Data

# Define list of sample names
samples <- c()

# Define list of sample names for each type of matrix
samples_gene <- c( "PDAC_TISSUE_16") 
samples_feature <- c("PDAC_TISSUE_1", "PDAC_TISSUE_2", "PDAC_TISSUE_3", "PDAC_TISSUE_4", "PDAC_TISSUE_5", 
             "PDAC_TISSUE_6", "PDAC_TISSUE_7", "PDAC_TISSUE_8", "PDAC_TISSUE_9", "PDAC_TISSUE_10", 
             "PDAC_TISSUE_11B", "PDAC_TISSUE_11A", "PDAC_TISSUE_12", "PDAC_TISSUE_13", 
             "PDAC_TISSUE_15", "AdjNorm_TISSUE_1", "AdjNorm_TISSUE_2", "AdjNorm_TISSUE_3")  

# Initialize an empty list to store the Seurat objects
seurat.list <- list()

# Loop through each sample for the gene matrices
for (i in 1:length(samples_gene)) {
    data.dir <- paste0("/Users/sarrahrose/Downloads/GSE155698_RAW/", samples_gene[i], "/filtered_gene_bc_matrices")
    counts <- Read10X(data.dir = data.dir)
    seurat.obj <- CreateSeuratObject(counts = counts, project = samples_gene[i])
    seurat.list[[samples_gene[i]]] <- seurat.obj
}

# Loop through each sample for the feature matrices
for (i in 1:length(samples_feature)) {
    data.dir <- paste0("/Users/sarrahrose/Downloads/GSE155698_RAW/", samples_feature[i], "/filtered_feature_bc_matrix")
    counts <- Read10X(data.dir = data.dir)
    seurat.obj <- CreateSeuratObject(counts = counts, project = samples_feature[i])
    seurat.list[[samples_feature[i]]] <- seurat.obj
}

# Specify the path to the .h5 file
data.dir <- "/Users/sarrahrose/Downloads/GSE155698_RAW/PDAC_TISSUE_14/filtered_feature_bc_matrix.h5"

# Read the count matrix
counts <- Read10X_h5(filename = data.dir)

# Create a Seurat object
seurat.obj <- CreateSeuratObject(counts = counts, project = "PDAC_TISSUE_14")

# Add the Seurat object to the list
seurat.list[["PDAC_TISSUE_14"]] <- seurat.obj


# Adding the "cancerous" tag onto metadata
# Define the list of cancerous sample names
cancer_samples <- c("PDAC_TISSUE_1", "PDAC_TISSUE_2", "PDAC_TISSUE_3", "PDAC_TISSUE_4", "PDAC_TISSUE_5", 
             "PDAC_TISSUE_6", "PDAC_TISSUE_7", "PDAC_TISSUE_8", "PDAC_TISSUE_9", "PDAC_TISSUE_10", 
             "PDAC_TISSUE_11B", "PDAC_TISSUE_11A", "PDAC_TISSUE_12", "PDAC_TISSUE_13", 
             "PDAC_TISSUE_14", "PDAC_TISSUE_15", "PDAC_TISSUE_16") 

# Add the 'cancerous' column to each Seurat object's metadata
for (i in 1:length(seurat.list)) {
    seurat.list[[i]]@meta.data$cancerous <- ifelse(seurat.list[[i]]@meta.data$orig.ident %in% cancer_samples, "Yes", "No")
}

#Visualise metadata head
head(seurat.list[[1]]@meta.data)

# Calculate the % of mitochondrial genes for each sample and visualize the distributions
for(i in 1:length(seurat.list)) {
  seurat.list[[i]][["percent.mt"]] <- PercentageFeatureSet(seurat.list[[i]], pattern = "^MT-")
  VlnPlot(seurat.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
}

# Filter each samples individually 
for(i in 1:length(seurat.list)) {
  seurat.list[[i]] <- subset(seurat.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
}

# Normalise data
for(i in 1:length(seurat.list)) {
  seurat.list[[i]] <- NormalizeData(seurat.list[[i]])
}

# Find variable features for each sample
for(i in 1:length(seurat.list)) {
  seurat.list[[i]] <- FindVariableFeatures(seurat.list[[i]])
}

for(i in 1:length(seurat.list)) {
  seurat.list[[i]] <- ScaleData(seurat.list[[i]])
  seurat.list[[i]] <- RunPCA(seurat.list[[i]])
}

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:20)

# Integrate data
integrated <- IntegrateData(anchorset = anchors, dims = 1:20)

# Run PCA on integrated data
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated)

# Run clustering and UMAP
integrated <- FindNeighbors(integrated, dims = 1:20)
integrated <- FindClusters(integrated, resolution = 0.5)
integrated <- RunUMAP(integrated, 
                      dims = 1:20,
                      n.neighbors = 15, 
                      min.dist = 0.35, 
                      n.components = 16, 
                      metric = “cosine”,
                      force.recalc = TRUE,
                      force.name = "cancerous")
# Visualize clusters
UMAPPlot(integrated, group.by = "cancerous")

# Load Muraro et al. dataset
library(scRNAseq)
library(SingleR)

refsc <- MuraroPancreasData()
refsc <- refsc[,!is.na(refsc$label) & refsc$label!="unclear"] 

# Copy the SingleCellExperiment to a new variable before cleaning up gene names
refcl <- refsc

# Clean up gene names in the new variable
rownames(refcl) <- gsub("__chr.*", "", rownames(refcl))

# Compute log-expression values for use in marker detection
library(scater)
refcl <- logNormCounts(refcl)

# Extract assay data
ref_assay <- logcounts(refcl)

# Apply SingleR() to dataset with Muraro dataset as reference
test_assay <- integrated@assays$integrated@counts 
pred_data <- SingleR(test=test_assay, ref=ref_assay, labels=refcl$label)

#Plot assignment scores 
plotScoreHeatmap(pred_data)

#Plot delta distribution
plotScoreHeatmap(pred_data)
