################################################
# Script Description ====
# This script uses the previously annotated full dataset and repeat the analysis only with myeloid cells
# MAIN STEPS:
# 1- Myeloid cells are extracted from the full dataset
# 2- Integration, Dimensional reduction and DEG testing is performed with Seurat
# 3- The annotated myeloid dataset is saved for further analysis

# Set Working Directory ====
dir <- "."
setwd(dir)

#-------------------------------------------------------------------------------------------------------------

# Load full dataset  and remove non-myeloid cells ====

# Load annotated full dataset 
path_seuobj_dir <- "../outputs/SeuratObjects"
ds_full <- readRDS(paste0(path_seuobj_dir, "/2-ds_annotated_full_dataset.rds"))

# Exclude non-myeloid cells 
excluded_clusters <- c(
  "B-cells 1",
  "B-cells 2",
  "B-cells 3",
  "T-cells 1",
  "T-cells 2",
  "T-cells 3",
  "T-cells 4",
  "T-cells 5",
  "T-cells 6",
  "T-cells 7",
  "Cycling Lymphocytes",
  "Low-Quality 1",
  "Low-Quality 2",
  "Low-Quality 3",
  "Doublets 1",
  "Doublets 2"
)

for (cluster in excluded_clusters) {
  print(paste0("Removing cells from cluster ", cluster))
  ds_full <- subset(ds_full, subset = General_annotation != cluster)
} # end for loop

# Set up a new SeuratObject to rerun the analysis
ds <- Seurat::CreateSeuratObject(counts = ds_full[["RNA"]]@counts, project = "ImmuBat")
ds@meta.data <- ds_full@meta.data #Get metadata from the full dataset
rm(ds_full)

#-------------------------------------------------------------------------------------------------------------

# Seurat integration of myeloid dataset ====
N <- 2000

ds_list <- Seurat::SplitObject(ds, split.by = "individual")
ds_list <- lapply(
  X = ds_list,
  FUN = function(x) {
    x <- Seurat::NormalizeData(x)
    x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = N)
  }
)

features <- Seurat::SelectIntegrationFeatures(object.list = ds_list, nfeatures = N)
ds_anchors <- Seurat::FindIntegrationAnchors(object.list = ds_list, anchor.features = features)
ds <- Seurat::IntegrateData(anchorset = ds_anchors)

ds@active.assay <- "integrated"
int.assay <- "integrated"
var.genes <- slot(ds@assays[[int.assay]], "var.features")


# Dimensional reduction and clustering ====
# Dimensional reduction
ds <- Seurat::ScaleData(ds, verbose = T)
ds <- Seurat::RunPCA(ds, npcs = 12, verbose = FALSE)
ds <- Seurat::RunUMAP(
  object =  ds,
  reduction = "pca",
  dims = 1:12,
  min.dist = 0.6,
  spread = 0.5,
  n.neighbors = 45
)

Seurat::DimPlot(ds, reduction = "umap")

# Clustering
ds <- Seurat::FindNeighbors(ds, reduction = "pca", dims = 1:12, k.param = 15) #Build KNN-graph

adj_matrix <- ds@graphs$integrated_snn #Extract and convert KNN-graph to an igraph
adj_matrix <- as(object = adj_matrix, Class = "dgCMatrix")
igraph <- igraph::graph_from_adjacency_matrix(adjmatrix = adj_matrix, weighted = TRUE)

ds@meta.data$Leiden_clusters <- as.factor(
  leidenbase::leiden_find_partition(
    igraph               = igraph,
    resolution_parameter = 0.003,
    seed                 = 29
  )$membership
)

Seurat::DimPlot(ds, reduction = "umap", group.by = "Leiden_clusters", label = T) + Seurat::NoLegend()
Seurat::DimPlot(ds, reduction = "umap", group.by = "General_annotation", label = T) + Seurat::NoLegend()

# Visualization of dataset 
Seurat::FeaturePlot(ds, reduction = "umap", features = "percent_mt", min.cutoff = 0, max.cutoff = 20)
Seurat::FeaturePlot(ds, reduction = "umap", features = "nCount_RNA", min.cutoff = 0)
Seurat::FeaturePlot(ds, reduction = "umap", features = "nFeature_RNA", min.cutoff = 0, max.cutoff = 5000)

Seurat::DimPlot(ds, reduction = "umap", group.by = "individual")
Seurat::DimPlot(ds, reduction = "umap", group.by = "status")
Seurat::DimPlot(ds, reduction = "umap", group.by = "type")
Seurat::DimPlot(ds, reduction = "umap", group.by = "age")

#-------------------------------------------------------------------------------------------------------------

# DEG testing to find marker genes ====
Seurat::Idents(ds) <- "Leiden_clusters"

deg_myelo <- Seurat::FindAllMarkers(
  ds, 
  assay = "RNA", 
  logfc.threshold = 0.5, 
  min.pct = 0.2,
  only.pos = T,
  return.thresh = 0.05,
  test.use = "wilcox"
)

# Cluster annotation ====
# Proposed annotation for all 9 myeloid clusters
clusters_annotation <-  c(
  "Neutrophils 1", #1     
  "CD16 Monocytes", #2
  "Eosinophils", #3
  "Neutrophils 2", #4
  "CD14 Monocytes", #5
  "cDC2", #6
  "cDC1", #7
  "Doublets 3", #8
  "pDC" #9
)

# Create lookup table Cluster/Annotation
list_clusters <- levels(ds$Leiden_clusters)

look_clusters <- data.frame(
  Annotation = clusters_annotation,
  Cluster = list_clusters
)

# Add Annotation as metadata
ds$Refined_annotation <- sapply(
  X = ds$Leiden_clusters,
  FUN = function(x) look_clusters$Annotation[match(x, look_clusters$Cluster)]
)

Seurat::DimPlot(ds, reduction = "umap", group.by = "Refined_annotation", label = T) + Seurat::NoLegend()

#-------------------------------------------------------------------------------------------------------------

# Save DEG testing results ====
# Add annotation to the DEG table
deg_myelo$Refined_annotation <- sapply(
  X = deg_myelo$cluster,
  FUN = function(x) look_clusters$Annotation[match(x, look_clusters$Cluster)]
)

# Set annotation order
annotation.levels <-  c(
  "Neutrophils 1",
  "Neutrophils 2",
  "Eosinophils",
  "CD16 Monocytes",
  "CD14 Monocytes",
  "cDC2",
  "cDC1",
  "pDC",
  "Doublets 3"
)


# Reorder DEG table 
deg_myelo$Refined_annotation <- factor(deg_myelo$Refined_annotation, levels = annotation.levels)
deg_myelo <- dplyr::arrange(deg_myelo, Refined_annotation) 
deg_myelo <- dplyr::select(deg_myelo, gene, cluster, Refined_annotation, p_val, p_val_adj, avg_log2FC, pct.1, pct.2) 
deg_myelo <- dplyr::filter(deg_myelo, p_val_adj <= 0.05) 

#Save table
path_deg_dir <- "../outputs/DEG_testing"
dir.create(path_deg_dir)

write.csv(deg_myelo, paste0(path_deg_dir, "/DEG_myeloid_dataset.csv"), row.names = F)

# Save Seurat Object with refined Annotation ====
saveRDS(ds, paste0(path_seuobj_dir, "/3-ds_annotated_myeloid_dataset.rds"))

# End of document 
################################################
