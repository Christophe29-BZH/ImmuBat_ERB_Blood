################################################
# Script Description ====
# This script uses the previously annotated full dataset and repeat the analysis only with Lymphoid cells
# MAIN STEPS:
# 1- Lymphoid cells are extracted from the full dataset
# 2- Integration, Dimensional reduction and DEG testing is performed with Seurat
# 3- The whole workflow is repeated again specifically on the T cells for subclustering
# 4- The annotated Lymphoid dataset is saved for further analysis

# Set Working Directory ====
dir <- "."
setwd(dir)

#-------------------------------------------------------------------------------------------------------------
# Load full dataset  and remove non-lymphoid cells ====

# Load annotated full dataset 
path_seuobj_dir <- "../outputs/SeuratObjects"
ds.full <- readRDS(paste0(path_seuobj_dir, "/2-ds_annotated_full_dataset.rds"))

# Exclude non-myeloid cells 
excluded_clusters <- c("Granulocytes 1",
                       "Granulocytes 2",
                       "Granulocytes 3",
                       "Granulocytes 4",
                       "Monocytes 1",
                       "Monocytes 2",
                       "Monocytes 3",
                       "DCs 1",
                       "DCs 2",
                       "Low-Quality 1",
                       "Low-Quality 2",
                       "Low-Quality 3",
                       "Doublets 1",
                       "Doublets 2")

for (cluster in excluded_clusters) {
  print(paste0("Removing cells from cluster ", cluster))
  ds.full <- subset(ds.full, subset = General_annotation != cluster)
}

# Set up a new SeuratObject to rerun the analysis
ds <- Seurat::CreateSeuratObject(counts = ds.full[["RNA"]]@counts, project = "ImmuBat")
ds@meta.data <- ds.full@meta.data #Get metadata from the full dataset
rm(ds.full)
#-------------------------------------------------------------------------------------------------------------

# Seurat integration of lymphoid dataset ====
N <- 3000

ds.list <- Seurat::SplitObject(ds, split.by = "individual")
ds.list <- lapply(X = ds.list, FUN = function(x) {
  x <- Seurat::NormalizeData(x)
  x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = N)
})
features <- Seurat::SelectIntegrationFeatures(object.list = ds.list, nfeatures = N)
ds.anchors <- Seurat::FindIntegrationAnchors(object.list = ds.list, anchor.features = features)
ds <- Seurat::IntegrateData(anchorset = ds.anchors)
ds@active.assay <- "integrated"
int.assay <- "integrated"
var.genes <- slot(ds@assays[[int.assay]], "var.features")


# Dimensional reduction and clustering ====
# Dimensional reduction
ds <- Seurat::ScaleData(ds, verbose = T)
ds <- Seurat::RunPCA(ds, npcs = 18, verbose = FALSE)
ds <- Seurat::RunUMAP(object =  ds,
                      reduction = "pca",
                      dims = 1:18,
                      min.dist = 0.35,
                      spread = 0.5,
                      n.neighbors = 50)
Seurat::DimPlot(ds, reduction = "umap")

# Clustering
ds <- Seurat::FindNeighbors(ds, reduction = "pca", dims = 1:18, k.param = 20) #Build KNN-graph

adj.matrix <- ds@graphs$integrated_snn #Extract and convert KNN-graph to an igraph
adj.matrix <- as(object = adj.matrix, Class = "dgCMatrix")
igraph <- igraph::graph_from_adjacency_matrix(adjmatrix = adj.matrix, weighted = TRUE)

ds@meta.data$Leiden_clusters <- as.factor(leidenbase::leiden_find_partition(
  igraph               = igraph,
  resolution_parameter = 0.002,
  seed                 = 29
)$membership
)

Seurat::DimPlot(ds, reduction = "umap", group.by = "Leiden_clusters", label = T) + Seurat::NoLegend()
Seurat::DimPlot(ds, reduction = "umap", group.by = "General_annotation", label = T) + Seurat::NoLegend()

# Visualization dataset ====
Seurat::FeaturePlot(ds, reduction = "umap", features = "percent.mt", min.cutoff = 0, max.cutoff = 20)
Seurat::FeaturePlot(ds, reduction = "umap", features = "nCount_RNA", min.cutoff = 0, max.cutoff = 10000)
Seurat::FeaturePlot(ds, reduction = "umap", features = "nFeature_RNA", min.cutoff = 0, max.cutoff = 3000)

Seurat::DimPlot(ds, reduction = "umap", group.by = "individual")
Seurat::DimPlot(ds, reduction = "umap", group.by = "status")
Seurat::DimPlot(ds, reduction = "umap", group.by = "type")
Seurat::DimPlot(ds, reduction = "umap", group.by = "age")
#-------------------------------------------------------------------------------------------------------------

# DEG testing to find marker genes ====
Seurat::Idents(ds) <- "Leiden_clusters"
deg_lymph <- Seurat::FindAllMarkers(ds, 
                                    assay = "RNA", 
                                    logfc.threshold = 0.5, 
                                    min.pct = 0.2,
                                    only.pos = T,
                                    return.thresh = 0.05,
                                    test.use = "wilcox")

# Cluster annotation ====
# Proposed annotation for all 12 lymphoid clusters
clusters_annotation <-  c("NKT-like 1", #1     
                          "LTB/B cells", #2
                          "PLAC8/B cells", #3
                          "T cells", #4
                          "T cells", #5
                          "T cells", #6
                          "NKT-like 2", #7
                          "ZBTB32/B cells", #8
                          "VPREB3/B cells", #9
                          "Doublets 4", #10
                          "Cycling Lymphocytes", #11
                          "T cells") #12

# Create lookup table Cluster/Annotation
list_clusters <- levels(ds$Leiden_clusters)
look_clusters <- data.frame(Annotation = clusters_annotation,
                            Cluster = list_clusters)

# Add Annotation as metadata
ds$Refined_annotation <- sapply(ds$Leiden_clusters,
                                function(x) look_clusters$Annotation[match(x, look_clusters$Cluster)])

Seurat::DimPlot(ds, reduction = "umap", group.by = "Refined_annotation", label = T) + Seurat::NoLegend()
#-------------------------------------------------------------------------------------------------------------

# T cell subclustering ====
# Set up a new SeuratObject to rerun the analysis
ds.tcells <- Seurat::CreateSeuratObject(counts = ds@assays$RNA@counts, meta.data = ds@meta.data)
ds.tcells <- subset(ds.tcells, subset = Refined_annotation == "T cells")

# Seurat integration of only T cells
N <- 1500

ds.tcells.list <- Seurat::SplitObject(ds.tcells, split.by = "individual")
ds.tcells.list <- lapply(X = ds.tcells.list, FUN = function(x) {
  x <- Seurat::NormalizeData(x)
  x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = N)
})
features <- Seurat::SelectIntegrationFeatures(object.list = ds.tcells.list, nfeatures = N)
ds.tcells.anchors <- Seurat::FindIntegrationAnchors(object.list = ds.tcells.list, anchor.features = features)
ds.tcells <- Seurat::IntegrateData(anchorset = ds.tcells.anchors)
ds.tcells@active.assay <- "integrated"
int.assay <- "integrated"
var.genes <- slot(ds.tcells@assays[[int.assay]], "var.features")

# Dimensional reduction
ds.tcells <- Seurat::ScaleData(ds.tcells, verbose = T)
ds.tcells <- Seurat::RunPCA(ds.tcells, npcs = 10, verbose = FALSE)
ds.tcells <- Seurat::RunUMAP(object =  ds.tcells,
                      reduction = "pca",
                      dims = 1:10,
                      min.dist = 0.4,
                      spread = 0.5,
                      n.neighbors = 30)
Seurat::DimPlot(ds.tcells, reduction = "umap")

# Clustering
ds.tcells <- Seurat::FindNeighbors(ds.tcells, reduction = "pca", dims = 1:10, k.param = 15) #Build KNN-graph

adj.matrix <- ds.tcells@graphs$integrated_snn #Extract and convert KNN-graph to an igraph
adj.matrix <- as(object = adj.matrix, Class = "dgCMatrix")
igraph <- igraph::graph_from_adjacency_matrix(adjmatrix = adj.matrix, weighted = TRUE)

ds.tcells@meta.data$Leiden_clusters <- as.factor(leidenbase::leiden_find_partition(
  igraph               = igraph,
  resolution_parameter = 0.025,
  seed                 = 29
)$membership
)

Seurat::DimPlot(ds.tcells, reduction = "umap", group.by = "Leiden_clusters", label = T) + Seurat::NoLegend()
Seurat::DimPlot(ds.tcells, reduction = "umap", group.by = "Refined_annotation", label = T) + Seurat::NoLegend()

# Visualization Tcells 
Seurat::FeaturePlot(ds.tcells, reduction = "umap", features = "percent.mt", min.cutoff = 0, max.cutoff = 20)
Seurat::FeaturePlot(ds.tcells, reduction = "umap", features = "nCount_RNA", min.cutoff = 0, max.cutoff = 10000)
Seurat::FeaturePlot(ds.tcells, reduction = "umap", features = "nFeature_RNA", min.cutoff = 0, max.cutoff = 3000)

Seurat::DimPlot(ds.tcells, reduction = "umap", group.by = "individual")
Seurat::DimPlot(ds.tcells, reduction = "umap", group.by = "status")
Seurat::DimPlot(ds.tcells, reduction = "umap", group.by = "type")
Seurat::DimPlot(ds.tcells, reduction = "umap", group.by = "age")

# Highlight population size of each cluster
data <- dplyr::tibble(cluster = names(summary(ds.tcells$Leiden_clusters)),
                      pop.size = summary(ds.tcells$Leiden_clusters)
                      )

data$cluster <- factor(data$cluster, levels = unique(data$cluster)) #order cluster properly

ggplot2::ggplot(data, ggplot2::aes(x = cluster, y = pop.size)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::ylab("Total cell number") +
  cowplot::theme_cowplot()

# Cluster annotation 
# Proposed annotation for all 14 T cells subclusters
clusters_annotation <-  c("CD8 EM", #1     
                          "Naive T cells", #2
                          "CD4 CM1", #3
                          "DN T cells", #4
                          "Treg", #5
                          "CD4 EM", #6
                          "CD8 CM", #7
                          "CD4 CM2", #8
                          "Doublets 5", #9
                          "Outlier_Tcell_clusters", #10
                          "Outlier_Tcell_clusters", #11
                          "Outlier_Tcell_clusters", #12
                          "Outlier_Tcell_clusters", #13
                          "Outlier_Tcell_clusters") #14

# Create lookup table Cluster/Annotation
list_clusters <- levels(ds.tcells$Leiden_clusters)
look_clusters <- data.frame(Annotation = clusters_annotation,
                            Cluster = list_clusters)

# Add Annotation as metadata
ds.tcells$Refined_annotation <- sapply(ds.tcells$Leiden_clusters,
                                function(x) look_clusters$Annotation[match(x, look_clusters$Cluster)])

Seurat::DimPlot(ds.tcells, reduction = "umap", group.by = "Refined_annotation", label = T) + Seurat::NoLegend()

#-------------------------------------------------------------------------------------------------------------

# Include results of Tcell subclustering into lymphoid dataset ====
# Extract labels and prepare lookup table
tcells_annotation <- dplyr::tibble(CBC = colnames(ds.tcells),
                                   annotation = ds.tcells$Refined_annotation)

lymph_annotation <- dplyr::tibble(CBC = colnames(ds),
                                  annotation = ds$Refined_annotation)

lymph_annotation <- lymph_annotation[-which(lymph_annotation$annotation == "T cells"), ]

clusters_annotation <- rbind(tcells_annotation, lymph_annotation)

#Add to lymphoid dataset
ds$Refined_annotation <- sapply(colnames(ds),
                                function(x) clusters_annotation$annotation[match(x, clusters_annotation$CBC)])

Seurat::DimPlot(ds, reduction = "umap", group.by = "Refined_annotation", label = T)
#-------------------------------------------------------------------------------------------------------------

# DEG testing to find marker genes ====
Seurat::Idents(ds) <- "Refined_annotation"
deg_lymph <- Seurat::FindAllMarkers(ds, 
                                    assay = "RNA", 
                                    logfc.threshold = 0.5, 
                                    min.pct = 0.2,
                                    only.pos = T,
                                    return.thresh = 0.05,
                                    test.use = "wilcox")


# Save DEG testing results ====
# Set annotation order
annotation_levels <-  c("NKT-like 2",
                        "NKT-like 1",
                        "CD8 EM",
                        "CD8 CM",
                        "CD4 EM",
                        "CD4 CM1",
                        "CD4 CM2",
                        "Treg",
                        "DN T cells",
                        "Naive T cells",
                        "ZBTB32/B cells",
                        "PLAC8/B cells",
                        "LTB/B cells",
                        "VPREB3/B cells",
                        "Cycling Lymphocytes",
                        "Doublets 4",
                        "Doublets 5",
                        "Outlier_Tcell_clusters")

# Reorder DEG table 
deg_lymph$cluster <- factor(deg_lymph$cluster, levels = annotation_levels)
deg_lymph <- dplyr::arrange(deg_lymph, cluster) #reorder rows based on annotations
deg_lymph <- dplyr::select(deg_lymph, gene, cluster, p_val, p_val_adj, avg_log2FC, pct.1, pct.2) #reorder columns
deg_lymph <- dplyr::filter(deg_lymph, p_val_adj <= 0.05) #filter out adjusted p-value higher than 0.05

#Save table
path_deg_dir <- "../outputs/DEG_testing"
dir.create(path_deg_dir)

write.csv(deg_lymph, paste0(path_deg_dir, "/DEG_lymphoid_dataset.csv"), row.names = F)

# Save Seurat Object with General Annotation ====
saveRDS(ds, paste0(path_seuobj_dir, "/4-ds_annotated_lymphoid_dataset.rds"))

# End of document 
################################################