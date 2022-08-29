################################################
# Script Description ====
# This script uses the downloaded and merged libraries to perform integration and explore the full dataset
# MAIN STEPS:
# 1- All metadata regarding the different libraries and individuals are added
# 2- Integration, Dimensional reduction and DEG testing is performed with Seurat
# 3- The annotated dataset is saved for further analysis

# Set Working Directory ====
dir <- "."
setwd(dir)

#-------------------------------------------------------------------------------------------------------------

# Load merged dataset ====
path_seuobj_dir <- "../outputs/SeuratObjects"
ds <- readRDS(paste0(path_seuobj_dir, "/1-ds_merged_libraries.rds"))

# Load library metadata table ====
lookup_table <- read.csv("../data/Library_info.csv")
#soup_cluster has been determined after running Souporcell for each library
#individual correspondence between A5 and A6 libraries has also been determined by Souporcell
#sex has been determined after checking expression of selected genes for each individual ( below in the script )

# Add library metadata for each single cell ====
# Function to extract supplementary metadata for the cells of one library
extract_metadata <- function(ds, lib, lookup, meta) {
  lookup <- lookup[lookup$library == lib, ] #retain only information about the considered library
  ds <- subset(ds, subset = library == lib)
  
  extracted <- list()
  extracted[["CBC"]] <- colnames(ds)
  
  for (info in meta) {
    col <- which(colnames(lookup) == info) #select col with metadata of interest
    data <- sapply(
      X = ds$soup_cluster,
      FUN = function(x) lookup[, col][match(stringr::str_sub(x, start = 1, end = 1), lookup$soup_cluster)]
    ) #str_sub is used for the case of doublets, to get the most probable individual
    extracted[[info]] <- data
  } # end for loop
  
  return(extracted)
}

# Extract current metadata of the SeuratObject
current_data <- cbind(
  rownames(ds@meta.data),
  ds@meta.data
  
)
colnames(current_data) <- c("CBC", colnames(ds@meta.data))

# Extract the supplementary metadata from library info
meta <- c("type", "individual", "age", "sex")

supp_data <- rbind(
  dplyr::as_tibble(extract_metadata(ds, "A1", lookup_table, meta)),
  dplyr::as_tibble(extract_metadata(ds, "A2", lookup_table, meta)),
  dplyr::as_tibble(extract_metadata(ds, "A5", lookup_table, meta)),
  dplyr::as_tibble(extract_metadata(ds, "A6", lookup_table, meta))
)

# Create and add the full metadata table
full_metadata <- dplyr::inner_join(current_data, supp_data, by = "CBC")
rownames(full_metadata) <- full_metadata$CBC
ds@meta.data <- full_metadata

# Add percentage of mitochondrial reads for each single cell ====
ds$percent_mt <- Seurat::PercentageFeatureSet(ds, pattern = "^MT-")

#-------------------------------------------------------------------------------------------------------------

# Check bat sex based on gene expression ====
both_genes <- c("TMSB4X","DDX3X")
male_genes <- c("LOC107506326", "LOC107502043")
female_genes <- c("LOC107502967")
genes <- c(both_genes, male_genes, female_genes)
#LOC107506326 found to be a candidate for TMSB4Y
#LOC107502043 found to be a candidate for DDX3Y
#LOC107502967 found to be a candidate for XIST (X Inactive Specific Transcript)

# Create plots
plot_list <- list()

for (i in genes) {
  data <- dplyr::tibble(
    individual = ds$individual,
    gene = ds[["RNA"]]@counts[which(rownames(ds)== i),],
    Souporcell_status = ds$status
  )
  
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = individual, y = gene, fill = Souporcell_status)) +
    ggplot2::geom_bar(stat = "identity") +
    cowplot::theme_cowplot() +
    ggplot2::ylab("Total UMI count") +
    ggplot2::ggtitle(i) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  plot_list[[i]] <- plot
} #end for loop

# Plotting and save
cowplot::plot_grid(
  plotlist = plot_list,
  labels = c("Both", "Both", "Male", "Male", "Female"),
  ncol = 2,
  nrow = 3
)

ggplot2::ggsave(
  filename = "../outputs/Sex_related_gene_expression.pdf",
  width    = 20,
  height   = 15
)

#-------------------------------------------------------------------------------------------------------------

# Seurat integration of full dataset ====
N <- 3000

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
ds <- Seurat::RunPCA(ds, npcs = 25, verbose = FALSE)
ds <- Seurat::RunUMAP(
  object =  ds,
  reduction = "pca",
  dims = 1:25,
  min.dist = 0.45,
  spread = 0.5,
  n.neighbors = 50
)

Seurat::DimPlot(ds, reduction = "umap")

# Clustering
ds <- Seurat::FindNeighbors(ds, reduction = "pca", dims = 1:25, k.param = 20) #Build KNN-graph

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

# Visualization of dataset 
Seurat::FeaturePlot(ds, reduction = "umap", features = "percent_mt", min.cutoff = 0, max.cutoff = 40)
Seurat::FeaturePlot(ds, reduction = "umap", features = "nCount_RNA", min.cutoff = 0)
Seurat::FeaturePlot(ds, reduction = "umap", features = "nFeature_RNA", min.cutoff = 0, max.cutoff = 5000)

Seurat::DimPlot(ds, reduction = "umap", group.by = "individual")
Seurat::DimPlot(ds, reduction = "umap", group.by = "status")
Seurat::DimPlot(ds, reduction = "umap", group.by = "type")
Seurat::DimPlot(ds, reduction = "umap", group.by = "age")

#-------------------------------------------------------------------------------------------------------------

# DEG testing to find marker genes ====
Seurat::Idents(ds) <- "Leiden_clusters"

deg_full <- Seurat::FindAllMarkers(
  ds, 
  assay = "RNA", 
  logfc.threshold = 0.5, 
  min.pct = 0.2,
  only.pos = T,
  return.thresh = 0.05,
  test.use = "wilcox"
)

# Cluster annotation ====
# Proposed annotation for all 25 clusters
clusters_annotation <-  c(
  "B-cells 1", #1     
  "T-cells 1", #2
  "T-cells 2", #3
  "B-cells 2", #4
  "Granulocytes 1", #5
  "Monocytes 1", #6
  "Granulocytes 2", #7
  "Granulocytes 3", #8
  "T-cells 3", #9
  "B-cells 3", #10
  "T-cells 4", #11
  "T-cells 5", #12
  "Low-Quality 1", #13
  "Granulocytes 4", #14
  "T-cells 6", #15
  "Monocytes 2", #16
  "Monocytes 3", #17
  "Low-Quality 2", #18
  "Doublets 1", #19
  "Low-Quality 3", #20
  "DCs 1", #21
  "Cycling Lymphocytes", #22
  "Doublets 2", #23
  "T-cells 7", #24
  "DCs 2" #25
) 

# Create lookup table Cluster/Annotation
list_clusters <- levels(ds$Leiden_clusters)

look_clusters <- data.frame(
  Annotation = clusters_annotation,
  Cluster = list_clusters
)

# Add Annotation as metadata
ds$General_annotation <- sapply(
  X = ds$Leiden_clusters,
  FUN = function(x) look_clusters$Annotation[match(x, look_clusters$Cluster)]
)

Seurat::DimPlot(ds, reduction = "umap", group.by = "General_annotation", label = T) + Seurat::NoLegend()

#-------------------------------------------------------------------------------------------------------------

# Save DEG testing results ====
# Add annotation to the DEG table
deg_full$General_annotation <- sapply(
  X = deg_full$cluster,
  FUN = function(x) look_clusters$Annotation[match(x, look_clusters$Cluster)]
)
# Set annotation order
annotation_levels <-  c(
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
  "Granulocytes 1",
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
  "Doublets 2"
)


# Reorder DEG table 
deg_full$General_annotation <- factor(deg_full$General_annotation, levels = annotation_levels)
deg_full <- dplyr::arrange(deg_full, General_annotation) 
deg_full <- dplyr::select(deg_full, gene, cluster, General_annotation, p_val, p_val_adj, avg_log2FC, pct.1, pct.2) 
deg_full <- dplyr::filter(deg_full, p_val_adj <= 0.05) 

# Save table
path_deg_dir <- "../outputs/DEG_testing"
dir.create(path_deg_dir)

write.csv(deg_full, paste0(path_deg_dir, "/DEG_full_dataset.csv"), row.names = F)

# Save Seurat Object with General Annotation ====
saveRDS(ds, paste0(path_seuobj_dir, "/2-ds_annotated_full_dataset.rds"))

# End of document 
################################################