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
look.table <- read.csv("../data/Library_info.csv")
#soup_cluster has been determined after running Souporcell for each library
#individual correspondence between A5 and A6 libraries has also been determined by Souporcell
#sex has been determined after checking expression of selected genes for each individual

# Add library metadata for each single cell ====
# Function to extract supplementary metadata for the cells of one library
meta <- c("type", "individual", "age", "sex")

extract.metadata <- function(ds, lookup, meta) {
  lib.lookup <- lookup[lookup$library == unique(ds$library),] #retain only information about the considered library
  extracted <- list()
  extracted[["CBC"]] <- colnames(ds)
  
  for (i in meta) {
    col <- which(colnames(lib.lookup)==i) #select col with metadata of interest
    data <- sapply(ds$soup_cluster,
                   function(x) lib.lookup[,col][match(stringr::str_sub(x, 1, 1), lib.lookup$soup_cluster)])
    extracted[[i]] <- data
  }
  return(extracted)
}

# Extract current metadata of the SeuratObject
current.data <- cbind(rownames(ds@meta.data),
                      ds@meta.data
                      )
colnames(current.data) <- c("CBC", colnames(ds@meta.data))

# Extract the supplementary metadata from library info
supp.data <- rbind(dplyr::as_tibble(extract.metadata(subset(ds, subset = library == "A1"), look.table, meta)),
                   dplyr::as_tibble(extract.metadata(subset(ds, subset = library == "A2"), look.table, meta)),
                   dplyr::as_tibble(extract.metadata(subset(ds, subset = library == "A5"), look.table, meta)),
                   dplyr::as_tibble(extract.metadata(subset(ds, subset = library == "A6"), look.table, meta))
                   )

# Create and add the full metadata table
meta.data <- dplyr::inner_join(current.data, supp.data, by = "CBC")
rownames(meta.data) <- meta.data$CBC
ds@meta.data <- meta.data

# Add percentage of mitochondrial reads for each single cell ====
ds$percent.mt <- Seurat::PercentageFeatureSet(ds, pattern = "^MT-")
#-------------------------------------------------------------------------------------------------------------

# Check bat sex based on gene expression ====
both.genes <- c("TMSB4X","DDX3X")
male.genes <- c("LOC107506326", "LOC107502043")
female.genes <- c("LOC107502967")
genes <- c(both.genes, male.genes, female.genes)
#LOC107506326 found to be a candidate for TMSB4Y
#LOC107502043 found to be a candidate for DDX3Y
#LOC107502967 found to be a candidate for XIST (X Inactive Specific Transcript)

# Create plots
plot.list <- list()
for (i in genes) {
  data <- dplyr::tibble(individual = ds$individual,
                        gene = ds[["RNA"]]@counts[which(rownames(ds)== i),],
                        Souporcell_status = ds$status)
  
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = individual, y = gene, fill = Souporcell_status)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::ylab("Total UMI count") +
    ggplot2::ggtitle(i) +
    cowplot::theme_cowplot()
  plot.list[[i]] <- plot
}

# Plotting and save
cowplot::plot_grid(plotlist = plot.list,
                   ncol = 3,
                   nrow = 2)

filename <- "../outputs/Sex_related_gene_expression.pdf"
ggplot2::ggsave(
  filename = filename,
  plot     = cowplot::plot_grid(plotlist = plot.list,
                                ncol = 2,
                                nrow = 3),
  width    = 20,
  height   = 15
)
#-------------------------------------------------------------------------------------------------------------

# Seurat integration of full dataset ====
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
ds <- Seurat::RunPCA(ds, npcs = 25, verbose = FALSE)
ds <- Seurat::RunUMAP(object =  ds,
              reduction = "pca",
              dims = 1:25,
              min.dist = 0.45,
              spread = 0.5,
              n.neighbors = 50)
Seurat::DimPlot(ds, reduction = "umap")

# Clustering
ds <- Seurat::FindNeighbors(ds, reduction = "pca", dims = 1:25, k.param = 20) #Build KNN-graph

adj.matrix <- ds@graphs$integrated_snn #Extract and convert KNN-graph to an igraph
adj.matrix <- as(object = adj.matrix, Class = "dgCMatrix")
igraph <- igraph::graph_from_adjacency_matrix(adjmatrix = adj.matrix, weighted = TRUE)

ds@meta.data$Leiden_clusters <- as.factor(leidenbase::leiden_find_partition(
  igraph               = igraph,
  resolution_parameter = 0.003,
  seed                 = 29
)$membership
)

Seurat::DimPlot(ds, reduction = "umap", group.by = "Leiden_clusters", label = T) + Seurat::NoLegend()

# Visualization dataset ====
Seurat::FeaturePlot(ds, reduction = "umap", features = "percent.mt", min.cutoff = 0, max.cutoff = 40)
Seurat::FeaturePlot(ds, reduction = "umap", features = "nCount_RNA", min.cutoff = 0)
Seurat::FeaturePlot(ds, reduction = "umap", features = "nFeature_RNA", min.cutoff = 0, max.cutoff = 5000)

Seurat::DimPlot(ds, reduction = "umap", group.by = "individual")
Seurat::DimPlot(ds, reduction = "umap", group.by = "status")
Seurat::DimPlot(ds, reduction = "umap", group.by = "type")
Seurat::DimPlot(ds, reduction = "umap", group.by = "age")
#-------------------------------------------------------------------------------------------------------------

# DEG testing to find marker genes ====
Seurat::Idents(ds) <- "Leiden_clusters"
deg.full <- Seurat::FindAllMarkers(ds, 
                                   assay = "RNA", 
                                   logfc.threshold = 0.5, 
                                   min.pct = 0.2,
                                   only.pos = T,
                                   return.thresh = 0.05,
                                   test.use = "wilcox")

# Cluster annotation ====
# Proposed annotation for all 25 clusters
clusters.annotation <-  c("B-cells 1", #1     
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
                          "DCs 2") #25

# Create lookup table Cluster/Annotation
list.clusters <- levels(ds$Leiden_clusters)
look.clusters <- data.frame(Annotation = clusters.annotation,
                            Cluster = list.clusters)

# Add Annotation as metadata
ds$General_annotation <- sapply(ds$Leiden_clusters,
                        function(x) look.clusters$Annotation[match(x, look.clusters$Cluster)])

Seurat::DimPlot(ds, reduction = "umap", group.by = "General_annotation", label = T) + Seurat::NoLegend()
#-------------------------------------------------------------------------------------------------------------

# Save DEG testing results ====
# Add annotation to the DEG table
deg.full$General_annotation <- sapply(deg.full$cluster,
                                      function(x) look.clusters$Annotation[match(x, look.clusters$Cluster)])
# Set annotation order
annotation.levels <-  c("B-cells 1",
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
                        "Doublets 2")


# Reorder DEG table 
deg.full$General_annotation <- factor(deg.full$General_annotation, levels = annotation.levels)
deg.full <- dplyr::arrange(deg.full, General_annotation) #reorder rows based on annotations
deg.full <- dplyr::select(deg.full, gene, cluster, General_annotation, p_val, p_val_adj, avg_log2FC, pct.1, pct.2) #reorder columns
deg.full <- dplyr::filter(deg.full, p_val_adj <= 0.05) #filter out adjusted p-value higher than 0.05

#Save table
path_deg_dir <- "../outputs/DEG_testing"
dir.create(path_deg_dir)

write.csv(deg.full, paste0(path_deg_dir, "/DEG_full_dataset.csv"), row.names = F)

# Save Seurat Object with General Annotation ====
saveRDS(ds, paste0(path_seuobj_dir, "/2-ds_annotated_full_dataset.rds"))

# End of document 
################################################
