################################################
# Script Description ====
# This script produces the plots used in Figure 1 of the manuscript

# Set Working Directory ====
dir <- "."
setwd(dir)

# Set Output Directory ====
path_plots_dir <- "../outputs/Plots/Fig1"
dir.create(path_plots_dir, recursive = T)

#-------------------------------------------------------------------------------------------------------------

# Fig 1A: UMAP full dataset ====

# Load data
path_plots_dir <- "../outputs/Plots/Fig1"
path_seuobj_dir <- "../outputs/SeuratObjects"

ds <- readRDS(paste0(path_seuobj_dir, "/2-ds_annotated_full_dataset.rds"))
Seurat::DimPlot(ds, reduction = "umap", group.by = "General_annotation", label = T) + Seurat::NoLegend()

# Order clusters 
ds$General_annotation <- factor(
  x = ds$General_annotation, 
  levels = c(
    "Granulocytes 1",
    "Granulocytes 2",
    "Granulocytes 3",
    "Granulocytes 4",
    "Monocytes 1",
    "Monocytes 2",
    "Monocytes 3",
    "DCs 1",
    "DCs 2",
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
)
  
# Build Color palette 
colors <- c(
  "#B30000", # Granulocytes 1
  "#B30000", # Granulocytes 2
  "#B30000", # Granulocytes 3
  "#B30000", # Granulocytes 4
  "#CC6600", # Monocytes 1
  "#CC6600", # Monocytes 2
  "#CC6600", # Monocytes 3
  "#E6E600", # DCs 1
  "#E6E600", # DCs 2
  "#99BBFF", # B-cells 1
  "#99BBFF", # B-cells 2
  "#99BBFF", # B-cells 3
  "#009933", # T-cells 1
  "#009933", # T-cells 2
  "#009933", # T-cells 3
  "#009933", # T-cells 4
  "#009933", # T-cells 5
  "#009933", # T-cells 6
  "#009933", # T-cells 7
  "#660066", # Cycling Lymphocytes
  "#8C8C8C", # Low-Quality 1
  "#8C8C8C", # Low-Quality 2
  "#8C8C8C", # Low-Quality 3
  "#CCCCCC", # Doublets 1
  "#CCCCCC" # Doublets 2 
) 

color_palette <- dplyr::tibble(
  General_annotation = levels(ds$General_annotation),
  color = colors
)

# Prepare data for plotting 
data <- dplyr::as_tibble(
  ds@reductions$umap@cell.embeddings
)
colnames(data) <- c("x", "y")

data$General_annotation <- ds$General_annotation
data$color <- sapply(
  X = data$General_annotation, 
  FUN = function (x) color_palette$color[match(x, color_palette$General_annotation)]
)

# Create UMAP plot
arrow <- grid::arrow(
  angle  = 30, 
  length = ggplot2::unit(0.25, "inches"), 
  type   = "closed"
)

ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x   = x,
    y   = y
  )
) +
  ggplot2::geom_point(size  = 0.2, col = data$color) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "") +
  ggplot2::coord_fixed() +
  ggplot2::annotate(
    geom  = "text", 
    x     = c( -8, -11, -6, 5, 0, 5, -2, -4.48),
    y     = c( -8, 3, -4, 8, 2, 7, 0, 3.66),
    label = c(
      "Granulocytes", 
      "Monocytes", 
      "Dendritic cells", 
      "B cells", 
      "T cells", 
      "Cycling Lymphocytes",
      "Low quality", 
      "Doublets"
    ),
    color = unique(color_palette$color), 
    size = 5
  ) +
  ggplot2::labs(x = "UMAP-1", y = "UMAP-2") +
  ggplot2::theme(axis.title.x = ggplot2::element_text(vjust = 0.1, hjust = 0.1, size = 20),
                 axis.title.y = ggplot2::element_text(vjust = 0.1, hjust = 0.1, size = 20, angle = 90)
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(angle  = 30, length = ggplot2::unit(0.25, "inches"), type   = "closed")
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(angle  = 30, length = ggplot2::unit(0.25, "inches"),type   = "closed")
  )

# Save Plot 
ggplot2::ggsave(
  filename = paste0(path_plots_dir, "/Fig1A_UMAP_full_dataset.eps"),
  width    = 6.5,
  height   = 7
)

#-------------------------------------------------------------------------------------------------------------

# Fig 1B: Grid of feature plots main markers ====
# Load data
path_plots_dir <- "../outputs/Plots/Fig1"
path_seuobj_dir <- "../outputs/SeuratObjects"

ds <- readRDS(paste0(path_seuobj_dir, "/2-ds_annotated_full_dataset.rds"))

# Select marker genes
Feature_list <- c("MS4A1", "CD19", "CD3D",
                  "CD3E", "CD3G", "ITGAM",
                  "SIRPA", "CD68", "CSF1R",
                  "CSF2RB", "CSF3R", "IRF8")

# Generate the grid of plots 
plot_list <- list() #initialize grid of plots

for (gene in Feature_list) {
  
  data <- dplyr::as_tibble(
    ds@reductions$umap@cell.embeddings
  )
  colnames(data) <- c("x", "y")
  
  data$Feature <- ds@assays$RNA@data[gene,]
  data <- dplyr::arrange(data, Feature) # sort by ascending expression order for visualization purpose
  
  
  arrow <- grid::arrow(
    angle  = 30, length = ggplot2::unit(0.25, "inches"), type   = "closed"
  )
  
  plot <- ggplot2::ggplot(
    data    = data,
    mapping = ggplot2::aes(
      x   = x,
      y   = y,
      col = Feature
    )
  ) +
    ggplot2::geom_point(size  = 0.05) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "",
      plot.title = ggplot2::element_text(family = "Helvetica", face = "bold", size = (15), hjust = 0.5)
    ) +
    ggplot2::scale_colour_viridis_c(option = "magma", direction = -1) +
    ggplot2::coord_fixed() +
    ggplot2::ggtitle(gene)
  
  plot_list[[gene]] <- plot
  
} #end for loop

# Plot 
cowplot::plot_grid(
  plotlist = plot_list,
  nrow = 4, 
  ncol = 3, 
  align = "hv"
)

# Save grid of plots 
ggplot2::ggsave(
  filename = paste0(path_plots_dir, "/Fig1B_FeaturePlot_main_markers.png"),
  width    = 14,
  height   = 14
)

#-------------------------------------------------------------------------------------------------------------

# Fig 1C: UMAP Myeloid compartment ====
# Load data 
path_plots_dir <- "../outputs/Plots/Fig1"
path_seuobj_dir <- "../outputs/SeuratObjects"

ds <- readRDS(paste0(path_seuobj_dir, "/3-ds_annotated_myeloid_dataset.rds"))

# Order Clusters 
ds$Refined_annotation <- factor(
  x = ds$Refined_annotation, 
  levels = c("Neutrophils 1",
             "Neutrophils 2",
             "Eosinophils",
             "CD16 Monocytes",
             "CD14 Monocytes",
             "cDC2",
             "cDC1",
             "pDC",
             "Doublets 3"
  )
)

# Define and build Color palette 
colors <- c(
  "#E60000", #Neutrophils 1
  "#FF4D4D", #Neutrophils 2
  "#B30000", #Eosinophils 
  "#CC6600", #CD16 Monocytes
  "#FFBF80", #CD14 Monocytes
  "#FF8000", #cDC2
  "#E6E600", #cDC1
  "#E60073", #pDC
  "#999999" # Doublets 3 
) 

color_palette <- dplyr::tibble(
  Refined_annotation = levels(ds$Refined_annotation),
  color = colors
)

# Prepare data for plotting
data <- dplyr::as_tibble(
  ds@reductions$umap@cell.embeddings
)
colnames(data) <- c("x", "y")

data$Refined_annotation <- ds$Refined_annotation
data$color <- sapply(
  X = data$Refined_annotation, 
  FUN = function (x) color_palette$color[match(x, color_palette$Refined_annotation)]
)

# Create UMAP plot
arrow <- grid::arrow(
  angle  = 30, 
  length = ggplot2::unit(0.25, "inches"), 
  type   = "closed"
)


ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x   = x,
    y   = y
  )
) +
  ggplot2::geom_point(size  = 0.8, col = data$color) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "") +
  ggplot2::coord_fixed() +
  ggplot2::annotate(
    geom  = "text", 
    x     = c( -6, 0, -3, 7, -7, 5, 5.5, 1, 2),
    y     = c( -4, -3, -8, -2, 7, 6, 9, 10, 1),
    label = c(
      "Neutrophils 1", 
      "Neutrophils 2", 
      "Eosinophils", 
      "CD16 Monocytes", 
      "CD14 Monocytes", 
      "cDC2",
      "cDC1", 
      "pDC", 
      "Doublets"
    ),
    color = unique(color_palette$color), 
    size = 5
  ) +
  ggplot2::labs(x = "UMAP-1", y = "UMAP-2") +
  ggplot2::theme(axis.title.x = ggplot2::element_text(vjust = 0.1, hjust = 0.1, size = 20),
                 axis.title.y = ggplot2::element_text(vjust = 0.1, hjust = 0.1, size = 20, angle = 90)
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(angle  = 30, length = ggplot2::unit(0.25, "inches"), type   = "closed")
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(angle  = 30, length = ggplot2::unit(0.25, "inches"), type   = "closed")
  )

# Save plot
ggplot2::ggsave(
  filename = paste0(path_plots_dir, "/Fig1C_UMAP_myeloid_dataset.eps"),
  width    = 6.5,
  height   = 7
)

#-------------------------------------------------------------------------------------------------------------

# Fig 1D: UMAP Lymphoid compartment ====
# Load data 
path_plots_dir <- "../outputs/Plots/Fig1"
path_seuobj_dir <- "../outputs/SeuratObjects"

ds <- readRDS(paste0(path_seuobj_dir, "/4-ds_annotated_lymphoid_dataset.rds"))

# Remove Outlier clusters identified during Tcell subclustering
ds <- subset(ds, subset = Refined_annotation != "Outlier_clusters")

# Order Clusters 
ds$Refined_annotation <- factor(
  x = ds$Refined_annotation, 
  levels = c(
    "NKT-like 2",
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
    "Doublets 5"
  )
)
# Define and build Color palette 
colors <- c(
  "#CC9966", # NKT-like 2
  "#86592D", # NKT-like 1
  "#888844", # CD8 EM
  "#BBBB77", # CD8 CM
  "#006622", # CD4 EM
  "#009933", # CD4 CM1
  "#53C68C", # CD4 CM2
  "#00E600", # Treg
  "#B3E6CC", # DN T cells
  "#66FF99", # Naive T cells
  "#002B80", # ZBTB32/B cells
  "#004DE6", # PLAC8/B cells
  "#4D88FF", # LTB/B cells
  "#99BBFF", # VPREB3/B cells
  "#660066", # Cycling Lymphocytes
  "#999999", # Doublets 4
  "#999999" # Doublets 5
)  

color_palette <- dplyr::tibble(
  Refined_annotation = levels(ds$Refined_annotation),
  color = colors
)

#Prepare data for UMAP 
data <- dplyr::as_tibble(
  ds@reductions$umap@cell.embeddings
)
colnames(data) <- c("x", "y")

data$Refined_annotation <- ds$Refined_annotation
data$color <- sapply(
  X = data$Refined_annotation, 
  FUN = function (x) color_palette$color[match(x, color_palette$Refined_annotation)]
)

#Create UMAP 
arrow <- grid::arrow(
  angle  = 30, 
  length = ggplot2::unit(0.25, "inches"), 
  type   = "closed"
)

ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x   = x,
    y   = y
  )
) +
  ggplot2::geom_point(size  = 0.8, col = data$color) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "") +
  ggplot2::coord_fixed() +
  ggplot2::annotate(
    geom  = "text", 
    x     = c( 1, 0, 0, 2, 5.5, 5, 1.5, 4.5, 1, 0.5, -7, -3.5, -8.5, -5, -0.5, -1),
    y     = c( 1, 3, 4, -1.5, -2, -3, -3.15, -5, -4.75, -6, 3, 3.25, -2, -3.5, 0, -1),
    label = c("NKT-like 2", "NKT-like 1", "CD8 EM", "CD8 CM", "CD4 EM", "CD4 CM1",
              "CD4 CM2", "Treg", "DN T cells", "Naive T cells", "ZBTB32/B cells",
              "PLAC8/B cells", "LTB/B cells", "VPREB3/B cells", "Cycling Lymphocytes",
              "Doublets"),
    color = unique(color_palette$color), 
    size = 5
  ) +
  ggplot2::labs(x = "UMAP-1", y = "UMAP-2") +
  ggplot2::theme(axis.title.x = ggplot2::element_text(vjust = 0.1, hjust = 0.1, size = 20),
                 axis.title.y = ggplot2::element_text(angle = 90, vjust = 0.1, hjust = 0.1, size = 20)
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(angle  = 30, length = ggplot2::unit(0.25, "inches"), type   = "closed")
  ) +
  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(angle  = 30, length = ggplot2::unit(0.25, "inches"), type   = "closed")
  )

# Print and save Plot 
ggplot2::ggsave(
  filename = paste0(path_plots_dir, "/Fig1D_UMAP_lymphoid_dataset.eps"),
  width    = 6.5,
  height   = 7
)

#-------------------------------------------------------------------------------------------------------------

# Fig 1E: Dotplot markers genes ====
# Load data
path_plots_dir <- "../outputs/Plots/Fig1"
path_seuobj_dir <- "../outputs/SeuratObjects"

ds <- readRDS(paste0(path_seuobj_dir, "/2-ds_annotated_full_dataset.rds"))

# Load annotations from myeloid and lymphoid cells 
extract_cell_annotation <- function(seuobj, path_dir) {
  ds <- readRDS(paste0(path_dir, "/", seuobj))
  labels <- dplyr::tibble(CBC = colnames(ds),
                          Annotation = ds$Refined_annotation)
  return(labels)
}

myeloid_labels <- extract_cell_annotation("3-ds_annotated_myeloid_dataset.rds", path_seuobj_dir)
lymphoid_labels <- extract_cell_annotation("4-ds_annotated_lymphoid_dataset.rds", path_seuobj_dir)

# Reshape as a single table
lookup <- rbind(myeloid_labels, lymphoid_labels)

# Combine together the 2 Neutrophil clusters
lookup$Annotation[grep("^Neutro", lookup$Annotation)] <- "Neutrophils"

# Add refined annotation for each single cell and exclude non-relevant cells
ds$Refined_annotation <- sapply(
  X = ds$CBC, 
  FUN = function(x) lookup$Annotation[match(x, lookup$CBC)]
)

excluded_clusters <- c(
  "Low-Quality 1",
  "Low-Quality 2",
  "Low-Quality 3",
  "Doublets 1",
  "Doublets 2",
  "Doublets 3",
  "Doublets 4",
  "Doublets 5",
  "Outlier_Tcell_clusters"
)

for (cluster in excluded_clusters) {
  print(paste0("Removing cells from cluster ", cluster))
  ds <- subset(ds, subset = Refined_annotation != cluster)
} # end for loop


# Select marker genes 
genes <- c("CSF3R", "CCRL2", "G0S2", #Neutrophils
           "LOC107516222", "CCR3", "MMP1", #Eosinophils
           "LOC107502476", "CX3CR1", #CD16 Mono
           "MAFB", "CD14", "F13A1", #CD14 Mono
           "IL18", "TCEA3", #cDC2
           "CLEC9A", #cDC1
           "IRF8","TCF4", #pDC
           "CD8A", "LOC107519796", #CD4 + CD8
           "CD160", "LOC107512753", "NCR1", #NKT-like 2
           "LOC107505367", "GNLY", "PRF1",  #NKT-like 1
           "TBX21", "CXCR6", "ZNF683", #CD8 + CD4 TEM
           "GZMM", "IL18R1", #CD8 CM
           "CXCR3", #
           "CD40LG", "STAMBPL1", #CD4 TCM1
           "BATF3", "CCR6", "CCR10", #CD4 TCM2
           "CMYA5", "FOXP3", #Treg
           "LIMA1", "ID4", "PI3", #T-DN
           "FST", "RFLNB", "TCF7", #Naive Tcells
           "ZBTB32", "BHLHE41", #Bcells 4 ZBTB32
           "PLAC8", "FAM149A", "DPF1", #Bcells 1 PLAC8
           "LTB", "LOC107502080", "CR2", #Bcells 2 LTB
           "CD72", "AKAP12", "VPREB3", #Bcells 3 VPREB3
           "PCLAF", "CLSPN", "KIF15") #Cycling Lymphocytes

# Prepare data for Dotplot
ids <- rownames(ds@assays$RNA@data) %in% genes
data <- t(as.matrix(ds@assays$RNA@data[ids, ]))
data <- dplyr::as_tibble(data)
data <- dplyr::select(data, genes) #reorder columns
data$Refined_annotation <- ds@meta.data$Refined_annotation

data$Refined_annotation <- factor(
  x = data$Refined_annotation, 
  levels = c("Neutrophils",
             "Eosinophils",
             "CD16 Monocytes",
             "CD14 Monocytes",
             "cDC2",
             "cDC1",
             "pDC",
             "NKT-like 2",
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
             "Cycling Lymphocytes")
)

data <- tidyr::gather(data, "Gene", "Expression", -Refined_annotation)

data$Gene <- factor(
  x      = data$Gene,
  levels = unique(data$Gene)
)

data$Pct <- data$Expression > 0
data$Freq <- rep(1, length(data$Pct))

data <- dplyr::group_by(data, Gene)
data <- dplyr::mutate(data, Scaled = scale(as.numeric(Expression))[, 1])
data <- dplyr::group_by(data, Refined_annotation, Gene)

data <- dplyr::summarise(
  data,
  Mean   = mean(as.numeric(Expression)),
  Scaled = mean(Scaled),
  Pct    = sum(Pct)/sum(Freq)*100
)

data$Scaled[data$Scaled > 1] <- 1
data$Mean[data$Mean > 2] <- 2

# Dotplot
colorscale <- RColorBrewer::brewer.pal(n = 9, name = "RdBu")

ggplot2::ggplot(
  data = data,
  mapping = ggplot2::aes(
    y = Refined_annotation,
    x = Gene,
    size = Pct,
    col = Scaled
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_color_gradient2(midpoint=0, 
                                 low=colorscale[length(colorscale)], 
                                 mid="white",
                                 high=colorscale[1]
  ) +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(
    axis.title       = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y      = ggplot2::element_text(
      face  = "italic"
    ),
    legend.position  = "right",
    legend.title     = ggplot2::element_blank()
  ) +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(
      barheight = 11, frame.colour = "black", ticks = FALSE, order = 1
    ),
    size  = ggplot2::guide_legend(label.position = "right")
  ) +
  ggplot2::scale_size_area(max_size = 8)

# Save Dotplot
ggplot2::ggsave(
  filename = paste0(path_plots_dir, "/Fig1E_Dotplot.eps"),
  width    = 16,
  height   = 7
)

# End of document 
################################################
