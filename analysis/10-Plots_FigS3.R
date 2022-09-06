################################################
# Script Description ====
# This script produces the plots used in Supplementary Figure 3 of the manuscript

# Set Working Directory ====
dir <- "."
setwd(dir)

# Set Output Directory ====
path_plots_dir <- "../outputs/Plots/FigS3"
dir.create(path_plots_dir, recursive = T)

#-------------------------------------------------------------------------------------------------------------

# Fig S3A: Additional marker genes for myeloid clusters ====

# Load data
path_plots_dir <- "../outputs/Plots/FigS3"
path_seuobj_dir <- "../outputs/SeuratObjects"

ds <- readRDS(paste0(path_seuobj_dir, "/3-ds_annotated_myeloid_dataset.rds"))

# Select marker genes 
Feature_list <- c(
  "PTGS2", "LTF", "MMP9", "GAPT",
  "OLIG1", "LTC4S", "GATA1", "CEBPE",
  "LOC107518414", "IRF7", "LOC107507331", "LOC107511999"
)
#LOC107518414 -> PRG3
#LOC107507331 -> DQA1
#LOC107511999 -> DOB


# Extract expression data
plot_list <- list() #initialize plot list

for (gene in Feature_list) {
  data <- dplyr::as_tibble(
    ds@reductions$umap@cell.embeddings
  )
  colnames(data) <- c("x", "y")
  
  data$Feature <- ds@assays$RNA@data[gene,]
  data <- dplyr::arrange(data, Feature) #sort ascending order, for plotting order
  
  
  arrow <- grid::arrow(
    angle  = 30, length = ggplot2::unit(0.25, "inches"), type   = "closed"
  )
  
  plot <- ggplot2::ggplot(
    data    = data,
    mapping = ggplot2::aes(
      x     = x,
      y     = y,
      col   = Feature
    )
  ) +
    ggplot2::geom_point(size  = 0.4) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(
        family = "Helvetica",
        face = "bold", 
        size = (15),
        hjust = 0.5
      )
    ) +
    ggplot2::scale_colour_viridis_c(option = "magma", direction = -1) +
    ggplot2::coord_fixed() +
    ggplot2::ggtitle(gene)
  
  plot_list[[gene]] <- plot
} #end for loop

# Create grid of expression plots
plot_grid <- cowplot::plot_grid(
  plotlist = plot_list,
  nrow = 3, 
  ncol = 4, 
  align = "hv"
)


# Extract legend
leg <- ggpubr::get_legend(
  plot_list[[1]] + ggplot2::theme(
    legend.position = "right", 
    legend.text = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank()
  )
)

# Add legend to plot_grid
cowplot::plot_grid(
  plot_grid,
  leg,
  nrow = 1,
  ncol = 2,
  rel_widths = c(3,1),
  align = "hv"
)


# Save plots 
ggplot2::ggsave(
  filename = paste0(path_plots_dir, "/FigS3A_FeaturePlots_myeloid_marker_genes.eps"),
  width    = 14,
  height   = 14
)


#-------------------------------------------------------------------------------------------------------------

# Fig S3B: Additional marker genes for lymphoid clusters ====

# Load data
path_plots_dir <- "../outputs/Plots/FigS3"
path_seuobj_dir <- "../outputs/SeuratObjects"

ds <- readRDS(paste0(path_seuobj_dir, "/4-ds_annotated_lymphoid_dataset.rds"))

# Select gene markers
Feature_list <- c(
  "LEF1", "CXCR6", "NKG7", "LOC107513518",
  "LOC107513519", "CCL5", "LOC107503311", "LOC107506475",
  "IFNG", "FCER1G", "KLRB1", "LOC107516556",
  "LOC107505617", "PRF1", "GZMM", "S100A4",
  "S100A10", "S100A11", "CD99", "CD44"
)
# LOC107513518 -> GZMB
# LOC107513519 -> GZMB
# LOC107503311 -> CCL4
# LOC107506475 -> CCL3
# LOC107516556 -> NKG2
# LOC107505617 -> NKG2

# Extract expression data
plot_list <- list() #initialize plot list

for (gene in Feature_list) {
  data <- dplyr::as_tibble(
    ds@reductions$umap@cell.embeddings
  )
  colnames(data) <- c("x", "y")
  
  data$Feature <- ds@assays$RNA@data[gene,]
  data <- dplyr::arrange(data, Feature) #sort ascending order, for plotting order
  
  
  arrow <- grid::arrow(
    angle  = 30, length = ggplot2::unit(0.25, "inches"), type   = "closed"
  )
  
  plot <- ggplot2::ggplot(
    data    = data,
    mapping = ggplot2::aes(
      x     = x,
      y     = y,
      col   = Feature
    )
  ) +
    ggplot2::geom_point(size  = 0.4) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(
        family = "Helvetica",
        face = "bold",
        size = (15),
        hjust = 0.5
      )
    ) +
    ggplot2::scale_colour_viridis_c(option = "magma", direction = -1) +
    ggplot2::coord_fixed() +
    ggplot2::ggtitle(gene)
  
  plot_list[[gene]] <- plot
} #end for loop

# Create grid of expression plots
plot_grid <- cowplot::plot_grid(
  plotlist = plot_list,
  nrow = 5, 
  ncol = 4, 
  align = "hv"
)

# Extract legend
leg <- ggpubr::get_legend(
  plot_list[[1]] + ggplot2::theme(
    legend.position = "right", 
    legend.text = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank()
  )
)

# Add legend to plot_grid
cowplot::plot_grid(
  plot_grid,
  leg,
  nrow = 1,
  ncol = 2,
  rel_widths = c(4,1),
  align = "hv"
)


# Save plots 
ggplot2::ggsave(
  filename = paste0(path_plots_dir, "/FigS3B_FeaturePlots_lymphoid_marker_genes.eps"),
  width    = 14,
  height   = 14
)

# End of document 
################################################