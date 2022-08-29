################################################
# Script Description ====
# This script produces the plots used in Supplementary Figure 2 of the manuscript

# Set Working Directory ====
dir <- "."
setwd(dir)

# Set Output Directory ====
path_plots_dir <- "../outputs/Plots/FigS2"
dir.create(path_plots_dir, recursive = T)

#-------------------------------------------------------------------------------------------------------------

# Fig S2B: Leukocyte viability ====
# Load data
path_plots_dir <- "../outputs/Plots/FigS2"
data <- read.csv("../data/Data_FlowCytometry_FigS2.csv")

# Add a GAP row data for plotting between Whole Blood and PBMC samples later on
data[nrow(data)+1,] <- c(rep("GAP", 3), 0, 0)

# Order Flow cytometry data for plotting 
data$Age <- factor(data$Age, levels = c("Adult", "Subadult", "Juvenile", "GAP"))
data$Type <- factor(data$Type, levels = c("WholeBlood", "GAP", "PBMC"))
data <- dplyr::arrange(data, Type, Age)
data$Sample <- factor(data$Sample, levels = unique(data$Sample))

# Generate plot
ann <- data.frame(
  x = c(5, 12),
  y = c(115, 115),
  label = c('Whole Blood', 'PBMCs')
)

ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
     x    = Sample, 
     y    = as.numeric(Viability....)
  )
) +
  cowplot::theme_cowplot() +
  ggplot2::geom_bar(
    stat = "identity", 
    color = "black", 
    width = 0.8,
    mapping = ggplot2::aes(fill = Age)
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "#F94040",
      "#009E73", 
      "#90BFF9", 
      "#FFFFFF"
    )
  ) +
  ggplot2::ylab("Leukocyte viability (%)") +
  ggplot2::xlab("") +
  ggplot2::scale_y_continuous(
    limits = c(0, 120), 
    expand = c(0,0), 
    breaks = c(0,25,50,75,100)
  ) +
  ggplot2::annotate(
    geom = "segment", 
    x = 1, 
    xend = 9, 
    y = 110, 
    yend = 110, 
    size = 1
  ) +
  ggplot2::annotate(
    geom = "segment",
    x = 11,
    xend = 13,
    y = 110,
    yend = 110,
    size = 1
  ) +
  ggplot2::theme(
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(100,100,100,100)
  ) +
  ggplot2::geom_text(data = ann, ggplot2::aes(x, y, label = label))

# Save plot
ggplot2::ggsave(
  filename = paste0(path_plots_dir, "/FigS2B_Viability_plot.eps"),
  width    = 14,
  height   = 14
)

#-------------------------------------------------------------------------------------------------------------

# Fig S2C: Recovery rate ====
# Load data
path_plots_dir <- "../outputs/Plots/FigS2"
data <- read.csv("../data/Data_FlowCytometry_FigS2.csv")

# Add a GAP row data for plotting between Whole Blood and PBMC samples later on
data[nrow(data)+1,] <- c(rep("GAP", 3), 0, 0)

# Order Flow cytometry data for plotting 
data$Age <- factor(data$Age, levels = c("Adult", "Subadult", "Juvenile", "GAP"))
data$Type <- factor(data$Type, levels = c("WholeBlood", "GAP", "PBMC"))
data <- dplyr::arrange(data, Type, Age)
data$Sample <- factor(data$Sample, levels = unique(data$Sample))

# Generate plot
ann <- data.frame(
  x = c(5, 12),
  y = c(115, 115),
  label = c('Whole Blood', 'PBMCs')
)

ggplot2::ggplot(
  data    = data, 
  mapping = ggplot2::aes(
     x    = Sample,
     y    = as.numeric(Recovery.rate....)
  )
) +
  cowplot::theme_cowplot() +
  ggplot2::geom_bar(
    stat = "identity", 
    color = "black", 
    width = 0.8,
    mapping = ggplot2::aes(fill = Age)
  ) +
  ggplot2::scale_fill_manual(
    values = c("#F94040",
               "#009E73", 
               "#90BFF9", 
               "#FFFFFF")
  ) +
  ggplot2::ylab("Recovery rate (%)") +
  ggplot2::xlab("") +
  ggplot2::scale_y_continuous(
    limits = c(0, 120),
    expand = c(0,0),
    breaks = c(0,25,50,75,100)
  ) +
  ggplot2::annotate(
    geom = "segment", 
    x = 1, 
    xend = 9, 
    y = 110, 
    yend = 110,
    size = 1
  ) +
  ggplot2::annotate(
    geom = "segment",
    x = 11,
    xend = 13,
    y = 110,
    yend = 110,
    size = 1
  ) +
  ggplot2::theme(
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(100,100,100,100)
  ) +
  ggplot2::geom_text(data = ann, ggplot2::aes(x, y, label = label))

# Save plot
ggplot2::ggsave(
  filename = paste0(path_plots_dir, "/FigS2C_Recovery_plot.eps"),
  width    = 14,
  height   = 14
)

#-------------------------------------------------------------------------------------------------------------

# Fig S2D: Total cell number ====
# Load data
path_plots_dir <- "../outputs/Plots/FigS2"
path_seuobj_dir <- "../outputs/SeuratObjects"
ds <- readRDS(paste0(path_seuobj_dir, "/2-ds_annotated_full_dataset.rds"))

# Extract data about cell number per individual
data <- as.data.frame(table(ds$individual, ds$type))
colnames(data) <- c("Individual", "Type", "Count")
data <- data[!data$Count==0,] #drop rows with 0 counts (i.e. counts of PBMCs for individual with only Whole Blood data )
data$Individual <- as.character(data$Individual) #to remove factors
data$Type <- as.character(data$Type) #to remove factors
data$Age <- substr(data$Individual, 1, 3)
data[data == "Adu"] <- "Adult"
data[data == "Sub"] <- "Subadult"
data[data == "Juv"] <- "Juvenile"
data$Sample <- paste0(data$Individual, "_", data$Type)

# Add a GAP row data for plotting between Whole Blood and PBMC samples 
data[nrow(data)+1,] <- c("GAP", "GAP", 0, "GAP", "GAP")

# Order cell count data for plotting 
data$Age <- factor(data$Age, levels = c("Adult", "Subadult", "Juvenile", "GAP"))
data$Type <- factor(data$Type, levels = c("WholeBlood", "GAP", "PBMC"))
data <- dplyr::arrange(data, Type, Age)
data$Sample <- factor(data$Sample, levels = unique(data$Sample))

# Generate plot 
ann <- data.frame(
  x = c(5, 12),
  y = c(2500, 2500),
  label = c('Whole Blood', 'PBMCs')
)

ggplot2::ggplot(
  data    = data, 
  mapping = ggplot2::aes(
      x   = Sample, 
      y   = as.numeric(Count)
  )
) +
  cowplot::theme_cowplot() +
  ggplot2::geom_bar(
    stat = "identity", 
    color = "black", 
    width = 0.8,
    mapping = ggplot2::aes(fill = Age)
  ) +
  ggplot2::scale_fill_manual(
    values = c("#F94040",
               "#009E73",
               "#90BFF9",
               "#FFFFFF")
  ) +
  ggplot2::ylab("Total cell number") +
  ggplot2::xlab("") +
  ggplot2::scale_y_continuous(
    limits = c(0, 2550),
    expand = c(0,0),
    breaks = c(0,500,1000,1500,2000,2500)
  ) +
  ggplot2::annotate(
    geom = "segment",
    x = 1,
    xend = 9,
    y = 2400,
    yend = 2400,
    size = 1
  ) +
  ggplot2::annotate(
    geom = "segment",
    x = 11,
    xend = 13,
    y = 2400,
    yend = 2400, 
    size = 1
  ) +
  ggplot2::theme(
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(100,100,100,100)
  ) +
  ggplot2::geom_text(data = ann, ggplot2::aes(x, y, label = label))

# Save plot
ggplot2::ggsave(
  filename = paste0(path_plots_dir, "/FigS2D_Total_cell_number_plot.eps"),
  width    = 14,
  height   = 14
)

#-------------------------------------------------------------------------------------------------------------

# Fig S2E: Doublet UMAPs ====
# Load data
path_plots_dir <- "../outputs/Plots/FigS2"
path_seuobj_dir <- "../outputs/SeuratObjects"

ds.full <- readRDS(paste0(path_seuobj_dir, "/2-ds_annotated_full_dataset.rds"))
ds.myelo <- readRDS(paste0(path_seuobj_dir, "/3-ds_annotated_myeloid_dataset.rds"))
ds.lymph <- readRDS(paste0(path_seuobj_dir, "/4-ds_annotated_lymphoid_dataset.rds"))

#Function to create UMAPs with doublets 
create_umap_doublet <- function(object, title) {
  data <- dplyr::as_tibble(
    object@reductions$umap@cell.embeddings
  )
  colnames(data) <- c("x", "y")
  
  data$status <- object$status
  data$status <- factor(data$status, levels = c("singlet", "doublet")) #order factors to get doublets plotted on the top of singlets
  
  arrow <- grid::arrow(
    angle  = 30, 
    length = ggplot2::unit(0.25, "inches"), 
    type   = "closed"
  )
  
ggplot2::ggplot(
  data    = dplyr::arrange(data, status),
  mapping = ggplot2::aes(
    x     = x,
    y     = y,
    col   = status
  )
) +
    ggplot2::ggtitle(label = title) +
    ggplot2::geom_point(size  = 0.5) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::scale_color_manual(values = c("#e6e6e6", "#000000")) +
    ggplot2::coord_fixed()
}

# Create the grid with UMAPs
grid_UMAP <- cowplot::plot_grid(
  create_umap_doublet(ds.full, "Blood leukocytes") + ggplot2::theme(legend.position = "none"),
  create_umap_doublet(ds.myelo, "Myeloid") + ggplot2::theme(legend.position = "none"),
  create_umap_doublet(ds.lymph, "Lymphoid") + ggplot2::theme(legend.position = "none"),
  nrow = 1, 
  ncol = 3, 
  align = "hv", 
  scale = 0.8
)

#Combine plots with legend and save
cowplot::plot_grid(
  grid_UMAP,
  cowplot::get_legend(create_umap_doublet(ds.full, "Blood leukocytes")),
  nrow = 2,
  ncol = 1,
  align = "hv",
  rel_heights = c(1, 0.2)
)

ggplot2::ggsave(
  filename = paste0(path_plots_dir, "/FigS2E_grid_UMAP_doublets.eps"),
  width    = 14,
  height   = 14
)

#-------------------------------------------------------------------------------------------------------------

# Fig S2F: Comparison Whole Blood vs PBMC samples ====
# Load  and prep data 
path_plots_dir <- "../outputs/Plots/FigS2"
path_seuobj_dir <- "../outputs/SeuratObjects"

ds <- readRDS(paste0(path_seuobj_dir, "/2-ds_annotated_full_dataset.rds"))

# Load annotations from myeloid and lymphoid cells 
extract_cell_annotation <- function(seuobj, path_dir) {
  ds <- readRDS(paste0(path_dir, "/", seuobj))
  labels <- dplyr::tibble(CBC = colnames(ds),
                          Annotation = ds$Refined_annotation)
  return(labels)
}

myeloid_labels <-extract_cell_annotation("3-ds_annotated_myeloid_dataset.rds", path_seuobj_dir)
lymphoid_labels <-extract_cell_annotation("4-ds_annotated_lymphoid_dataset.rds", path_seuobj_dir)

# Reshape as a single table
lookup <- rbind(myeloid_labels, lymphoid_labels)

# Combine together the 2 Neutrophil clusters
lookup$Annotation[grep("^Neutro", lookup$Annotation)] <- "Neutrophils"

# Add refined annotation for each single cell and exclude non-relevant cells
ds$Refined_annotation <- sapply(ds$CBC, function(x) lookup$Annotation[match(x, lookup$CBC)])

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

# Keep only the individuals with WholeBlood and PBMCs data
ds <- subset(ds, subset = individual == "Adu4" | individual == "Sub1" | individual == "Sub2")

# Extract data for plotting 
data <- dplyr::as_tibble(
  ds@meta.data
)

data <- dplyr::select(data, c(individual, Refined_annotation, type))

data$Refined_annotation <- factor(
  ds$Refined_annotation, 
  levels = c(
    "Neutrophils",
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
    "Cycling Lymphocytes", 
    "GAP"
  )
)

data$sample <- paste0(data$individual, "_", data$type)

# Order samples appropriately
data$sample <- factor(
  data$sample, 
  levels = c(
    "Adu4_WholeBlood", 
    "Sub1_WholeBlood", 
    "Sub2_WholeBlood", 
    "Adu4_PBMC", 
    "Sub1_PBMC", 
    "Sub2_PBMC"
  )
)

data$type <- factor(data$type, levels = c("WholeBlood", "PBMC"))

# Prepare Color palette
colors <- c(
  "#E60000", #Neutrophils
  "#B30000", #Eosinophils 
  "#CC6600", #CD16 Monocytes
  "#FFBF80", #CD14 Monocytes
  "#FF8000", #cDC2
  "#E6E600", #cDC1
  "#E60073", #pDC
  "#CC9966", #NKT-like 2
  "#86592D", #NKT-like 1
  "#888844", #CD8 EM
  "#BBBB77", #CD8 CM
  "#006622", #CD4 EM
  "#009933", #CD4 CM1
  "#53C68C", #CD4 CM2
  "#00E600", #Treg
  "#B3E6CC", #DN T cells
  "#66FF99", #Naive T cells
  "#002B80", #ZBTB32/B cells
  "#004DE6", #PLAC8/B cells
  "#4D88FF", #LTB/B cells
  "#99BBFF", #VPREB3/B cells
  "#660066", #Cycling Lymphocytes
  "#FFFFFF" # GAP
)

color_palette <- dplyr::tibble(
  Refined_annotation = levels(data$Refined_annotation),
  color = colors
)

data$colors <- sapply(
  X = data$Refined_annotation, 
  FUN = function(x) color_palette$color[match(x, color_palette$Refined_annotation)]
)

# Create Barplot 
ggplot2::ggplot(
  data    = data, 
  mapping = ggplot2::aes(
     x    = individual,
     fill = Refined_annotation
  )
) +
  cowplot::theme_cowplot() +
  ggplot2::geom_bar(position = "fill") +
  ggplot2::scale_fill_manual(values = colors) +
  ggplot2::facet_wrap(~type) +
  ggplot2::ylab("Fraction") +
  ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      vjust = 1)
  ) 

# Save plot
ggplot2::ggsave(
  filename = paste0(path_plots_dir, "/FigS2F_Comparison_WholeBlood_PBMC.eps"),
  width    = 14,
  height   = 14
)

#-------------------------------------------------------------------------------------------------------------

# Fig S2G: Number of gene detected per cell type ====

# Load data
path_plots_dir <- "../outputs/Plots/FigS2"
path_seuobj_dir <- "../outputs/SeuratObjects"

ds <- readRDS(paste0(path_seuobj_dir, "/2-ds_annotated_full_dataset.rds"))

# Load annotations from myeloid and lymphoid cells 
extract_cell_annotation <- function(seuobj, path_dir) {
  ds <- readRDS(paste0(path_dir, "/", seuobj))
  labels <- dplyr::tibble(CBC = colnames(ds),
                          Annotation = ds$Refined_annotation)
  return(labels)
}

myeloid_labels <-extract_cell_annotation("3-ds_annotated_myeloid_dataset.rds", path_seuobj_dir)
lymphoid_labels <-extract_cell_annotation("4-ds_annotated_lymphoid_dataset.rds", path_seuobj_dir)

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
}

# Order clusters 
ds$Refined_annotation <- factor(
  ds$Refined_annotation, 
  levels = c(
    "Neutrophils", 
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
    "Cycling Lymphocytes"
  )
)

# Prepare Color palette
colors <- c(
  "#E60000", # Neutrophils 1
  "#FF4D4D", # Neutrophils 2
  "#B30000", # Eosinophils 
  "#CC6600", # CD16 Monocytes
  "#FFBF80", # CD14 Monocytes
  "#FF8000", # cDC2
  "#E6E600", # cDC1
  "#E60073", # pDC
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
  "#660066" # Cycling Lymphocytes
) 

# Generate and save plot 
Seurat::VlnPlot(
  ds, 
  features = "nFeature_RNA", 
  group.by = "Refined_annotation", 
  cols = colors,
  pt.size = 0
) +
  ggplot2::geom_violin(draw_quantiles = 0.5,  scale = "width") +
  ggplot2::ylab("Number of genes detected per cell") +
  ggplot2::xlab("") +
  ggplot2::scale_y_continuous(
    breaks=c(0, 1000, 2000, 3000, 4000, 5000, 6000, 7000), 
    limits = c(0,7000)
    ) +
  ggplot2::theme(legend.position="none") +
  ggplot2::ggtitle("")

ggplot2::ggsave(
  filename = paste0(path_plots_dir, "/FigS2G_number_genes_celltype.eps"),
  width    = 10,
  height   = 8
)

#-------------------------------------------------------------------------------------------------------------

# Fig S2H: Host genes related to Marburg virus ====

# Load data
path_plots_dir <- "../outputs/Plots/FigS2"
path_seuobj_dir <- "../outputs/SeuratObjects"

ds <- readRDS(paste0(path_seuobj_dir, "/2-ds_annotated_full_dataset.rds"))

# Load annotations from myeloid and lymphoid cells 
extract_cell_annotation <- function(seuobj, path_dir) {
  ds <- readRDS(paste0(path_dir, "/", seuobj))
  labels <- dplyr::tibble(CBC = colnames(ds),
                          Annotation = ds$Refined_annotation)
  return(labels)
}

myeloid_labels <-extract_cell_annotation("3-ds_annotated_myeloid_dataset.rds", path_seuobj_dir)
lymphoid_labels <-extract_cell_annotation("4-ds_annotated_lymphoid_dataset.rds", path_seuobj_dir)

# Reshape as a single table
lookup <- rbind(myeloid_labels, lymphoid_labels)

# Combine together the 2 Neutrophil clusters
lookup$Annotation[grep("^Neutro", lookup$Annotation)] <- "Neutrophils"

# Add refined annotation for each single cell and exclude non-relevant cells
ds$Refined_annotation <- sapply(ds$CBC, function(x) lookup$Annotation[match(x, lookup$CBC)])

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
}

# Selected marker genes
genes <- c(
  "CLEC4G", "HAVCR1", "TIMD4", "ITGB1",
  "TYRO3", "CTSL", "CTSB",
  "NPC1", "TPCN2", "TOP1",
  "DYNLL1", "SEC61A1", "HSP90AA1",
  "GSPT1", "UPF1", "IQGAP1",
  "SEC24C", "TSG101", "NEDD4",
  "SOCS3", "TLR4", "BST2"
)

# Extract and prepare for Dotplot 
ids <- rownames(ds@assays$RNA@data) %in% genes
data <- t(as.matrix(ds@assays$RNA@data[ids, ]))
data <- tidyr::as_tibble(data)
data <- dplyr::select(data, genes) #reorder columns
data$Refined_annotation <- ds@meta.data$Refined_annotation

data$Refined_annotation <- factor(
  data$Refined_annotation, 
  levels = c(
    "Neutrophils",
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
    "Cycling Lymphocytes"
  )
)

data <- tidyr::gather(data, "Gene", "Expression", -Refined_annotation)
data$Gene <- factor(
  x      = data$Gene,
  levels = unique(data$Gene)
)

data$Pct <- data$Expression > 0
data$Freq <- rep(1, length(data$Pct))

data <- dplyr::group_by(data, Gene)
data <- dplyr::mutate(data, Scaled = scale(Expression)[, 1])

data <- dplyr::group_by(data, Refined_annotation, Gene)
data <- dplyr::summarise(
  data,
  Mean   = mean(Expression),
  Scaled = mean(Scaled),
  Pct    = sum(Pct)/sum(Freq)*100
)
data$Scaled[data$Scaled > 1] <- 1
data$Mean[data$Mean > 2] <- 2

# Generate Dotplot 
colorscale <- RColorBrewer::brewer.pal(n = 9, name = "RdBu")

ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    y     = Refined_annotation,
    x     = Gene,
    size  = Pct,
    col   = Scaled
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_color_gradient2(
    midpoint=0, 
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
      barheight = 11, 
      frame.colour = "black", 
      ticks = FALSE, 
      order = 1
    ),
    size  = ggplot2::guide_legend(
      label.position = "right"
    )
  ) +
  ggplot2::scale_size_area(max_size = 8)

# Save Dotplot 
ggplot2::ggsave(
  filename = paste0(path_plots_dir, "/FigS2H_Host_genes_related_MARV.eps"),
  width    = 16,
  height   = 7
)

# End of document 
################################################
