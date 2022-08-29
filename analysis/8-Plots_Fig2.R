################################################
# Script Description ====
# This script produces the plots used in Figure 2 of the manuscript

# Set Working Directory ====
dir <- "."
setwd(dir)

# Set Output Directory ====
path_plots_dir <- "../outputs/Plots/Fig2"
dir.create(path_plots_dir, recursive = T)

#-------------------------------------------------------------------------------------------------------------

# Fig 2A: UMAPs highlighthing Adult, Subadult and Juvenile cells ====

# Load data
path_plots_dir <- "../outputs/Plots/Fig2"
path_seuobj_dir <- "../outputs/SeuratObjects"

ds_myelo <- readRDS(paste0(path_seuobj_dir, "/3-ds_annotated_myeloid_dataset.rds"))
ds_lymph <- readRDS(paste0(path_seuobj_dir, "/4-ds_annotated_lymphoid_dataset.rds"))

# Function to create UMAPs and highlight age
create_umap_age <- function(object, ordered_levels, color) {
  data <- dplyr::as_tibble(
    object@reductions$umap@cell.embeddings
  )
  colnames(data) <- c("x", "y")
  
  data$Age <- object$age
  data$Age <- factor(data$Age, levels = ordered_levels) #important, to have cells of age of interest plotted on top (last)
  
  arrow <- grid::arrow(
    angle  = 30, 
    length = ggplot2::unit(0.25, "inches"), 
    type   = "closed"
  )
  
  
  ggplot2::ggplot(
    data    = dplyr::arrange(data, Age),
    mapping = ggplot2::aes(
      x   = x,
      y   = y,
      col = Age
    )
  ) +
    ggplot2::geom_point(size  = 0.5) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "") +
    ggplot2::scale_color_manual(values = c("#E6E6E6", "#E6E6E6", color)) +
    ggplot2::coord_fixed()

}

# Create UMAPS 
plot_list <- list()
plot_list[[1]] <- create_umap_age(ds_myelo, c("Juvenile", "Subadult", "Adult"), "#F94040")
plot_list[[2]] <- create_umap_age(ds_myelo, c("Juvenile", "Adult", "Subadult"), "#009E73")
plot_list[[3]] <- create_umap_age(ds_myelo, c("Subadult", "Adult", "Juvenile"), "#90BFF9")

plot_list[[4]] <- create_umap_age(ds_lymph, c("Juvenile", "Subadult", "Adult"), "#F94040")
plot_list[[5]] <- create_umap_age(ds_lymph, c("Juvenile", "Adult", "Subadult"), "#009E73")
plot_list[[6]] <- create_umap_age(ds_lymph, c("Subadult", "Adult", "Juvenile"), "#90BFF9")

# Save UMAPs as a grid 
cowplot::plot_grid(
  plotlist = plot_list,
  nrow = 2, 
  ncol = 3, 
  align = "hv",
  labels = c("Adult", "Subadult", 'Juvenile'),
  label_x = c(0.3, 0.2, 0.25)
) +
  cowplot::draw_label(
    "Myeloid", 
    x = 0, 
    y =0.75, 
    vjust = 1.5, 
    angle =90, 
    fontface = "bold"
    
  ) +
  cowplot::draw_label(
    "Lymphoid", 
    x = 0, 
    y =0.25, 
    vjust= 1.5, 
    angle =90, 
    fontface = "bold"
  )


ggplot2::ggsave(
  filename = paste0(path_plots_dir, "/Fig2A_grid_UMAP_age.eps"),
  width    = 14,
  height   = 14
)

#-------------------------------------------------------------------------------------------------------------

# Fig 2B: Boxplot from scCODA results ====

# Load data and merge all results as a single table
path_plots_dir <- "../outputs/Plots/Fig2"
path_scCODA_dir <- "../outputs/scCODA"

fdr_val <- c(0.1, 0.05)
fdr_val <- fdr_val[order(fdr_val, decreasing = F)] # important to order correctly because of the filtering step for lowest FDRs

merged_table <- list()
for (fdr in fdr_val) {
  list_files <- list.files(path_scCODA_dir, pattern = as.character(fdr))
  
  list_tables <- list()
  for (file in list_files) {
    list_tables[[file]] <- read.csv(file = paste0(path_scCODA_dir,"/",file))
  } #end for loop file
  table <- data.table::rbindlist(list_tables)
  table$FDR <- as.character(fdr)
  merged_table[[as.character(fdr)]] <- table
} #end for loop fdr

merged_table <- data.table::rbindlist(merged_table)

# Function to extract only significant results from scCODA
extract_scCODA_results <- function(merged_table, fdr_val) {
  table_list <- list()
  for (fdr in fdr_val) {
    data <- merged_table[merged_table$FDR == fdr, ]
    data$Ref_model <- stringi::stri_extract_all_regex(data$Covariate, "(?<=').*?(?=')") #extract the age used as reference for scCODA
    data$Comparison <-  stringi::stri_extract_all_regex(data$Covariate, "(?<=\\[T.)[^{}]+(?=\\])") #extract the age compared to the scCODA reference
    data$ID_Comp <- paste0(data$Cell.Type, "-", data$Ref_model, "-", data$Comparison)
    data$ID_Comp_reverse <- paste0(data$Cell.Type, "-", data$Comparison, "-", data$Ref_model)
    data$Final.Parameter.reverse <- data$Final.Parameter[match(data$ID_Comp_reverse, data$ID_Comp)]
    
    data <- data[data$Final.Parameter == "True" & data$Final.Parameter.reverse == "True", ] # difference must be significant independently of the age used as reference
    
    # Clean up table to remove significant results reported twice
    maxi <- nrow(data)
    
    for (i in data$ID_Comp) {
      counter <- 1
      if (counter > maxi) {
        break
      } else if (data$ID_Comp_reverse[counter] %in% data$ID_Comp) {
        idx <- grep(data$ID_Comp_reverse[counter], data$ID_Comp_reverse)
        data <- data[-idx, ]
      } else {
        counter <- counter + 1
      }
    } #end i for loop
    table_list[[as.character(fdr)]] <- data
  } # end fdr for loop
  data.table::rbindlist(table_list)
}

# Call function and keep only lowest possible FDRs
scCODA_results <- extract_scCODA_results(merged_table, fdr_val)
scCODA_results <- scCODA_results[-which(duplicated(scCODA_results$ID_Comp)), ]

print("Here are the significant differences based on scCODA: ")
print(dplyr::select(scCODA_results, ID_Comp, FDR))

# Load and reshape data for plotting as a boxplot
data <- dplyr::as_tibble(read.csv("../outputs/scCODA/Input_table_scCODA.csv"))
colnames(data) <- gsub("[.]+", " ", colnames(data)) #cleanup colnames
data[, -c(1,2)] <- prop.table(as.matrix(data[,-c(1,2)]), margin = 1) * 100 #convert counts to percentages

annotations <- colnames(data)[-c(1,2)]
data <- tidyr::pivot_longer(data, cols = annotations, names_to = "Annotation") # reshape table in long format
data$Annotation <- factor(data$Annotation, levels = annotations)
data$Age <- factor(data$Age, levels = c("Adult", "Subadult", "Juvenile"))
  
# Create boxplot
# Significant bars are based on scCODA_results
print(dplyr::select(scCODA_results, ID_Comp, FDR))

ggplot2::ggplot(
  data    = data,
  mapping = ggplot2::aes(
    x     = Annotation,
    y     = value,
    fill  = Age
  )
) +
  cowplot::theme_cowplot() +
  ggplot2::geom_boxplot(coef=10) + # increase coef (default =1.5) so that outlier values are included in the whiskers
  ggplot2::scale_fill_manual(values = c("#F94040", "#009E73", "#90BFF9")) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      vjust = 1) 
  ) +
  ggplot2::xlab("") +
  ggplot2::ylab("Percent") +
  ggplot2::scale_y_continuous(
    expand = c(0, 0), 
    limits = c(0,65), 
    breaks = c(0, 10, 20, 30, 40, 50, 60)
  ) +
  ggpubr::geom_signif(
    y_position = c(62,59,27,27,27,23,23), 
    xmin = c(0.75,0.75,18.75,19.75,20.75,18.75,19.75), 
    xmax = c(1.25,1.00,19.25,20.25,21.25,19,20),
    annotation=c("**","*","**","**","*","*","**"), 
    size = 1, 
    tip_length=0, 
    col = c("black")
  )

# Save Boxplot
ggplot2::ggsave(
  filename = paste0(path_plots_dir, "/Fig2B_scCODA_boxplot.eps"),
  width    = 14,
  height   = 7
)

# End of document 
################################################