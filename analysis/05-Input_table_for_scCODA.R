################################################
# Script Description ====
# This script produces the input data table for compositional analysis with scCODA
# MAIN STEPS:
# 1- Load annotations of each single-cell from the myeloid and lymphoid analysis
# 2- Produce a table with cell counts of each cell type for each bat individual

# Set Working Directory ====
dir <- "."
setwd(dir)

#-------------------------------------------------------------------------------------------------------------

# Extract annotations and individuals for single cells relevant for scCODA analysis ====
# Load annotations from myeloid and lymphoid cells 
path_seuobj_dir <- "../outputs/SeuratObjects"

extract_cell_annotation <- function(seuobj, path_dir) {
  ds <- readRDS(paste0(path_dir, "/", seuobj))
  labels <- dplyr::tibble(CBC = colnames(ds),
                          Individual = ds$individual,
                          Age = ds$age,
                          Type = ds$type,
                          Annotation = ds$Refined_annotation)
  return(labels)
}
  
myeloid_labels <-extract_cell_annotation("3-ds_annotated_myeloid_dataset.rds", path_seuobj_dir)
lymphoid_labels <-extract_cell_annotation("4-ds_annotated_lymphoid_dataset.rds", path_seuobj_dir)

# Reshape as a single table
data <- rbind(myeloid_labels, lymphoid_labels)

# Combine together the 2 Neutrophil clusters
data$Annotation[grep("^Neutro", data$Annotation)] <- "Neutrophils"

# Remove Doublets, Low quality cells and Outlier Tcell clusters 
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
  data <- subset(data, subset = Annotation != cluster)
} # end for loop

# Keep only cells from Whole Blood samples, as we don't have data across the 3 age classes for the PBMC samples
data <- subset(data, subset = Type == "WholeBlood")

#-------------------------------------------------------------------------------------------------------------

# Reshape data to have total cell counts per individual ====
# Extract total count per cell type per individual 
table <- dplyr::data_frame()

for (ind in unique(data$Individual)) {
  tab <- data[which(data$Individual == ind),] #select a specific indiviudal
  for (ann in unique(tab$Annotation)) {
    count <- length(which(tab$Annotation == ann))
    line <- c(ind, ann, count)
    table <- rbind(table, line)
  } #end for loop ind
} #end for loop ann

# Reshape table and add information about Age
colnames(table) <- c("Individual", "Annotation", "Count")
table <- reshape2::dcast(table, Individual ~ Annotation, value.var = "Count")

look_age <- read.csv("../data/Library_info.csv")
table$Age <- sapply(
  X = table$Individual, 
  FUN = function(x) look_age$age[match(x, look_age$individual)]
)

# Reorder the columns
table <- table[, c(
  "Individual", 
  "Age", 
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
)]


table$Age <- factor(table$Age, levels = c("Adult", "Subadult", "Juvenile")) #order ages properly
table <- dplyr::arrange(table, Age)

# Safety: Replace NA counts (i.e. if cluster absent from one individual) by 0
table[, -c(1,2)] <- apply(
  X = table[,-c(1,2)], 
  MARGIN = c(1,2), 
  FUN = function(x) if(is.na(x)) 0 else (x)
) 

# Save table for scCODA
scCODA_dir <- "../outputs/scCODA"
dir.create(scCODA_dir)

write.csv(table, paste0(scCODA_dir, "/Input_table_scCODA.csv"), row.names = F)

# End of document 
################################################
