################################################
# Script Description ====
# This script downloads the count matrices from the GEO submission and prepare the data as a single dataset
# for further analysis
# MAIN STEPS:
# 1- Download count matrices from GEO
# 2- Add metadata obtained from Souporcell (demultiplexing of the cell pools)
# 3- Merge all matrices together to get the full dataset as a single object

# Set Working Directory ====
dir <- "."
setwd(dir)

#-------------------------------------------------------------------------------------------------------------

# Download GEO data and prepare files for each library  ==== 
# Download the data
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE183925&format=file" 
path_geo_dir <- "../data/GEO_download"
geo_file <- paste0(path_geo_dir, "/GSE183925_RAW.tar")
dir.create(path_geo_dir)

options(timeout=600) # To avoid being timeout during the download. Increase if necessary.
download.file(url, geo_file)

# Untar and extract the data from each library in dedicated folders
libraries <- c("A1", "A2", "A5", "A6")
list_tar_files <- untar(tarfile = geo_file,
                        list = T)

for (library in libraries) {
  untar(tarfile = geo_file,
        files = list_tar_files[grep(library, list_tar_files)],
        exdir = paste0(path_geo_dir, "/library_", library))
  
  old_filenames <- list.files(paste0(path_geo_dir, "/library_", library))
  new_filenames <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
  
  file.rename(
    paste0(path_geo_dir, "/library_", library, "/", old_filenames),
    paste0(path_geo_dir, "/library_", library, "/", new_filenames)
  )
} # end for loop 

#-------------------------------------------------------------------------------------------------------------

# Create SeuratObjects with metadata from downloaded GEO data ====
data_A1 <- Seurat::Read10X(data.dir = paste0(path_geo_dir, "/library_A1"))
data_A2 <- Seurat::Read10X(data.dir = paste0(path_geo_dir, "/library_A2"))
data_A5 <- Seurat::Read10X(data.dir = paste0(path_geo_dir, "/library_A5"))
data_A6 <- Seurat::Read10X(data.dir = paste0(path_geo_dir, "/library_A6"))

ds_A1 <- Seurat::CreateSeuratObject(counts = data_A1, project = "ImmuBat", min.cells = 3, min.features = 200)
ds_A2 <- Seurat::CreateSeuratObject(counts = data_A2, project = "ImmuBat", min.cells = 3, min.features = 200)
ds_A5 <- Seurat::CreateSeuratObject(counts = data_A5, project = "ImmuBat", min.cells = 3, min.features = 200)
ds_A6 <- Seurat::CreateSeuratObject(counts = data_A6, project = "ImmuBat", min.cells = 3, min.features = 200)

rm(data_A1, data_A2, data_A5, data_A6)

# Load Souporcell results per library
path_soup_dir <- "../data/Souporcell"

soup_A1 <- read.delim(paste0(path_soup_dir, "/clusters_libA1.tsv"))
soup_A2 <- read.delim(paste0(path_soup_dir, "/clusters_libA2.tsv"))
soup_A5 <- read.delim(paste0(path_soup_dir, "/clusters_libA5.tsv"))
soup_A6 <- read.delim(paste0(path_soup_dir, "/clusters_libA6.tsv"))

# Add assignment and singlet/doublet data to the SeuratObjects 
ds_A1@meta.data$soup_cluster <- sapply(
  X = rownames(ds_A1@meta.data), 
  FUN = function(x) soup_A1$assignment[match(x, soup_A1$barcode)]
)
ds_A2@meta.data$soup_cluster <- sapply(
  X = rownames(ds_A2@meta.data), 
  FUN = function(x) soup_A2$assignment[match(x, soup_A2$barcode)]
)
ds_A5@meta.data$soup_cluster <- sapply(
  X = rownames(ds_A5@meta.data), 
  FUN = function(x) soup_A5$assignment[match(x, soup_A5$barcode)]
)
ds_A6@meta.data$soup_cluster <- sapply(
  X = rownames(ds_A6@meta.data), 
  FUN = function(x) soup_A6$assignment[match(x, soup_A6$barcode)]
)

ds_A1@meta.data$status <- sapply(
  X = rownames(ds_A1@meta.data), 
  FUN = function(x) soup_A1$status[match(x, soup_A1$barcode)]
)
ds_A2@meta.data$status <- sapply(
  X = rownames(ds_A2@meta.data), 
  FUN = function(x) soup_A2$status[match(x, soup_A2$barcode)]
)
ds_A5@meta.data$status <- sapply(
  X = rownames(ds_A5@meta.data), 
  FUN = function(x) soup_A5$status[match(x, soup_A5$barcode)]
)
ds_A6@meta.data$status <- sapply(
  X = rownames(ds_A6@meta.data), 
  FUN = function(x) soup_A6$status[match(x, soup_A6$barcode)]
)

# Remove cells reported as unassigned by souporcell
ds_A1 <- subset(ds_A1, subset = status != "unassigned")
ds_A2 <- subset(ds_A2, subset = status != "unassigned")
ds_A5 <- subset(ds_A5, subset = status != "unassigned")
ds_A6 <- subset(ds_A6, subset = status != "unassigned")

# Add library name as metadata
ds_A1$library <- "A1"
ds_A2$library <- "A2"
ds_A5$library <- "A5"
ds_A6$library <- "A6"

#-------------------------------------------------------------------------------------------------------------

# Merge all datasets and save as a single Seurat Object ====
path_seuobj_dir <- "../outputs/SeuratObjects"
dir.create(path_seuobj_dir, recursive = T)

ds_merged <- merge(
  x = ds_A1, 
  y = c(ds_A2, ds_A5, ds_A6), 
  add.cell.ids = c("A1", "A2", "A5", "A6"), 
  project = "ImmuBat"
)

saveRDS(
  object = ds_merged,
  file = paste0(path_seuobj_dir, "/1-ds_merged_libraries.rds")
)

# End of document 
################################################