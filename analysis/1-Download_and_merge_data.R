################################################
# Script Description ====
# This script downloads the count matrices from the GEO submission and prepare the data as a single dataset
# for further analysis
# MAIN STEPS:
# 1- Download count matrices from GEO
# 2- Add metadata obtained from Souporcell (demultiplexing of the cell pools)
# 3- Merge all matrices together to get the full dataset as a single object

#Set Working Directory ====
dir <- "."
setwd(dir)

#-------------------------------------------------------------------------------------------------------------
# Download GEO data and prepare files for each library  ==== 
# As GEO submission is still not publicly released, count matrices are provided in the data folder: GSE183925_RAW.tar

# Download the data
# to be used when publicly released on GEO
# url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE183925&format=file" 
url <- "https://nubes.helmholtz-berlin.de/s/PmK3zz8eLntfGKL/download/GSE183925_RAW.tar"
path_geo_dir <- "../data/GEO_download"
geo_file <- paste0(path_geo_dir, "/GSE183925_RAW.tar")

dir.create(path_geo_dir)
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
  file.rename(paste0(path_geo_dir, "/library_", library, "/", old_filenames),
              paste0(path_geo_dir, "/library_", library, "/", new_filenames))
} 

#-------------------------------------------------------------------------------------------------------------
# Create SeuratObjects with metadata from downloaded GEO data ====
data.A1 <- Seurat::Read10X(data.dir = paste0(path_geo_dir, "/library_A1"))
data.A2 <- Seurat::Read10X(data.dir = paste0(path_geo_dir, "/library_A2"))
data.A5 <- Seurat::Read10X(data.dir = paste0(path_geo_dir, "/library_A5"))
data.A6 <- Seurat::Read10X(data.dir = paste0(path_geo_dir, "/library_A6"))

ds.A1 <- Seurat::CreateSeuratObject(counts = data.A1, project = "ImmuBat", min.cells = 3, min.features = 200)
ds.A2 <- Seurat::CreateSeuratObject(counts = data.A2, project = "ImmuBat", min.cells = 3, min.features = 200)
ds.A5 <- Seurat::CreateSeuratObject(counts = data.A5, project = "ImmuBat", min.cells = 3, min.features = 200)
ds.A6 <- Seurat::CreateSeuratObject(counts = data.A6, project = "ImmuBat", min.cells = 3, min.features = 200)

rm(data.A1, data.A2, data.A5, data.A6)

# Add Souporcell results as metadata
path_soup_dir <- "../data/Souporcell"

# Load Souporcell results per library
soup.A1 <- read.delim(paste0(path_soup_dir, "/clusters_libA1.tsv"))
soup.A2 <- read.delim(paste0(path_soup_dir, "/clusters_libA2.tsv"))
soup.A5 <- read.delim(paste0(path_soup_dir, "/clusters_libA5.tsv"))
soup.A6 <- read.delim(paste0(path_soup_dir, "/clusters_libA6.tsv"))

# Add assignment and singlet/doublet data to the SeuratObjects 
ds.A1@meta.data$soup_cluster <- sapply(rownames(ds.A1@meta.data), 
                                       function(x) soup.A1$assignment[match(x, soup.A1$barcode)])
ds.A2@meta.data$soup_cluster <- sapply(rownames(ds.A2@meta.data), 
                                       function(x) soup.A2$assignment[match(x, soup.A2$barcode)])
ds.A5@meta.data$soup_cluster <- sapply(rownames(ds.A5@meta.data), 
                                       function(x) soup.A5$assignment[match(x, soup.A5$barcode)])
ds.A6@meta.data$soup_cluster <- sapply(rownames(ds.A6@meta.data), 
                                       function(x) soup.A6$assignment[match(x, soup.A6$barcode)])

ds.A1@meta.data$status <- sapply(rownames(ds.A1@meta.data), 
                                 function(x) soup.A1$status[match(x, soup.A1$barcode)])
ds.A2@meta.data$status <- sapply(rownames(ds.A2@meta.data), 
                                 function(x) soup.A2$status[match(x, soup.A2$barcode)])
ds.A5@meta.data$status <- sapply(rownames(ds.A5@meta.data), 
                                 function(x) soup.A5$status[match(x, soup.A5$barcode)])
ds.A6@meta.data$status <- sapply(rownames(ds.A6@meta.data), 
                                 function(x) soup.A6$status[match(x, soup.A6$barcode)])

# Remove cells reported as unassigned by souporcell
ds.A1 <- subset(ds.A1, subset = status != "unassigned")
ds.A2 <- subset(ds.A2, subset = status != "unassigned")
ds.A5 <- subset(ds.A5, subset = status != "unassigned")
ds.A6 <- subset(ds.A6, subset = status != "unassigned")

# Add library name as metadata
ds.A1$library <- "A1"
ds.A2$library <- "A2"
ds.A5$library <- "A5"
ds.A6$library <- "A6"
#-------------------------------------------------------------------------------------------------------------

# Merge all datasets and save as a single Seurat Object ====
#path_outputs_dir <- "../outputs"
#dir.create(path_outputs_dir)

path_seuobj_dir <- "../outputs/SeuratObjects"
dir.create(path_seuobj_dir, recursive = T)

ds.merged <- merge(ds.A1, 
                     y = c(ds.A2, ds.A5, ds.A6), 
                     add.cell.ids = c("A1", "A2", "A5", "A6"), 
                     project = "ImmuBat")
saveRDS(ds.merged, paste0(path_seuobj_dir, "/1-ds_merged_libraries.rds"))

# End of document 
################################################